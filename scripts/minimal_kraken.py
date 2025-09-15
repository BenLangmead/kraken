#!/usr/bin/env python3

"""
"Minimal Kraken": a simple Python script recreating Kraken's
functionality to the point of generating identical output (on
everything we've tested so far).

The databse does _not_ use the same format as Kraken, indexes cannot be
swapped between this version and the official version.  This script is
also somewhat more flexible in how it takes inputs; they do not need to
be pre-arranged in a directory structure in the same was as they do for
the original kraken.

Copyright Ben Langmead <blangme2@jhu.edu>, 2025

This file is part of the Kraken taxonomic sequence classification
system, by Derrick Wood and Steven Salzberg.

Kraken is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Kraken is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Kraken.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import argparse
import pickle
import os
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

DEF_KMER_LEN = 31

def print_progress_bar(current: int, total: int, width: int = 50, prefix: str = "Progress"):
    """Print a simple progress bar"""
    percent = 100 if total == 0 else min(100, (current * 100) // total)
    filled = (current * width) // total if total > 0 else width
    bar = '█' * filled + '░' * (width - filled)
    
    # Progress bar w/ carriage return overwrite
    print(f'\r{prefix}: |{bar}| {percent}% ({current}/{total})', end='', flush=True, file=sys.stderr)
    if current >= total:
        print(file=sys.stderr)

class KmerScanner:
    """Extract k-mers from DNA sequences, similar to C++ KmerScanner"""
    __slots__ = ['seq', 'k', 'start', 'finish', 'curr_pos', 'kmer_mask', 'nuc_encoding']
    
    def __init__(self, seq: str, k: int, start: int = 0, finish: Optional[int] = None):
        self.seq = seq.upper()
        self.k = k
        self.start = start
        self.finish = finish if finish is not None else len(seq)
        self.curr_pos = start
        self.kmer_mask = (1 << (k * 2)) - 1
        # Pre-compute nucleotide encoding for faster lookup
        ne = self.nuc_encoding = [0] * 256  # ASCII lookup table
        ne[ord('A')], ne[ord('C')], ne[ord('G')], ne[ord('T')] = 0, 1, 2, 3
        
    def next_kmer(self) -> Optional[Tuple[int, bool]]:
        """Return (kmer_encoding, has_ambiguous) or None if exhausted"""
        if self.curr_pos + self.k > self.finish:
            return None
        
        kmer, has_ambiguous = 0, False
        seq_bytes = self.seq.encode('ascii')  # Convert to bytes for faster access
        for i in range(self.k):
            char_code = seq_bytes[self.curr_pos + i]
            kmer <<= 2
            nuc_val = self.nuc_encoding[char_code]
            if nuc_val == 0 and char_code not in (ord('A'), ord('C'), ord('G'), ord('T')):
                has_ambiguous = True
                continue
            kmer |= nuc_val
                
        self.curr_pos += 1
        return (kmer & self.kmer_mask, has_ambiguous)

def reverse_complement(kmer: int, k: int) -> int:
    """Compute reverse complement of k-mer encoding"""
    rev_comp = 0
    for _ in range(k):
        nt = kmer & 3
        rev_comp = (rev_comp << 2) | (3 - nt)
        kmer >>= 2
    return rev_comp

def canonical_representation(kmer: int, k: int) -> int:
    """Return lexicographically smaller of kmer and its revcomp"""
    # Pre-compute reverse complement only if needed
    rev_comp = reverse_complement(kmer, k)
    return kmer if kmer < rev_comp else rev_comp

def kmer_to_string(kmer: int, k: int) -> str:
    """Convert encoded k-mer back to string"""
    result = []
    for i in range(k):
        result.append('ACGT'[kmer & 3])
        kmer >>= 2
    return ''.join(reversed(result))

class TaxonomyHandler:
    """Handle taxonomy tree operations"""
    
    def __init__(self, nodes_file: Optional[str] = None, names_file: Optional[str] = None):
        self.parent_map: Dict[int, int] = {}
        self.names_map: Dict[int, str] = {}
        self.root_id: int = 1  # Default, will be updated by load_nodes
        if nodes_file:
            self.load_nodes(nodes_file)
        if names_file:
            self.load_names(names_file)
    
    def load_nodes(self, filename: str):
        """Load taxonomy tree from nodes.dmp file"""
        roots = []
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    break
                parts = line.split('\t|\t')
                if len(parts) >= 2:
                    node_id = int(parts[0])
                    parent_id = int(parts[1])
                    if node_id == parent_id:  # root has self as parent
                        roots.append(node_id)                    
                    self.parent_map[node_id] = parent_id
        
        if len(roots) == 0:
            raise ValueError(f"No root node found in {filename}. Root node should be its own parent.")
        elif len(roots) > 1:
            raise ValueError(f"Multiple root nodes (having self as parent) found in {filename}: {roots}. There should be exactly one root.")
        
        self.root_id = roots[0]
        print(f"Taxonomy root node: {self.root_id}", file=sys.stderr)
    
    def load_names(self, filename: str):
        """Load taxon names from names.dmp file"""
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    break
                parts = line.split('\t|\t')
                if len(parts) >= 3:
                    taxon_id = int(parts[0])
                    name, name_type = parts[1], parts[2]
                    # Only store scientific names (type "scientific name")
                    if name_type == "scientific name":
                        self.names_map[taxon_id] = name
    
    def get_taxon_name(self, taxon_id: int) -> str:
        """Get the scientific name for a taxon ID"""
        if taxon_id in self.names_map:
            return self.names_map[taxon_id]
        else:
            return str(taxon_id)  # Fallback to taxon ID if name not found
    
    def lca(self, a: int, b: int) -> int:
        """Find lowest common ancestor of two taxa"""
        if a == 0 or b == 0:
            return a if a else b
            
        # Build ascending path from a
        a_path = set()
        while a != 0 and a in self.parent_map:
            a_path.add(a)
            a = self.parent_map.get(a, 0)
            # Stop if we reach the root (node that is its own parent)
            if a == self.root_id:
                a_path.add(a)
                break
            
        # Walk ascending path from b; check if we hit a's path
        while b != 0 and b in self.parent_map:
            if b in a_path:
                return b
            b = self.parent_map.get(b, 0)
            # Stop if we reach the root
            if b == self.root_id:
                return self.root_id
            
        return self.root_id  # Default to root
    
    def resolve_tree(self, hit_counts: Dict[int, int]) -> int:
        """Resolve classification using root-to-leaf path heuristic"""
        if not hit_counts:
            return 0
            
        max_taxa = set()
        max_taxon, max_score = 0, 0
        
        # Sum each taxon's leaf-to-root path
        for taxon in hit_counts:
            node, score = taxon, 0
            while node != 0 and node in self.parent_map:
                score += hit_counts.get(node, 0)
                node = self.parent_map.get(node, 0)
                # Stop if we reach the root
                if node == self.root_id:
                    score += hit_counts.get(node, 0)
                    break
                
            if score > max_score:
                max_taxa.clear()
                max_score, max_taxon = score, taxon
            elif score == max_score:
                if not max_taxa:
                    max_taxa.add(max_taxon)
                max_taxa.add(taxon)
        
        # If tied, return LCA of all tied taxa
        if len(max_taxa) > 1:
            taxa_list = list(max_taxa)
            result = taxa_list[0]
            for taxon in taxa_list[1:]:
                result = self.lca(result, taxon)
            return result
            
        return max_taxon

class KrakenDB:
    """Simple dictionary structure for canonical-k-mer-to-taxon map"""
    __slots__ = ['k', 'kmer_to_taxon', 'taxonomy', '_canonical_cache']
    
    def __init__(self, k: int = DEF_KMER_LEN, nodes_file: Optional[str] = None):
        self.k = k
        self.kmer_to_taxon: Dict[int, int] = {}
        self.taxonomy = TaxonomyHandler(nodes_file)
        # Cache for canonical representations to avoid recomputation
        self._canonical_cache: Dict[int, int] = {}
    
    def add_kmer(self, kmer: int, taxon: int):
        """Add or update k-mer mapping"""
        canonical = self._get_canonical(kmer)
        if canonical in self.kmer_to_taxon:
            existing = self.kmer_to_taxon[canonical]
            self.kmer_to_taxon[canonical] = self.taxonomy.lca(existing, taxon)
        else:
            self.kmer_to_taxon[canonical] = taxon
    
    def query_kmer(self, kmer: int) -> Optional[int]:
        """Query k-mer and return taxon ID"""
        canonical = self._get_canonical(kmer)
        return self.kmer_to_taxon.get(canonical)
    
    def _get_canonical(self, kmer: int) -> int:
        """Get canonical representation with caching"""
        if kmer in self._canonical_cache:
            return self._canonical_cache[kmer]
        canonical = canonical_representation(kmer, self.k)
        self._canonical_cache[kmer] = canonical
        return canonical
    
    def save(self, filename: str):
        """Save database to file using pickle"""
        data = {
            'k': self.k,
            'kmer_to_taxon': self.kmer_to_taxon
        }
        with open(filename, 'wb') as f:
            pickle.dump(data, f)
    
    def load(self, filename: str):
        """Load database from file using pickle"""
        with open(filename, 'rb') as f:
            data = pickle.load(f)
        self.k = data['k']
        self.kmer_to_taxon = data['kmer_to_taxon']

def read_fasta(filename: str):
    """Simple FASTA reader"""
    sequences = []
    current_seq, current_id = "", ""
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append((current_id, current_seq))
                current_id = line[1:].split()[0]  # Take first word after >
                current_seq = ""
            else:
                current_seq += line
                
    if current_seq:
        sequences.append((current_id, current_seq))
        
    return sequences

def build_database(fasta_files: List[str], taxon_map: Dict[str, int], 
                  nodes_file: Optional[str], k: int = DEF_KMER_LEN) -> KrakenDB:
    """Build Kraken database from FASTA files"""
    db = KrakenDB(k, nodes_file)
    
    print(f"Building database with k={k}...", file=sys.stderr)
    
    # Calculate total file size for progress estimation
    total_size = 0
    for fasta_file in fasta_files:
        if os.path.exists(fasta_file):
            total_size += os.path.getsize(fasta_file)
        else:
            raise FileNotFoundError(f"Library FASTA input not found: {fasta_file}")
    
    print(f"Total input size: {total_size:,} bytes", file=sys.stderr)
    
    processed_bytes, kmer_count = 0, 0
    
    for fasta_file in fasta_files:
        file_size = os.path.getsize(fasta_file)
        print(f"\nProcessing {os.path.basename(fasta_file)} ({file_size:,} bytes)...", file=sys.stderr)
        
        sequences = read_fasta(fasta_file)
        file_processed = 0
        
        for seq_id, sequence in sequences:
            taxon = taxon_map.get(seq_id, 0)
            if taxon == 0:
                print(f"Warning: No taxon mapping for {seq_id}", file=sys.stderr)
                continue
            scanner = KmerScanner(sequence, k)
            total_kmer_count = 0  # Count all k-mers, not just non-ambiguous ones
            
            while True:
                result = scanner.next_kmer()
                if result is None:
                    break
                kmer, has_ambiguous = result
                total_kmer_count += 1
                
                if not has_ambiguous:
                    db.add_kmer(kmer, taxon)
                    kmer_count += 1
                
                if total_kmer_count % 10000 == 0:
                    # Estimate progress based on sequence length processed
                    estimated_bytes = min(file_processed + len(sequence), file_size)
                    print_progress_bar(estimated_bytes, file_size, prefix=f"  {os.path.basename(fasta_file)}")
            
            file_processed += len(sequence)
        
        processed_bytes += file_size
        print_progress_bar(file_size, file_size, prefix=f"  {os.path.basename(fasta_file)}")
    
    print(f"\nDatabase built with {len(db.kmer_to_taxon):,} unique k-mers from {kmer_count:,} total k-mers", file=sys.stderr)
    return db

def classify_sequence(sequence: str, db: KrakenDB, taxonomy: TaxonomyHandler) -> Tuple[int, List]:
    """Classify a sequence"""
    hit_counts = defaultdict(int)
    taxa = []
    scanner = KmerScanner(sequence, db.k)
    
    seq_len = len(sequence) - db.k + 1
    if seq_len <= 0:
        return 0, []
    
    taxa = [None] * seq_len
    taxa_idx = 0
    
    for result in iter(scanner.next_kmer, None):
        kmer, has_ambiguous = result
        if has_ambiguous:
            taxa[taxa_idx] = "A"
        else:
            taxon = db.query_kmer(kmer)
            if taxon:
                hit_counts[taxon] += 1
            taxa[taxa_idx] = taxon or 0
        taxa_idx += 1
    
    # Resolve classification using root-to-leaf path heuristic
    call = taxonomy.resolve_tree(hit_counts)
    return call, taxa

def main():
    parser = argparse.ArgumentParser(description='Minimal Kraken Implementation')
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    build_parser = subparsers.add_parser('build', help='Build database')
    build_parser.add_argument('--db', required=True, help='Database name')
    build_parser.add_argument('--kmer-len', type=int, default=DEF_KMER_LEN, help='K-mer length')
    build_parser.add_argument('--nodes', help='Taxonomy nodes.dmp file')
    build_parser.add_argument('--names', help='Taxonomy names.dmp file')
    build_parser.add_argument('--library', nargs='+', help='FASTA files to include')
    build_parser.add_argument('--taxon-map', help='File mapping sequence IDs to taxon IDs')
    
    classify_parser = subparsers.add_parser('classify', help='Classify sequences')
    classify_parser.add_argument('--db', required=True, help='Database name')
    classify_parser.add_argument('--nodes', help='Taxonomy nodes.dmp file')
    classify_parser.add_argument('--names', help='Taxonomy names.dmp file')
    classify_parser.add_argument('input', help='Input FASTA file')
    
    args = parser.parse_args()
    
    if args.command == 'build':
        taxon_map = {}
        if args.taxon_map:
            with open(args.taxon_map, 'r') as f:
                for line in f:
                    # Allow for slight differences depending on whether
                    # this is prelim_map.txt or taxon_map.txt
                    parts = line.strip().split()
                    if parts[0] == 'TAXID':
                        taxon_map[parts[1]] = int(parts[2])
                    elif len(parts) >= 2:
                        taxon_map[parts[0]] = int(parts[1])
                    else:
                        raise ValueError(f"Invalid taxon map line: {line}")

        db = build_database(args.library or [], taxon_map, args.nodes, args.kmer_len)
        db_name = args.db
        if not db_name.endswith('.db'):
            db_name += '.db'
        db.save(db_name)
        print(f"Database saved to {db_name}", file=sys.stderr)
        
    elif args.command == 'classify':
        db = KrakenDB()
        db_name = args.db
        if not db_name.endswith('.db'):
            db_name += '.db'
        db.load(db_name)
        
        taxonomy = TaxonomyHandler(args.nodes, args.names)
        sequences, classified = read_fasta(args.input), 0
        for seq_id, sequence in sequences:
            call, taxa = classify_sequence(sequence, db, taxonomy)
            # Collapse consecutive identical taxon IDs into runs w/ counts
            if not taxa:
                print(f"U\t{seq_id}\t0\t{len(sequence)}\t")
                continue
                
            collapsed_parts = []
            current_taxon, current_count = taxa[0], 1
            
            for i in range(1, len(taxa)):
                if taxa[i] == current_taxon:
                    current_count += 1
                else:
                    collapsed_parts.append(f"{current_taxon}:{current_count}")
                    current_taxon, current_count = taxa[i], 1
            
            collapsed_parts.append(f"{current_taxon}:{current_count}")
            
            status = 'C' if call else 'U'
            taxa_str = str(call) if call else "0"
            print(f"{status}\t{seq_id}\t{taxa_str}\t{len(sequence)}\t{' '.join(collapsed_parts)}")
            if call:
                classified += 1
        
        print(f"Classified {classified}/{len(sequences)} sequences", file=sys.stderr)
        
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
