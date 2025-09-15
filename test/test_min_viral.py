#!/usr/bin/env python3

"""
Copyright Ben Langmead <blangme2@jhu.edu>

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


Minimal Viral Kraken Test - Comprehensive comparison between the
minimal_kraken.py script and the original Kraken.  Goal is to use a
non-trivial database and non-trivial set of reads to test whether the
two give identical results.

To run:
1. make -C src
2. python3 test/test_min_viral.py --num-reads 1000 --download-viral

"""

import os
import sys
import tempfile
import argparse
import shutil
from typing import List

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.minimal_kraken import *

# Import required modules for utility functions
import subprocess
import random
from pathlib import Path
from typing import List, Tuple, Dict, Optional

# Utility functions (previously from test_utils)

def run_command(cmd, check=True, capture_output=True, cwd=None):
    """Run a command and return the result."""
    if isinstance(cmd, list):
        cmd_str = ' '.join(cmd)
    else:
        cmd_str = cmd
    print(f"Running: {cmd_str}")
    result = subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=capture_output, text=True, cwd=cwd)
    if check and result.returncode != 0:
        print(f"Command failed with return code {result.returncode}")
        if result.stdout:
            print(f"stdout: {result.stdout}")
        if result.stderr:
            print(f"stderr: {result.stderr}")
        sys.exit(1)
    return result

def load_sequences(fasta_file: str) -> List[Tuple[str, str]]:
    """Load sequences from FASTA file"""
    sequences = []
    current_seq = ""
    current_id = ""
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append((current_id, current_seq))
                current_id = line[1:].split()[0]
                current_seq = ""
            else:
                current_seq += line
                
    if current_seq:
        sequences.append((current_id, current_seq))
        
    return sequences

def mutate_sequence(sequence: str, mutation_rate: float) -> str:
    """Apply random mutations to a sequence"""
    bases = ['A', 'C', 'G', 'T']
    mutated = list(sequence)
    
    for i in range(len(mutated)):
        if random.random() < mutation_rate:
            # Random substitution
            original = mutated[i]
            new_bases = [b for b in bases if b != original]
            mutated[i] = random.choice(new_bases)
    
    return ''.join(mutated)

def generate_simulated_reads(sequences: List[Tuple[str, str]], num_reads: int = 50) -> List[Tuple[str, str, str, float]]:
    """Generate simulated reads with different mutation rates"""
    reads = []
    
    # Random reads (should not classify)
    print(f"Generating {num_reads//4} random reads...")
    for i in range(num_reads // 4):
        random_seq = ''.join(random.choices(['A', 'C', 'G', 'T'], k=100))
        reads.append((f"random_{i}", random_seq, "random", 0.0))
    
    # Reads from sequences with different mutation rates
    mutation_rates = [0.01, 0.05, 0.10]  # 1%, 5%, 10%
    reads_per_rate = num_reads // (4 * len(mutation_rates))
    
    for rate in mutation_rates:
        print(f"Generating {reads_per_rate} reads with {rate*100}% mutation rate...")
        for i in range(reads_per_rate):
            # Pick a random sequence
            seq_id, sequence = random.choice(sequences)
            
            # Extract a random 100bp read
            if len(sequence) < 100:
                read_seq = sequence
            else:
                start = random.randint(0, len(sequence) - 100)
                read_seq = sequence[start:start+100]
            
            # Apply mutations
            mutated_seq = mutate_sequence(read_seq, rate)
            
            reads.append((f"viral_{rate*100}pct_{i}", mutated_seq, seq_id, rate))
    
    return reads

def write_reads_to_fasta(reads: List[Tuple[str, str, str, float]], filename: str):
    """Write reads to FASTA file"""
    with open(filename, 'w') as f:
        for read_id, sequence, source, mutation_rate in reads:
            f.write(f">{read_id} source={source} mutation_rate={mutation_rate}\n")
            f.write(f"{sequence}\n")

def parse_kraken_output(output_lines: List[str]) -> Dict[str, Dict[str, str]]:
    """Parse Kraken output into a dictionary with all fields"""
    results = {}
    for line in output_lines:
        if line.strip() and '\t' in line:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                read_id = parts[1]
                results[read_id] = {
                    'status': parts[0],           # C or U
                    'seq_id': parts[1],           # Sequence ID
                    'taxon_id': parts[2],         # Taxon ID or 0
                    'length': parts[3],           # Sequence length
                    'lca_mapping': parts[4]       # Colon-separated taxon IDs
                }
    return results

def compare_classifications_detailed(minimal_results: Dict[str, Dict[str, str]], 
                                   original_results: Dict[str, Dict[str, str]], 
                                   output_file: Optional[str] = None) -> Dict[str, int]:
    """Compare classification results with detailed field-by-field analysis"""
    comparison = {
        'total_reads': 0,
        'exact_matches': 0,
        'status_differences': 0,
        'taxon_id_differences': 0,
        'length_differences': 0,
        'lca_mapping_differences': 0,
        'minimal_only': 0,
        'original_only': 0
    }
    
    differences = []  # Store differences for output file
    
    all_reads = set(minimal_results.keys()) | set(original_results.keys())
    comparison['total_reads'] = len(all_reads)
    
    for read_id in all_reads:
        minimal_data = minimal_results.get(read_id, {})
        original_data = original_results.get(read_id, {})
        
        # Handle reads present in only one result set
        if not minimal_data:
            comparison['original_only'] += 1
            if output_file:
                differences.append(('original', original_data))
            continue
        elif not original_data:
            comparison['minimal_only'] += 1
            if output_file:
                differences.append(('minimal', minimal_data))
            continue
        
        # Compare all fields
        status_match = minimal_data.get('status') == original_data.get('status')
        taxon_match = minimal_data.get('taxon_id') == original_data.get('taxon_id')
        length_match = minimal_data.get('length') == original_data.get('length')
        lca_match = minimal_data.get('lca_mapping') == original_data.get('lca_mapping')
        
        if status_match and taxon_match and length_match and lca_match:
            comparison['exact_matches'] += 1
        else:
            # Track specific differences
            if not status_match:
                comparison['status_differences'] += 1
            if not taxon_match:
                comparison['taxon_id_differences'] += 1
            if not length_match:
                comparison['length_differences'] += 1
            if not lca_match:
                comparison['lca_mapping_differences'] += 1
            
            # Store differences for output file
            if output_file:
                differences.append(('min', minimal_data))
                differences.append(('org', original_data))
    
    # Write differences to output file if requested
    if output_file:
        with open(output_file, 'w') as f:
            if differences:
                f.write("# Differences between minimal and original Kraken\n")
                f.write("# Format: tool\tstatus\tseq_id\ttaxon_id\tlength\tlca_mapping\n")
                for tool, data in differences:
                    if data:  # Only write if data exists
                        f.write(f"{tool}\t{data.get('status', 'N/A')}\t{data.get('seq_id', 'N/A')}\t"
                            f"{data.get('taxon_id', 'N/A')}\t{data.get('length', 'N/A')}\t"
                            f"{data.get('lca_mapping', 'N/A')}\n")
            else:
                f.write("# No differences found\n")
    
    return comparison

def print_tree(path, prefix=""):
    if os.path.isfile(path):
        print(f"{prefix}{os.path.basename(path)}")
    elif os.path.isdir(path):
        print(f"{prefix}{os.path.basename(path)}/")
        entries = sorted(os.listdir(path))
        for i, entry in enumerate(entries):
            full_path = os.path.join(path, entry)
            is_last = (i == len(entries) - 1)
            new_prefix = prefix + ("    " if is_last else "│   ")
            print_tree(full_path, prefix + ("└── " if is_last else "├── "))

def create_viral_database(db_inp: str):
    """Download viral database from NCBI and write to .fna file """
    os.makedirs(db_inp, exist_ok=True)
    
    # Check if files already exist
    combined_fasta = os.path.join(db_inp, 'viral_genomes.fna')
    prelim_map_file = os.path.join(db_inp, "prelim_map.txt")
    nodes_file = os.path.join(db_inp, "nodes.dmp")
    names_file = os.path.join(db_inp, "names.dmp")
    
    required_files = [combined_fasta, prelim_map_file, nodes_file, names_file]
    all_files_exist = all(os.path.exists(f) for f in required_files)
    
    if all_files_exist:
        print(f"Viral database files already exist in {db_inp}, skipping download...")
        return
    
    print("Downloading & writing viral references from NCBI...")
    
    # You can copy and paste these (or upload test_min_viral_taxa.txt) into
    # https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi
    # And that gives the tree found in test_min_viral_taxtree.txt
    viral_genomes = {
        "NC_001416.1": "2681611",  # Enterobacteria phage lambda
        "NC_001422.1": "2886930",  # Escherichia phage phiX174
        "NC_001604.1": "10760",    # Enterobacteria phage T7
        "NC_001405.1": "129951",   # Human adenovirus C
        "NC_001357.1": "333761",   # Human papillomavirus 18
        "NC_001612.1": "138948",   # Human enterovirus A
        "NC_001802.1": "11676",    # HIV-1
        "NC_001479.1": "12104",    # Encephalomyocarditis virus
        "NC_001489.1": "12092",    # Hepatitis A virus
        "NC_001477.1": "11053",    # Dengue virus 1
        "NC_001474.1": "11060",    # Dengue virus 2
        "NC_001475.1": "11069",    # Dengue virus 3
        "NC_001417.1": "12022",    # Escherichia phage MS2
    }
    
    URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    with open(combined_fasta, 'w') as fasta_fh, open(prelim_map_file, 'w') as prelim_fh:
        for accession, taxon_id in viral_genomes.items():
            print(f"  Downloading {accession} (taxon {taxon_id})...")
            cmd = [
                "curl", "-s", 
                f"{URL}?db=nuccore&id={accession}&rettype=fasta&retmode=text"
            ]
            result = run_command(cmd)
            if result.returncode != 0:
                raise RuntimeError(f"Failed to download {accession}")
            fasta_fh.write(result.stdout)
            prelim_fh.write(f"TAXID\t{accession}\t{taxon_id}\n")
    
    # Create taxonomy files if they don't exist
    for fname, srcname in [("nodes.dmp", "test_min_viral_nodes.dmp"),
                           ("names.dmp", "test_min_viral_names.dmp")]:
        fpath = os.path.join(db_inp, fname)
        if not os.path.exists(fpath):
            src = os.path.join(os.path.dirname(__file__), 'test_minimal_files', srcname)
            with open(src, 'r') as s, open(fpath, 'w') as f:
                for line in s:
                    f.write(line.replace(' | ', '\t|\t').replace(' |', '\t|'))

def use_existing_viral_database(viral_db_path: str):
    """Use existing viral database"""
    print(f"Using existing viral database: {viral_db_path}")
    
    combined_fasta = os.path.join(viral_db_path, 'viral_genomes.fna')
    nodes_file = os.path.join(viral_db_path, "nodes.dmp")
    names_file = os.path.join(viral_db_path, "names.dmp")
    
    if (not os.path.exists(combined_fasta) or 
        not os.path.exists(nodes_file) or 
        not os.path.exists(names_file)):
        raise RuntimeError(f"Required files not found in {viral_db_path}")

def setup_viral_database(viral_db_path: str, download_fresh: bool = False):
    """Set up viral database - either use existing or download fresh"""
    if download_fresh:
        create_viral_database(viral_db_path)
    else:
        use_existing_viral_database(viral_db_path)

def setup_library_structure(db_inp: str, working_dir: str) -> str:
    """Set up proper library structure for original Kraken"""
    print(f"Setting up library structure for db_dir={working_dir}...")
    
    if not os.path.exists(db_inp):
        raise RuntimeError(f"Input dir not found: {db_inp}")

    src_fasta = os.path.join(db_inp, 'viral_genomes.fna')
    if not os.path.exists(src_fasta):
        raise RuntimeError(f"Source library file not found: {src_fasta}")

    db_dir = os.path.join(working_dir, 'original_kraken_db')
    dst_library_dir = os.path.join(db_dir, "library")
    dst_taxonomy_dir = os.path.join(db_dir, "taxonomy")
    os.makedirs(dst_library_dir, exist_ok=True)
    os.makedirs(dst_taxonomy_dir, exist_ok=True)
    # Please don't change this to taxon_map.txt, no matter how tempting it is
    for fn in ['viral_genomes.fna', 'prelim_map.txt']:
        shutil.copy2(os.path.join(db_inp, fn), os.path.join(dst_library_dir, fn))

    # Please don't change this to taxon_map.txt, no matter how tempting it is
    for fn in ['nodes.dmp', 'names.dmp', 'prelim_map.txt']:
        src = os.path.join(db_inp, fn)
        if not os.path.exists(src):
            raise RuntimeError(f"File not found: {src}")
        shutil.copy2(src, os.path.join(dst_taxonomy_dir, fn))

def build_minimal_kraken_database(temp_dir: str, db_inp: str, k: int = 31) -> str:
    """Build database using minimal Kraken"""
    print("Building database with minimal Kraken...")
    
    db_name = os.path.join(temp_dir, 'minimal_kraken.db')
    combined_fasta = os.path.join(db_inp, 'viral_genomes.fna')
    nodes_file = os.path.join(db_inp, "nodes.dmp")
    names_file = os.path.join(db_inp, "names.dmp")
    map_file = os.path.join(db_inp, "prelim_map.txt")
    
    # Get absolute path to minimal_kraken.py
    script_dir = os.path.dirname(os.path.abspath(__file__))
    minimal_kraken_path = os.path.join(os.path.dirname(script_dir), "scripts", "minimal_kraken.py")
    if not os.path.exists(minimal_kraken_path):
        print(f"minimal_kraken.py not found at {minimal_kraken_path}")
        return None
    
    cmd = [
        sys.executable, minimal_kraken_path, "build",
        "--db", db_name, "--library", combined_fasta,
        "--taxon-map", map_file,
        "--nodes", nodes_file, "--names", names_file,
        "--kmer-len", str(k)
    ]
    
    result = run_command(cmd)
    if result.returncode != 0:
        raise RuntimeError(f"Minimal Kraken build failed: {result.stderr}")
    
    print(f"Minimal Kraken database built: {db_name}.db")
    return db_name

def build_original_kraken_database(db_inp: str, working_dir: str, k: int = 31, minimizer_len: int = 10) -> str:
    """Build database using original Kraken with proper library structure"""
    print("Building database with original Kraken...")
    
    db_name = os.path.join(working_dir, 'original_kraken_db')
    
    # Set up library structure
    setup_library_structure(db_inp, working_dir)
    
    # Get absolute path to kraken-build
    script_dir = os.path.dirname(os.path.abspath(__file__))
    kraken_build_path = os.path.join(os.path.dirname(script_dir), "scripts", "kraken-build")
    
    # Check if kraken-build exists
    if not os.path.exists(kraken_build_path):
        raise RuntimeError(f"kraken-build not found at {kraken_build_path}")

    print(f"Building database with k={k}, minimizer_len={minimizer_len}, db={db_name}...")
    cmd = [
        kraken_build_path, "--db", db_name,
        "--build", "--kmer-len", str(k),
        "--minimizer-len", str(minimizer_len)
    ]
    
    result = run_command(cmd)
    if result.returncode != 0:
        raise RuntimeError(f"Original Kraken build failed: {result.stderr}")
    
    print(f"Original Kraken database built: {db_name}")
    return db_name

def classify_with_minimal_kraken(db_name: str, db_inp: str, reads_file: str) -> List[str]:
    """Classify reads using minimal Kraken"""
    # Get absolute path to minimal_kraken.py
    script_dir = os.path.dirname(os.path.abspath(__file__))
    minimal_kraken_path = os.path.join(os.path.dirname(script_dir), "scripts", "minimal_kraken.py")

    # I'm leaving nodes.dmp and names.dmp in the db_input directory;
    # would also be reasonable to bundle the minimal kraken db in a
    # directory like a standard kraken db and copy them into that
    nodes_file = os.path.join(db_inp, "nodes.dmp")
    names_file = os.path.join(db_inp, "names.dmp")
    if not db_name.endswith('.db'):
        db_name += '.db'
    
    cmd = [
        sys.executable, minimal_kraken_path, "classify",
        "--db", db_name, "--nodes", nodes_file, "--names", names_file,
        reads_file
    ]
    
    result = run_command(cmd)
    if result.returncode != 0:
        print(f"Minimal Kraken classification failed: {result.stderr}")
        return []
    
    return result.stdout.strip().split('\n')

def classify_with_original_kraken(db_name: str, reads_file: str) -> List[str]:
    """Classify reads using original Kraken"""
    print("Classifying reads with original Kraken...")
    
    # Get absolute path to classify (the main Kraken binary)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    kraken_path = os.path.join(os.path.dirname(script_dir), "src", "classify")
    
    cmd = [
        kraken_path, "-d", os.path.join(db_name, "database.kdb"),
        "-i", os.path.join(db_name, "database.idx"),
        "-n", os.path.join(db_name, "taxonomy", "nodes.dmp"),
        "-t", "1",
        reads_file
    ]
    
    result = run_command(cmd)
    if result.returncode != 0:
        print(f"Original Kraken classification failed: {result.stderr}")
        return []
    
    return result.stdout.strip().split('\n')


def main():
    parser = argparse.ArgumentParser(description='Minimal Viral Kraken Test script')
    parser.add_argument('--num-reads', type=int, default=1000, help='Number of simulated reads')
    parser.add_argument('--kmer-len', type=int, default=31, help='K-mer length')
    parser.add_argument('--minimizer-len', type=int, default=10, help='Minimizer length (default: 10)')
    parser.add_argument('--viral-db', type=str, default='../data/viral_library', help='Path to viral database')
    parser.add_argument('--download-viral', action='store_true', help='Download fresh viral genomes from NCBI')
    parser.add_argument('--skip-original', action='store_true', help='Skip original Kraken (faster)')
    parser.add_argument('--timeout', type=int, default=300, help='Timeout for original Kraken build (seconds)')
    parser.add_argument('--differences-file', type=str, help='Output file for detailed differences between minimal and original Kraken')
    parser.add_argument('--work-dir', type=str, help='Custom working directory (overrides temporary directory)')
    args = parser.parse_args()
    
    print("=== Minimal Viral Kraken Test ===")
    print(f"# reads: {args.num_reads}")
    print(f"K-mer length: {args.kmer_len}")
    print(f"Minimizer length: {args.minimizer_len}")
    print(f"Existing DB: {args.viral_db}")
    print(f"Download fresh?: {args.download_viral}")
    print(f"Skip original Kraken: {args.skip_original}")
    print(f"Custom work directory: {args.work_dir if args.work_dir else 'Using temporary directory'}")
    
    # Set up working directory
    temp_dir_context = None
    try:
        if args.work_dir:
            # Use custom working directory
            temp_dir = args.work_dir
            os.makedirs(temp_dir, exist_ok=True)
            print(f"Working directory: {temp_dir}")
        else:
            # Use temporary directory (will be cleaned up automatically)
            temp_dir_context = tempfile.TemporaryDirectory()
            temp_dir = temp_dir_context.name
            print(f"Working directory: {temp_dir}")
        
        print("Step 1: Set up database preliminaries")
        if args.download_viral:
            db_inp = os.path.join(temp_dir, "db_inputs")
            setup_viral_database(db_inp, download_fresh=True)
        else:
            db_inp = args.viral_db
            setup_viral_database(args.viral_db, download_fresh=False)

        # Step 2: Load sequences and generate reads
        print("Step 2: Loading viral sequences")
        sequences, nfiles = [], 0
        for fname in os.listdir(db_inp):
            if fname.endswith('.fa') or fname.endswith('.fna'):
                nfiles += 1
                fpath = os.path.join(db_inp, fname)
                sequences.extend(load_sequences(fpath))
        # TODO: Actually, we should insist that the db_inp directory
        # have a single combined .fna file
        if len(sequences) == 0:
            raise RuntimeError("No sequences loaded")
        print(f"Loaded {len(sequences)} viral sequences from {nfiles} files")
        
        print("Step 3: Generating simulated reads")
        reads = generate_simulated_reads(sequences, args.num_reads)
        reads_file = os.path.join(temp_dir, "simulated_reads.fa")
        write_reads_to_fasta(reads, reads_file)
        print(f"Generated {len(reads)} simulated reads")
        
        print("Step 4: Build minimal Kraken database")
        minimal_db = build_minimal_kraken_database(temp_dir, db_inp, args.kmer_len)
        if not minimal_db:
            raise RuntimeError("Failed to build minimal Kraken database")
        
        print(f"Step 5: Minimal Kraken classification")
        minimal_output = classify_with_minimal_kraken(minimal_db, db_inp, reads_file)
        minimal_results = parse_kraken_output(minimal_output)
        print(f"Minimal Kraken classified {len(minimal_results)} reads")
        
        # Save minimal Kraken output to file
        minimal_output_file = os.path.join(temp_dir, "minimal_kraken_output.txt")
        with open(minimal_output_file, 'w') as f:
            f.write('\n'.join(minimal_output))
        print(f"Minimal Kraken output saved to: {minimal_output_file}")
        
        # Step 5: Try to classify with original Kraken (if not skipped)
        original_results = {}
        if not args.skip_original:
            print("Directory structure before building original Kraken database:")
            print_tree(temp_dir)
            print("Attempting to build original Kraken database...")
            try:
                original_db = build_original_kraken_database(
                    db_inp, temp_dir, args.kmer_len, args.minimizer_len
                )
                if original_db:
                    original_output = classify_with_original_kraken(original_db, reads_file)
                    original_results = parse_kraken_output(original_output)
                    print(f"Original Kraken classified {len(original_results)} reads")
                    
                    # Save original Kraken output to file
                    original_output_file = os.path.join(temp_dir, "original_kraken_output.txt")
                    with open(original_output_file, 'w') as f:
                        f.write('\n'.join(original_output))
                    print(f"Original Kraken output saved to: {original_output_file}")
                else:
                    print("Original Kraken database build failed, continuing with minimal Kraken only")
            except Exception as e:
                print(f"Original Kraken failed with exception: {e}")
                print("Continuing with minimal Kraken only")
        
        # Step 6: Compare results
        if original_results:
            # Use detailed comparison with optional differences file
            comparison = compare_classifications_detailed(minimal_results, original_results, args.differences_file)
            
            print("\n=== Detailed Classification Comparison ===")
            print(f"Total reads: {comparison['total_reads']}")
            print(f"Exact matches: {comparison['exact_matches']}")
            print(f"Status differences: {comparison['status_differences']}")
            print(f"Taxon ID differences: {comparison['taxon_id_differences']}")
            print(f"Length differences: {comparison['length_differences']}")
            print(f"LCA mapping differences: {comparison['lca_mapping_differences']}")
            print(f"Minimal only: {comparison['minimal_only']}")
            print(f"Original only: {comparison['original_only']}")
            
            exact_match_rate = (comparison['exact_matches'] / comparison['total_reads'] * 100) if comparison['total_reads'] > 0 else 0
            print(f"\nExact match rate: {exact_match_rate:.1f}%")
            
            if comparison['taxon_id_differences'] > 0:
                taxon_diff_rate = (comparison['taxon_id_differences'] / comparison['total_reads'] * 100)
                print(f"Taxon ID difference rate: {taxon_diff_rate:.1f}%")
            
            if args.differences_file:
                print(f"\nDetailed differences written to: {args.differences_file}")
        else:
            # Analyze minimal Kraken results only
            print("\n=== Minimal Kraken Results ===")
            classified = sum(1 for v in minimal_results.values() if v.get('status') != 'U')
            unclassified = len(minimal_results) - classified
            print(f"Classified: {classified}")
            print(f"Unclassified: {unclassified}")
            
            # Analyze by read type
            random_classified = sum(1 for read_id, result in minimal_results.items() 
                                  if read_id.startswith('random_') and result.get('status') != 'U')
            viral_classified = sum(1 for read_id, result in minimal_results.items() 
                                 if read_id.startswith('viral_') and result.get('status') != 'U')
            
            total_random = sum(1 for read_id in minimal_results.keys() if read_id.startswith('random_'))
            total_viral = sum(1 for read_id in minimal_results.keys() if read_id.startswith('viral_'))
            
            print(f"\n=== Classification Analysis ===")
            if total_random > 0:
                print(f"Random reads: {random_classified}/{total_random} classified ({random_classified/total_random*100:.1f}%)")
            if total_viral > 0:
                print(f"Viral reads: {viral_classified}/{total_viral} classified ({viral_classified/total_viral*100:.1f}%)")
            
            if viral_classified / total_viral >= 0.5 if total_viral > 0 else True:
                print("GOOD: Reasonable classification rate for viral reads")
            
            print("PASS: Minimal Kraken classification completed")
            return 0
    
    finally:
        # Clean up temporary directory only if we created one
        if temp_dir_context is not None:
            temp_dir_context.cleanup()
            print(f"Cleaned up temporary directory: {temp_dir}")
        elif args.work_dir:
            print(f"Preserved custom working directory: {temp_dir}")

if __name__ == '__main__':
    sys.exit(main())
