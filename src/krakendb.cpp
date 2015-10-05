/*
 * Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken taxonomic sequence classification system.
 *
 * Kraken is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Kraken is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kraken.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "kraken_headers.hpp"
#include "krakendb.hpp"
#include "quickfile.hpp"

using std::string;
using std::vector;

namespace kraken {

// File type code for Jellyfish/Kraken DBs
static const char * DATABASE_FILE_TYPE = "JFLISTDN";

// File type code on Kraken DB index
// Next byte determines # of indexed nt
static const char * KRAKEN_INDEX_STRING = "KRAKIDX";

// File type code for Kraken DB index (v2)
// v2 index corresponds to DB sorted by scrambled order
// Next byte determines # of indexed nt
static const char * KRAKEN_INDEX2_STRING = "KRAKIX2";

// Basic constructor
KrakenDB::KrakenDB() {
  fptr = NULL;
  index_ptr = NULL;
  key_ct = 0;
  val_len = 0;
  key_len = 0;
  key_bits = 0;
  k = 0;
}

// Assumes ptr points to start of a readable mmap'ed file
KrakenDB::KrakenDB(char *ptr) {
  index_ptr = NULL;
  fptr = ptr;
  if (strncmp(ptr, DATABASE_FILE_TYPE, strlen(DATABASE_FILE_TYPE)))
    errx(EX_DATAERR, "database in improper format");
  memcpy(&key_bits, ptr + 8, 8);
  memcpy(&val_len, ptr + 16, 8);
  memcpy(&key_ct, ptr + 48, 8);
  if (val_len != 4)
    errx(EX_DATAERR, "can only handle 4 byte DB values");
  k = key_bits / 2;
  key_len = key_bits / 8 + !! (key_bits % 8);
}

// Creates an index, indicating starting positions of each bin
// Bins contain k-mer/taxon pairs with k-mers that share a bin key
void KrakenDB::make_index(string index_filename, uint8_t nt) {
  uint64_t entries = 1ull << (nt * 2);
  vector<uint64_t> bin_counts(entries);
  char *ptr = get_pair_ptr();
  #pragma omp parallel for schedule(dynamic,400)
  for (uint64_t i = 0; i < key_ct; i++) {
    uint64_t kmer = 0;
    memcpy(&kmer, ptr + i * pair_size(), key_len);
    uint64_t b_key = bin_key(kmer, nt);
    #pragma omp atomic
    bin_counts[b_key]++;
  }

  uint64_t *bin_offsets = new uint64_t[ entries + 1 ];
  bin_offsets[0] = 0;
  for (uint64_t i = 1; i <= entries; i++)
    bin_offsets[i] = bin_offsets[i-1] + bin_counts[i-1];

  QuickFile idx_file(index_filename, "w",
    strlen(KRAKEN_INDEX2_STRING) + 1 + sizeof(*bin_offsets) * (entries + 1));
  char *idx_ptr = idx_file.ptr();
  memcpy(idx_ptr, KRAKEN_INDEX2_STRING, strlen(KRAKEN_INDEX2_STRING));
  idx_ptr += strlen(KRAKEN_INDEX2_STRING);
  memcpy(idx_ptr++, &nt, 1);
  memcpy(idx_ptr, bin_offsets, sizeof(*bin_offsets) * (entries + 1));
}

// Simple accessor
char *KrakenDB::get_ptr() {
  return fptr;
}

// Returns start of k-mer/taxon pair array (skips header)
char *KrakenDB::get_pair_ptr() {
  return fptr == NULL ? NULL : fptr + header_size();
}

// Simple accessor
KrakenDBIndex *KrakenDB::get_index() {
  return index_ptr;
}

// Associates the index with this database
void KrakenDB::set_index(KrakenDBIndex *i_ptr) {
  index_ptr = i_ptr;
}

// Simple accessors/convenience methods
uint8_t KrakenDB::get_k() { return k; }
uint64_t KrakenDB::get_key_bits() { return key_bits; }
uint64_t KrakenDB::get_key_len() { return key_len; }
uint64_t KrakenDB::get_val_len() { return val_len; }
uint64_t KrakenDB::get_key_ct() { return key_ct; }
uint64_t KrakenDB::pair_size() { return key_len + val_len; }
size_t KrakenDB::header_size() { return 72 + 2 * (4 + 8 * key_bits); }

KrakenDBIndex::KrakenDBIndex() {
  fptr = NULL;
  idx_type = 1;
  nt = 0;
}

KrakenDBIndex::KrakenDBIndex(char *ptr) {
  fptr = ptr;
  idx_type = 1;
  if (strncmp(ptr, KRAKEN_INDEX_STRING, strlen(KRAKEN_INDEX_STRING))) {
    idx_type = 2;
    if (strncmp(ptr, KRAKEN_INDEX2_STRING, strlen(KRAKEN_INDEX2_STRING)))
      errx(EX_DATAERR, "illegal Kraken DB index format");
  }
  ptr += strlen(KRAKEN_INDEX_STRING);
  memcpy(&nt, ptr, 1);
}

// Index version (v2 uses different minimizer sort order)
uint8_t KrakenDBIndex::index_type() {
  return idx_type;
}

// How long are bin keys (i.e., what is minimizer length in bp?)
uint8_t KrakenDBIndex::indexed_nt() {
  return nt;
}

// Return start of index array (skips header)
uint64_t *KrakenDBIndex::get_array() {
  return (uint64_t *) (fptr + strlen(KRAKEN_INDEX_STRING) + 1);
}

} // namespace
