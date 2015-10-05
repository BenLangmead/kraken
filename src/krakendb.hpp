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

#ifndef KRAKENDB_HPP
#define KRAKENDB_HPP

#include "kraken_headers.hpp"

// XOR mask for minimizer bin keys (allows for better distribution)
// scrambles minimizer sort order
static const uint64_t INDEX2_XOR_MASK = 0xe37e28c4271b5a2dULL;

namespace kraken {
  class KrakenDBIndex {
    public:
    KrakenDBIndex();
    // ptr points to mmap'ed existing file opened in read or read/write mode
    KrakenDBIndex(char *ptr);

    uint8_t index_type();
    uint8_t indexed_nt();
    uint64_t *get_array();
    uint64_t at(uint64_t idx);

    private:
    uint8_t idx_type;
    char *fptr;
    uint8_t nt;
  };

  class KrakenDB {
    public:

    char *get_ptr();            // Return the file pointer
    char *get_pair_ptr();       // Return pointer to start of pairs
    KrakenDBIndex *get_index(); // Return ptr to assoc'd index obj
    uint8_t get_k();            // how many nt are in each key?
    uint64_t get_key_bits();    // how many bits are in each key?
    uint64_t get_key_len();     // how many bytes does each key occupy?
    uint64_t get_val_len();     // how many bytes does each value occupy?
    uint64_t get_key_ct();      // how many key/value pairs are there?
    uint64_t pair_size();       // how many bytes does each pair occupy?

    size_t header_size();  // Jellyfish uses variable header sizes
    uint32_t *kmer_query(uint64_t kmer);  // return ptr to pair w/ kmer

    // perform search over last range to speed up queries
    uint32_t *kmer_query(uint64_t kmer, uint64_t *last_bin_key,
                         int64_t *min_pos, int64_t *max_pos,
                         bool retry_on_failure=true);
    
    // return "bin key" for kmer, based on index
    // If idx_nt not specified, use index's value
    inline uint64_t bin_key(uint64_t kmer, uint64_t idx_nt) {
        uint8_t nt = idx_nt;
        uint64_t xor_mask = INDEX2_XOR_MASK;
        uint64_t mask = 1 << (nt * 2);
        mask--;
        xor_mask &= mask;
        uint64_t min_bin_key = ~0;
        for (uint64_t i = 0; i < key_bits / 2 - nt + 1; i++) {
            uint64_t temp_bin_key = xor_mask ^ canonical_representation(kmer & mask, nt);
            if (temp_bin_key < min_bin_key)
                min_bin_key = temp_bin_key;
            kmer >>= 2;
        }
        return min_bin_key;
    }

    inline uint64_t bin_key(uint64_t kmer) {
        uint8_t nt = index_ptr->indexed_nt();
        uint8_t idx_type = index_ptr->index_type();
        uint64_t xor_mask = idx_type == 1 ? 0 : INDEX2_XOR_MASK;
        uint64_t mask = 1 << (nt * 2);
        mask--;
        xor_mask &= mask;
        uint64_t min_bin_key = ~0;
        for (uint64_t i = 0; i < key_bits / 2 - nt + 1; i++) {
            uint64_t temp_bin_key = xor_mask ^ canonical_representation(kmer & mask, nt);
            if (temp_bin_key < min_bin_key)
                min_bin_key = temp_bin_key;
            kmer >>= 2;
        }
        return min_bin_key;
    }

    // Code from Jellyfish, rev. comp. of a k-mer with n nt.
    // If n is not specified, use k in DB, otherwise use first n nt in kmer
    inline uint64_t reverse_complement(uint64_t kmer, uint8_t n) {
        kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
        kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
        kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
        kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
        kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
        return (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - (n << 1));
    }

    inline uint64_t reverse_complement(uint64_t kmer) {
        kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
        kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
        kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
        kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
        kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
        return (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - (k << 1));
    }

    // Return lexicographically smallest of kmer/revcom(kmer)
    // If n is not specified, use k in DB, otherwise use first n nt in kmer
    inline uint64_t canonical_representation(uint64_t kmer, uint8_t n) {
        uint64_t revcom = reverse_complement(kmer, n);
        return kmer < revcom ? kmer : revcom;
    }

    inline uint64_t canonical_representation(uint64_t kmer) {
       uint64_t revcom = reverse_complement(kmer, k);
       return kmer < revcom ? kmer : revcom;
    }

    void make_index(std::string index_filename, uint8_t nt);

    void set_index(KrakenDBIndex *i_ptr);

    // Null constructor
    KrakenDB();

    // ptr points to start of mmap'ed DB in read or read/write mode
    KrakenDB(char *ptr);

    private:

    char *fptr;
    KrakenDBIndex *index_ptr;
    uint8_t k;
    uint64_t key_bits;
    uint64_t key_len;
    uint64_t val_len;
    uint64_t key_ct;
  };
}

#endif
