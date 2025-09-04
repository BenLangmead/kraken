/*
 * Copyright 2025, Ben Langmead <blangme2@jhu.edu>
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

// This is a short program that converts a KMC database to a Jellyfish
// database, with the goal of allowing KMC to be used instead of
// Jellyfish in the Kraken database building process.  Later, this
// could be leveraged to allow Kraken to support k-mers longer than
// 31-mers.

#include "kraken_headers.hpp"
#include "quickfile.hpp"
#include "krakendb.hpp"
#include "krakenutil.hpp"
#include "seqreader.hpp"
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;
using namespace kraken;

void parse_command_line(int argc, char **argv);
void usage(int exit_code=EX_USAGE);
void convert_kmc_to_jellyfish();
string get_default_temp_directory();

string KMC_prefix, Output_filename;
uint32_t Kmer_length = 0;
bool Verbose = false;
string Temp_directory;

int main(int argc, char **argv) {
  parse_command_line(argc, argv);
  convert_kmc_to_jellyfish();
  return 0;
}

void usage(int exit_code) {
  cerr << "Usage: kmc_to_jellyfish [options] <kmc_prefix> <output.jdb>" << endl
       << endl
       << "Options: (*mandatory)" << endl
       << "  -k <len>     *K-mer length" << endl
       << "  -t <dir>     Temporary directory (default: system temp dir)" << endl
       << "  -v           Verbose output" << endl
       << "  -h           This message" << endl
       << endl
       << "Converts KMC3 database format to Jellyfish v1 database format." << endl
       << "The KMC prefix should point to the .kmc_pre and .kmc_suf files." << endl;
  exit(exit_code);
}

void parse_command_line(int argc, char **argv) {
  int opt;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);

  while ((opt = getopt(argc, argv, "k:t:v")) != -1) {
    switch (opt) {
      case 'k':
        Kmer_length = atoi(optarg);
        if (Kmer_length == 0 || Kmer_length > 31) {
          errx(EX_USAGE, "k-mer length must be between 1 and 31");
        }
        break;
      case 't':
        Temp_directory = optarg;
        break;
      case 'v':
        Verbose = true;
        break;
      default:
        usage();
    }
  }

  if (optind + 2 != argc) {
    errx(EX_USAGE, "must specify KMC prefix and output filename");
  }

  if (Kmer_length == 0) {
    errx(EX_USAGE, "must specify k-mer length with -k");
  }

  KMC_prefix = argv[optind];
  Output_filename = argv[optind + 1];

  // Check if KMC files exist
  const char* suffixes[] = {".kmc_pre", ".kmc_suf"};
  for (const char* suf : suffixes) {
    string fname = KMC_prefix + suf;
    if (access(fname.c_str(), R_OK) != 0) {
      err(EX_NOINPUT, "can't read KMC %s file %s", suf, fname.c_str());
    }
  }

  // Set up temporary directory
  if (Temp_directory.empty()) {
    Temp_directory = get_default_temp_directory();
  }
  
  // Ensure the temporary directory exists and is writable
  if (access(Temp_directory.c_str(), W_OK) != 0) {
    err(EX_CANTCREAT, "temporary directory %s is not writable", Temp_directory.c_str());
  }
}

// Convert DNA string to integer representation
uint64_t string_to_kmer(const string& seq) {
  uint64_t kmer = 0;
  for (size_t i = 0; i < seq.length(); i++) {
    kmer <<= 2;
    switch (seq[i]) {
      case 'A': case 'a': break;
      case 'C': case 'c': kmer |= 1; break;
      case 'G': case 'g': kmer |= 2; break;
      case 'T': case 't': kmer |= 3; break;
      default: errx(EX_DATAERR, "invalid character in sequence: %c", seq[i]);
    }
  }
  return kmer;
}

// Convert integer to DNA string
string kmer_to_string(uint64_t kmer, uint32_t len) {
  string seq;
  for (uint32_t i = 0; i < len; i++) {
    uint64_t base = (kmer >> (2 * (len - 1 - i))) & 3;
    switch (base) {
      case 0: seq += 'A'; break;
      case 1: seq += 'C'; break;
      case 2: seq += 'G'; break;
      case 3: seq += 'T'; break;
    }
  }
  return seq;
}

// Get reverse complement of a k-mer.  It's true that I could have
// reused the implementation in KrakenDB but we want to keep this
// program simple.
uint64_t reverse_complement(uint64_t kmer, uint32_t len) {
  uint64_t rc = 0;
  for (uint32_t i = 0; i < len; i++) {
    rc <<= 2;
    uint64_t base = (kmer >> (2 * i)) & 3;
    rc |= (3 - base);  // A<->T, C<->G
  }
  return rc;
}

// Get canonical representation (lexicographically smaller of k-mer and
// its reverse complement)
uint64_t canonical_representation(uint64_t kmer, uint32_t len) {
  uint64_t rc = reverse_complement(kmer, len);
  return (kmer < rc) ? kmer : rc;
}

// Get a portable default temporary directory
string get_default_temp_directory() {
  // Try environment variables first
  const char* temp_env_vars[] = {"TMPDIR", "TEMP", "TMP"};
  for (const char* var : temp_env_vars) {
    const char* value = getenv(var);
    if (value && strlen(value) > 0) {
      return string(value);
    }
  }
  
  // Fallback to common system directories
  const char* fallback_dirs[] = {"/tmp", "/var/tmp", "/usr/tmp"};
  for (const char* dir : fallback_dirs) {
    if (access(dir, W_OK) == 0) {
      return string(dir);
    }
  }
  
  // Last resort: current directory
  return ".";
}

void convert_kmc_to_jellyfish() {
  if (Verbose) {
    cerr << "Converting KMC database " << KMC_prefix << " to Jellyfish format " << Output_filename << endl;
  }

  // Use kmc_tools to dump the KMC database to text
  string temp_dump = Temp_directory + "/kmc_dump_" + to_string(getpid()) + ".txt";
  string cmd = "kmc_tools transform " + KMC_prefix + " dump " + temp_dump;
  
  if (Verbose) {
    cerr << "Running: " << cmd << endl;
  }
  
  int ret = system(cmd.c_str());
  if (ret != 0) {
    errx(EX_OSERR, "failed to dump KMC database");
  }

  // Convert the dump format to FASTA format for Jellyfish
  string temp_fasta = Temp_directory + "/jellyfish_fasta_" + to_string(getpid()) + ".fa";
  ifstream dump_file(temp_dump);
  if (!dump_file) {
    err(EX_NOINPUT, "can't open temporary dump file %s", temp_dump.c_str());
  }

  ofstream fasta_file(temp_fasta);
  if (!fasta_file) {
    err(EX_CANTCREAT, "can't create FASTA file %s", temp_fasta.c_str());
  }
  
  if (Verbose) {
    cerr << "Created FASTA file: " << temp_fasta << endl;
  }

  string line;
  while (getline(dump_file, line)) {
    if (Verbose) {
      cerr << "Processing line: '" << line << "'" << endl;
    }
    
    istringstream iss(line);
    string kmer_str;
    uint32_t count;
    iss >> kmer_str >> count;
    
    if (Verbose) {
      cerr << "Parsed kmer: '" << kmer_str << "' count: " << count << endl;
    }
    
    if (kmer_str.length() != Kmer_length) {
      errx(EX_DATAERR, "k-mer length mismatch: expected %u, got %zu", Kmer_length, kmer_str.length());
    }
    
    // Write in FASTA format, repeating the k-mer 'count' times
    for (uint32_t i = 0; i < count; i++) {
      fasta_file << ">" << kmer_str << "_" << i << endl;
      fasta_file << kmer_str << endl;
    }
  }
  dump_file.close();
  fasta_file.close();

  if (Verbose) {
    cerr << "Closed files. Checking FASTA file size..." << endl;
    ifstream check_file(temp_fasta);
    check_file.seekg(0, ios::end);
    streamsize size = check_file.tellg();
    cerr << "FASTA file size: " << size << " bytes" << endl;
    check_file.close();
  }

  // Clean up temporary dump file
  unlink(temp_dump.c_str());

  // Use Jellyfish count to create the database
  string temp_prefix = Temp_directory + "/jellyfish_temp_" + to_string(getpid());
  string count_cmd = "jellyfish count -m " + to_string(Kmer_length) + 
                     " -s 1000 -C -o " + temp_prefix + " " + temp_fasta;
  
  if (Verbose) {
    cerr << "Running: " << count_cmd << endl;
  }
  
  ret = system(count_cmd.c_str());
  if (ret != 0) {
    errx(EX_OSERR, "failed to create Jellyfish database (return code: %d)", ret);
  }

  // Check if Jellyfish created any files
  string single_file = temp_prefix + "_0";
  if (access(single_file.c_str(), R_OK) == 0) {
    // Single file created, copy it to the output location
    string copy_cmd = "cp " + single_file + " " + Output_filename;
    ret = system(copy_cmd.c_str());
    if (ret != 0) {
      errx(EX_OSERR, "failed to copy Jellyfish database file from %s to %s", 
           single_file.c_str(), Output_filename.c_str());
    }
  } else {
    // Try to merge multiple files
    string merge_cmd = "jellyfish merge -o " + Output_filename + " " + temp_prefix + "_*";
    ret = system(merge_cmd.c_str());
    if (ret != 0) {
      errx(EX_OSERR, "failed to merge Jellyfish database");
    }
  }

  // Clean up temporary files
  unlink(temp_fasta.c_str());
  string cleanup_cmd = "rm -f " + temp_prefix + "_*";
  system(cleanup_cmd.c_str());

  if (Verbose) {
    cerr << "Successfully converted KMC database to Jellyfish format" << endl;
    cerr << "Output file: " << Output_filename << endl;
  }
}
