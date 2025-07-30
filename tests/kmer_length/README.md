# Kraken k-mer length test data

This directory contains a minimal taxonomy and DNA sequences for
testing Kraken with different k-mer lengths.

Two species are defined under a common genus:

* TaxID 3: Test_species_A (sequence in `library/speciesA.fa`)
* TaxID 4: Test_species_B (individual 31-mers in `library/speciesB.fa`)

The file `reads.fa` contains a single read that matches the genome of
Test_species_A exactly.

When a database is built with **k=31**, each 31-mer from the read exists
in the library of both species. Classification will therefore assign the
read to their lowest common ancestor (TaxID 2, the genus).

When a database is built with **k=63**, only Test_species_A contains the
required 63-mer. The read should then be classified as TaxID 3.

To build the database and run Kraken, use the commands in the k31.sh and
k63.sh files.  They should be run from the same directory as this
README.md file.

The first command should label `read1` with TaxID 2, while the second
should label it with TaxID 3.  The scripts themselves check for these
outcomes and print PASSED with exitlevel 0 if successful, or exit with
nonzero exitlevel otherwise.
