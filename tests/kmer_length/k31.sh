#!/bin/bash

set -ex

for fn in library/speciesA.fa library/speciesB.fa ; do
    ../../dist/kraken-build --db db31 --add-to-library ${fn} \
            --jellyfish-hash-size 1M --kmer-len 31 --minimizer-len 15
done

cp -r taxonomy db31/
KRAKEN_DIR=/code/dist && cd db31 && ../../../dist/scan_fasta_file.pl library/species* > library/prelim_map.txt && cd ..
../../dist/kraken-build --db db31 --jellyfish-hash-size 1M --build --kmer-len 31 --minimizer-len 15
../../dist/kraken --db db31 reads.fa
