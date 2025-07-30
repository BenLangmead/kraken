#!/bin/bash

set -ex

K=63
EXPECTED_TAXON=3

for fn in library/speciesA.fa library/speciesB.fa ; do
    ../../dist/kraken-build --db db${K} --add-to-library ${fn} \
            --jellyfish-hash-size 1M --kmer-len ${K} --minimizer-len 15
done

cp -r taxonomy db${K}/
export KRAKEN_DIR=/code/dist && cd db${K} && ../../../dist/scan_fasta_file.pl library/species* > library/prelim_map.txt && cd ..
../../dist/kraken-build --db db${K} --jellyfish-hash-size 1M --build --kmer-len ${K} --minimizer-len 15
../../dist/kraken --db db${K} reads.fa > k${K}_out.txt
grep -q "read1.${EXPECTED_TAXON}" k${K}_out.txt && echo "PASSED" && exit 0
echo "FAILED"
exit 1


