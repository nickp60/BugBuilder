get_genomes.py -q NC_002951.2 -o ./tests/references/
art_illumina -ss HS25 -sam -i ./tests/references/NC_002951.2.fasta -p -l 100 -f 20 -m 200 -s 10 -o tests/references/NC_002951_
#THE reference\n  get_genomes.py -q AP017922.1 -o ./tests/references/
art_illumina -ss HS25 -sam -i ./tests/references/NC_002951.2.fasta -p -l 100 -f 20 -m 200 -s 10 -o tests/references/full_
