#!/bin/bash
set -e


if [ -z "$refdir" ]
then
    refdir=$1
    echo "saving results to $refdir"
else
    refdir=./tests/references/
fi


if command -v art_illumina  2>/dev/null
then
    echo "art found"
else
    conda install art
fi

mkdir $refdir
cd $refdir
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/985/GCF_000009985.1_ASM998v1/GCF_000009985.1_ASM998v1_genomic.fna.gz
gunzip GCF_000009985.1_ASM998v1_genomic.fna.gz

echo ">chrom1" >> 2chrom.fasta
tail -n+2  GCF_000009985.1_ASM998v1_genomic.fna | head -n 250 >> 2chrom.fasta
echo ">chrom2" >> 2chrom.fasta
tail -n 250  GCF_000009985.1_ASM998v1_genomic.fna >> 2chrom.fasta
mv GCF_000009985.1_ASM998v1_genomic.fna AP017923.1.fasta

art_illumina -ss HS25 -sam -i 2chrom.fasta -p -l 100 -f 10 -m 200 -s 10 -o 2chrom

art_illumina -ss HS25 -i AP017923.1.fasta -p -l 150 -f 10 -m 200 -s 10 -o AP017923.1_reads



bwa index 2chrom.fasta

bwa mem -M  2chrom.fasta 2chrom1.fq 2chrom2.fq | samtools sort - > mapped.sam

samtools view -q 10 -b mapped.sam | samtools sort - > mapped.bam


# get_genomes.py -q AP007255.1 -o ./tests/references/

# art_illumina -ss HS25 -sam -i ./tests/references/AP007255.1.fasta -p -l 100 -f 20 -m 200 -s 10 -o tests/references/NC_002951_

#THE reference\n  get_genomes.py -q AP017922.1 -o ./tests/references/
# art_illumina -ss HS25 -sam -i ./tests/references/NC_002951.2.fasta -p -l 100 -f 20 -m 200 -s 10 -o tests/references/full_
# AP017923.1.fasta
# bwa mem -M  tests/references/AP017923.1.fasta tests/references/AP017923.1_reads1.fq tests/references/AP017923.1_reads2.fq | samtools sort - > mapped.sam
# bwa index tests/references/AP017923.1.fasta
# art_illumina -ss HS25 -i AP017923.1.fasta -p -l 150 -f 20 -m 200 -s 10 -o ./AP017923.1_reads
# samtools view -q 10 -b mapped.sam | samtools sort - > mapped.bam
# art_illumina -ss HS25 -i AP017923.1.fasta -p -l 150 -f 10 -m 200 -s 10 -o ./AP017923.1_reads
# split the AP017923.1 into two sequences to better test on multi chrom sequences
# art_illumina -ss HS25 -sam -i ./tests/references/2chrom.fasta -p -l 100 -f 10 -m 200 -s 10 -o tests/references/2chrom
