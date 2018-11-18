# test data
We have a tricky problem --  we want the test data to be deterministic, but assemblers are not.  Here we have committed the results of the following analyses, with the understanding that the generating proceses are dynamic and may change in the future

## creating a test genome and reads
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/985/GCF_000009985.1_ASM998v1/GCF_000009985.1_ASM998v1_genomic.fna.gz
gunzip GCF_000009985.1_ASM998v1_genomic.fna.gz

echo ">chrom1" >> 2chrom.fasta
tail -n+2  GCF_000009985.1_ASM998v1_genomic.fna | head -n 250 >> 2chrom.fasta
echo ">chrom2" >> 2chrom.fasta
tail -n 250  GCF_000009985.1_ASM998v1_genomic.fna >> 2chrom.fasta
mv GCF_000009985.1_ASM998v1_genomic.fna AP017923.1.fasta

art_illumina -ss HS25 -sam -i 2chrom.fasta -p -l 100 -f 10 -m 200 -s 10 -o 2chrom

# art_illumina -ss HS25 -i AP017923.1.fasta -p -l 150 -f 10 -m 200 -s 10 -o AP017923.1_reads
```

We delete the .aln files.
## generating mapping files
```
bwa index 2chrom.fasta

bwa mem -M  2chrom.fasta 2chrom1.fq 2chrom2.fq | samtools sort - > mapped.sam

samtools view -q 10 -b mapped.sam | samtools sort - > mapped.bam
```

Note that we also deleted the index files (.amp .ann .bwt .pac .sa) for the reference
# generating a little assembly
```
# k is small so we have lots of contigs for the example
spades.py -o ./assembly/ -1 2chrom1.fq -2 2chrom2.fq --careful -k 17
```
We then deleted some of the intermediate files from the assembly to keep the repo small.

```
./tests/references/
├── 2chrom.fasta
├── 2chrom.sam
├── AP017923.1.fasta
├── 2chrom1.fq
├── 2chrom2.fq
├── README.md
├── already_assembled_test
│   ├── contigs.fasta
│   └── scaffolds.fasta
├── assembly
│   ├── K17
│   ├── assembly_graph.fastg
│   ├── assembly_graph_with_scaffolds.gfa
│   ├── before_rr.fasta
│   ├── contigs.fasta
│   ├── contigs.paths
│   ├── dataset.info
│   ├── input_dataset.yaml
│   ├── params.txt
│   ├── scaffolds.fasta
│   ├── scaffolds.paths
│   └── spades.log
├── configs
│   ├── broken_config.yaml
│   ├── empty_config.yaml
│   ├── semicomplete_config.yaml
│   └── static_config.yaml
├── contigs_split_ori.fasta
├── contigs_to_scaffold.fasta
├── distant_contigs.fasta
├── encoding
│   ├── illumina13.fq
│   ├── illumina15.fq
│   ├── illumina18_sanger.fq
│   ├── sanger.fq
│   └── solexa.fq
├── mapped.bam
├── mapped.sam
├── needs_renaming.fq
├── origin
│   ├── delta-filter.log
│   ├── nucmer.log
│   ├── ori.coords
│   ├── ori.delta
│   ├── ori.filter
│   └── show-coords.log
├── picard_insert.txt
├── renamed_ref.fq
├── scaffs.fasta
├── sickle.log
└── sickle_bad.log

6 directories, 44 files
```
`contigs_split_ori.fasta` has the results of the origin splitting script, checked with Mauve
The `encoding` dir has mini fastqs manually created to conform to the fastq quality encodings of the desired file types.

`insert_stats_new.txt` and `picard_insert_old.txt` are the results of the picard tool to check the insert size.  This syntax changed in a recent version, so we included both version to test the parsing method.

`needs_renaming.fq` and `renamed.fq` show the before and after of the renaming mothod, which looks for reads ending in -1 or -2 changing them to /1 and /2

`sickle.log` is , well, a log from a successful run sickle. `sickle_bad.log` shows the results of a run where too much trimming occured, which should result in a warning.
