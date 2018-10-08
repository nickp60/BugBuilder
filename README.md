BugBuilder NRW

----------
Updates to the BugBuilder pipeline incorporating two major chanages:
1. Ability to submit --trusted-contigs to SPAdes
2. Ability to restart an analaysis after the (timeconsuming) assembly step

Additionally, various bug fixes were incorportated to make it work with more
recent versions of some of the programs.

==========

Hands-free microbial genome assembly

BugBuilder is a pipeline to facilitate assembly and annotation of database submission-ready,
draft-quality microbial genomes from high-throughput sequencing data. It utilises a range of
existing tools to produce an annotated genome assembly (optionally scaffolded against a reference
genome), outputting EMBL formatted records describing the assembled contigs, and AGP v2 files
describing the scaffolds. Required inputs are fastq format reads from either a fragmment of
paired-end sequencing run, and optionally a reference genome sequence in fasta format.

BugBuilder is implemented in Perl and is distributed under the Artistic License v2.0

Installation
============

Everything is now available to install via conda:

```
conda create --name bugbuilder abyss amos aragorn arrow barrnap biopython blast bwa cgview curl fastqc gmp gnuplot minced mummer mkl  parallel pilon prokka picard prodigal pyyaml samtools sis seqtk sickle-trim spades tabulate tbl2asn hmmer vcflib
source activate bugbuilder
```
then, clone this repo, and install:
```
python setup.py install
```
To set the configuration, run the following command:
```
BugBuilder --configure
```

and you are ready to go!  Have a try with the test data first:


```
BugBuilder --platform illumina -o output --assemblers spades --fastq1 ./tests/references/2chrom1.fq --fastq2 ./tests/references/2chrom2.fq --reference ./tests/references/2chrom.fasta  --scaffolder sis --finisher pilon
```

A virtual machine image preconfigured with the freely redistributable
prerequistite packages is available from
http://web.bioinformatics.ic.ac.uk/BugBuilder/BugBuilder_current.vdi

The latest version of the software can be downloaded by running:
<code>git clone git://github.com/jamesabbott/BugBuilder/</code>

Full installation and configuration instructions are available in the user
guide, which can be obtained from
http://www.imperial.ac.uk/bioinformatics-support-service/resources/software/BugBuilder/,
or if you have latex installed on your machine you can build the documentation
locally:

<code>
[jamesa@codon ~]$ cd BugBuilder/doc
[jamesa@codon ~]$ ./build.sh
</code>



conda create -n bug python==3.5 spades==3.9 fastqc sickle-trim
