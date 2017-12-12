#!/usr/bin/env python3
#-*- coding: utf-8 -*-

################################################################################
#
# $HeadURL: https://bss-srv4.bioinformatics.ic.ac.uk/svn/BugBuilder/trunk/bin/BugBuilder $
# $Author: jamesa $
# $Revision: 181 $
# $Date: 2016-03-13 16:00:41 +0000 (Sun, 13 Mar 2016) $
#
# This file is part of BugBuilder (https://github.com/jamesabbott/BugBuilder)
# and is distributed under the Artistic License 2.0 - see accompanying LICENSE
# file for details
#
################################################################################

__config_data__ = \
    """STATUS: COMPLETE
######################################################################
#
# BugBuilder configuration in YAML format
#
# This file defines the BugBuilder configuration. See the BugBuilder
# User Guide for details of the dependencies which need to be installed.
#
######################################################################
# tmp_dir specifies the location on the machine where working directories will be created
tmp_dir: null
# java specifies the java binary
java: null
# number of parallel threads to run
threads: 4

# Definition of assembly categories, and platforms
# These are used for automated assembler selection based on assesment of the
# provided reads. These should ideally not overlap or the choice of category
# may become a bit random
assembler_categories:
  - name: 'short_illumina'
    min_length: 25
    max_length: 100
    single_fastq: 'optional'
    paired_fastq: 'optional'
    platforms:
      - 'illumina'
    assemblers:
      - spades
      - abyss
    scaffolders:
      - mauve
      - SIS
      - sspace
  - name: 'long_illumina'
    min_length:  75
    max_length: 250
    single_fastq: 'optional'
    paired_fastq: 'optional'
    platforms:
      - 'illumina'
    assemblers:
      - spades
      - celera
    scaffolders:
      - mauve
      - SIS
      - sspace
  - name: '454_IonTorrent'
    min_length: 100
    max_length: 1000
    single_fastq: 'optional'
    paired_fastq: 'optional'
    platforms:
      - '454'
      - 'iontorrent'
    assemblers:
      - celera
    scaffolders:
      - mauve
      - SIS
  - name: 'long'
    min_length: 500
    max_length: 500000
    long_fastq: 'required'
    platforms:
      - 'PacBio'
      - 'MinION'
    assemblers:
      - PBcR
  - name: 'hybrid'
    min_length: 75
    max_length: 50000
    platforms:
      - hybrid
    paired_fastq: 'required'
    long_fastq: 'required'
    assemblers:
      - masurca
      - spades
    scaffolders:
      - mauve
  - name: 'de_fere'
    min_length: 75
    max_length: 50000
    platforms:
      - de_fere
    paired_fastq: 'required'
    de_fere_contigs: 'required'
    assemblers:
      - masurca
      - spades
    scaffolders:
      - SIS
      - mauve

#Assembler configuration
assemblers:
   - name: abyss
     create_dir: 1
     min_length: null
     max_length: null
     command_pe: run_abyss -C __TMPDIR__/abyss/  in='__FASTQ1__  __FASTQ2__' name=abyss
     command_se: run_abyss -C __TMPDIR__/abyss/  in=__FASTQ1__ name=abyss
     contig_output: __TMPDIR__/abyss/abyss-contigs.fa
     scaffold_output: __TMPDIR__/abyss/abyss-scaffolds.fa
     default_args: k=55
     insert_size_required: 0
     downsample_reads: 1
   - name: spades
     create_dir: 0
     min_length: null
     max_length: 303
     command_se: run_spades --memory __MEMORY__ -t __THREADS__ -s __FASTQ1__ -o __TMPDIR__/spades
     command_pe: run_spades --memory __MEMORY__ -t __THREADS__ -1 __FASTQ1__ -2 __FASTQ2__ -o __TMPDIR__/spades
     command_hybrid: run_spades --memory __MEMORY__ -t __THREADS__ -1 __FASTQ1__ -2 __FASTQ2__ --pacbio __LONGFASTQ__ -o __TMPDIR__/spades
     command_de_fere: run_spades --memory __MEMORY__ -t __THREADS__ -1 __FASTQ1__ -2 __FASTQ2__ --trusted-contigs __DE_FERE_CONTIGS__ -o __TMPDIR__/spades
     contig_output: __TMPDIR__/spades/contigs.fasta
     scaffold_output: __TMPDIR__/spades/scaffolds.fasta
     default_args: -t __THREADS__ --careful
     insert_size_required: 0
     downsample_reads: 1
   - name: riboSeed
     create_dir: 0
     min_length: null
     max_length: null
     command_se: ribo run -r __REFERENCE__ -S1 __FASTQ1__ -o __TMPDIR__/riboSeed/ --memory __MEMORY__ --cores __THREADS__
     command_pe: ribo run -r __REFERENCE__ -F __FASTQ1__ --R __FASTQ2__ -o __TMPDIR__/riboSeed/ --memory __MEMORY__ --cores __THREADS__
     contig_output: __TMPDIR__/seed/final_de_fere_novo_assebmly/contigs.fasta
     scaffold_output: __TMPDIR__/seed/final_de_fere_novo_assebmly/scaffolds.fasta
     default_args: null
     insert_size_required: 0
     downsample_reads: 1
   - name: celera
     create_dir: 1
     min_length: 75
     max_length: null
     command_se: run_celera --fastq1 __FASTQ1__ --tmpdir __TMPDIR__ --category __CATEGORY__ --encoding __ENCODING__ --genome_size __GENOME_SIZE__
     command_pe: run_celera --fastq1 __FASTQ1__ --fastq2 --tmpdir __TMPDIR__ --category __CATEGORY__ --encoding __ENCODING__ --genome_size __GENOME_SIZE__
     contig_output: __TMPDIR__/celera/output/9-terminator/BugBuilder.ctg.fasta
     scaffold_output: __TMPDIR__/celera/output/9-terminator/BugBuilder.scf.fasta
     downsample_reads: 0
     insert_size_required: 0
     # masurca works best with untrimmed reads, so use __ORIG_FASTQ1__ nad __ORIG_FASTQ2__
   - name: masurca
     create_dir: 1
     min_length: null
     max_length: null
     command_pe: run_masurca --fastq1 __ORIG_FASTQ1__ --fastq2 __ORIG_FASTQ2__ --tmpdir __TMPDIR__ --category __CATEGORY__ --insert_size __INSSIZE__ --insert_stddev __INSSD__
     command_hybrid: run_masurca --fastq1 __ORIG_FASTQ1__ --fastq2 __ORIG_FASTQ2__ --longfastq __LONGFASTQ__ --tmpdir __TMPDIR__ --category __CATEGORY__ --insert_size __INSSIZE__ --insert_stddev __INSSD__
     contig_output: __TMPDIR__/masurca/contigs.fasta
     scaffold_output: __TMPDIR__/masurca/scaffolds.fasta
     default_args: --threads __THREADS__
     downsample_reads: 0
     insert_size_required: 1

scaffolders:
   - name: SIS
     ref_needed: single
     linkage_evidence: align_genus
     command: run_sis --reference __REFERENCE__ --contigs __CONTIGS__ --tmpdir __TMPDIR__ --scaff_dir __SCAFFDIR__
     scaffold_output: scaffolds.fasta
     unscaffolded_output: unplaced_contigs.fasta
     default_args: null
     create_dir: 1
     priority: 2
     multiple_refs: 0
   - name: ragout
     ref_needed: multiple
     linkage_evidence: align_genus
     command: run_ragout --reference __REFERENCE__ --contigs __CONTIGS__ --tmpdir __TMPDIR__ --scaff_dir __SCAFFDIR__
     scaffold_output: scaffolds.fasta
     unscaffolded_output: unplaced_contigs.fasta
     default_args: null
     create_dir: 1
     priority: 2
   - name: mauve
     ref_needed: single
     linkage_evidence: align_genus
     command: run_mauve --reference __REFERENCE__ --run __RUN__ --contigs __CONTIGS__ --tmpdir __TMPDIR__ --scaff_dir __SCAFFDIR__
     default_args: null
     create_dir: 1
     priority: 1
     scaffold_output: scaffolds.fasta
   - name: sspace
     ref_needed: null
     linkage_evidence: paired-ends
     command: run_sspace --tmpdir __TMPDIR__ --scaff_dir __SCAFFDIR__ --contigs __CONTIGS__ --insert_size __INSSIZE__ --insert_sd __INSSD__
     scaffold_output: BugBuilder.scaffolds.fasta
     default_args: null
     create_dir: 1
     priority: 3

merge_tools:
   - name: gfinisher
     command: run_gfinisher --tmpdir __TMPDIR__  --assembler __ASSEMB1__ --assembler __ASSEMB2__ --reference __REFERENCE__
     contig_output: renamed.fasta
     create_dir: 1
     priority: 1
     allow_scaffolding: 1
   - name: minimus
     command: run_minimus --tmpdir __TMPDIR__  --assembler __ASSEMB1__ --assembler __ASSEMB2__
     contig_output: renumbered.fasta
     create_dir: 1
     priority: 2
     allow_scaffolding: 1

finishers:
   - name: gapfiller
     command: run_gapfiller --tmpdir __TMPDIR__ --insert_size __INSSIZE__ --insert_sd __INSSD__ --threads __THREADS__
     create_dir: 1
     ref_required: 0
     paired_reads: 1
     output_scaffolds: pilon.fasta
     needs: variable_gaps
     priority: 2
       - SIS
       - SSPADES
   - name: abyss-sealer
     command: run_abyss-sealer --tmpdir __TMPDIR__ --encoding __ENCODING__ --threads __THREADS__
     create_dir: 1
     output_scaffolds: sealed_scaffold.fa
     needs: fixed_gaps
     ref_required: 0
     paired_reads: 1
     priority: 3
   - name: pilon
     command: run_pilon --tmpdir __TMPDIR__ --threads __THREADS__
     output_scaffolds: pilon.fasta
     create_dir: 1
     needs: fixed_gaps
     ref_required: 0
     paired_reads: 1
     priority: 1
   - name: fgap
     command: run_pilon --tmpdir __TMPDIR__ --threads __THREADS__
     output_scaffolds: fgap.final.fasta
     create_dir: 1
     needs: fixed_gaps
     ref_required: 0
     paired_reads: 0
     priority: 1

varcallers:
   - name: pilon
     command: run_pilon_var --tmpdir __TMPDIR__ --threads __THREADS__
     ref_required: 1
     create_dir: 1
     priority: 1

"""
