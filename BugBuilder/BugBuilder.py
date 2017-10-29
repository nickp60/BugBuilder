#!/usr/bin/env python3

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
---
# tmp_dir specifies the location on the machine where working directories will be created
tmp_dir: <TMPL_VAR NAME=TMP_DIR>
# java specifies the java binary
java: <TMPL_VAR NAME=JAVA>
# number of parallel threads to run
threads: <TMPL_VAR NAME=THREADS>
<TMPL_VAR NAME=INSTALL_PATHS>

<TMPL_VAR NAME=APPEND_PATH>

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
    max_length: 50000
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
     max_length: 200
     command_pe: __BUGBUILDER_BIN__/run_abyss --tmpdir __TMPDIR__ --fastq1 __FASTQ1__ --fastq2 __FASTQ2__ --read_length __READ_LENGTH__
     contig_output: __TMPDIR__/abyss/abyss-contigs.fa
     scaffold_output: __TMPDIR__/abyss/abyss-scaffolds.fa
     downsample_reads: 1
   - name: spades
     create_dir: 0
     max_length: 303
     command_se: __ASMDIR__/spades.py -s __FASTQ1__ -o __TMPDIR__/spades
     command_pe: __ASMDIR__/spades.py -1 __FASTQ1__ -2 __FASTQ2__ -o __TMPDIR__/spades
     command_hybrid: __ASMDIR__/spades.py -1 __FASTQ1__ -2 __FASTQ2__ --pacbio __LONGFASTQ__ -o __TMPDIR__/spades
     command_de_fere: __ASMDIR__/spades.py -1 __FASTQ1__ -2 __FASTQ2__ --trusted-contigs __DE_FERE_CONTIGS__ -o __TMPDIR__/spades
     contig_output: __TMPDIR__/spades/contigs.fasta
     scaffold_output: __TMPDIR__/spades/scaffolds.fasta
     default_args: -t __THREADS__ --careful
     downsample_reads: 1
   - name: celera
     create_dir: 1
     min_length: 75
     command_se: __BUGBUILDER_BIN__/run_celera --fastq1 __FASTQ1__ --tmpdir __TMPDIR__ --category __CATEGORY__ --encoding __ENCODING__ --genome_size __GENOME_SIZE__
     command_pe: __BUGBUILDER_BIN__/run_celera --fastq1 __FASTQ1__ --fastq2 --tmpdir __TMPDIR__ --category __CATEGORY__ --encoding __ENCODING__ --genome_size __GENOME_SIZE__
     contig_output: __TMPDIR__/celera/output/9-terminator/BugBuilder.ctg.fasta
     scaffold_output: __TMPDIR__/celera/output/9-terminator/BugBuilder.scf.fasta
     downsample_reads: 0
   - name: PBcR
     create_dir: 1
     min_length: 500
     command_se: __BUGBUILDER_BIN__/run_PBcR --fastq __LONGFASTQ__ --tmpdir __TMPDIR__ --genome_size __GENOME_SIZE__ --platform __PLATFORM__
     contig_output: __TMPDIR__/PBcR/BugBuilder/9-terminator/asm.ctg.fasta
     scaffold_output: __TMPDIR__/PBcR/BugBuilder/9-terminator/asm.scf.fasta
     downsample_reads: 0
     # masurca works best with untrimmed reads, so use __ORIG_FASTQ1__ nad __ORIG_FASTQ2__
   - name: masurca
     create_dir: 1
     command_pe: __BUGBUILDER_BIN__/run_masurca --fastq1 __ORIG_FASTQ1__ --fastq2 __ORIG_FASTQ2__ --tmpdir __TMPDIR__ --category __CATEGORY__ --insert_size __INSSIZE__ --insert_stddev __INSSD__
     command_hybrid: __BUGBUILDER_BIN__/run_masurca --fastq1 __ORIG_FASTQ1__ --fastq2 __ORIG_FASTQ2__ --longfastq __LONGFASTQ__ --tmpdir __TMPDIR__ --category __CATEGORY__ --insert_size __INSSIZE__ --insert_stddev __INSSD__
     contig_output: __TMPDIR__/masurca/contigs.fasta
     scaffold_output: __TMPDIR__/masurca/scaffolds.fasta
     default_args: --threads __THREADS__
     downsample_reads: 0
     insert_size_required: 1

scaffolders:
   - name: SIS
     linkage_evidence: align_genus
     command: __BUGBUILDER_BIN__/run_sis --reference __REFERENCE__ --contigs __CONTIGS__ --tmpdir __TMPDIR__ --scaff_dir __SCAFFDIR__
     scaffold_output: scaffolds.fasta
     unscaffolded_output: unplaced_contigs.fasta
     create_dir: 1
     priority: 2
   - name: mauve
     linkage_evidence: align_genus
     command: __BUGBUILDER_BIN__/run_mauve --reference __REFERENCE__ --run __RUN__ --contigs __CONTIGS__ --tmpdir __TMPDIR__ --scaff_dir __SCAFFDIR__
     create_dir: 1
     priority: 1
     scaffold_output: scaffolds.fasta
   - name: sspace
     linkage_evidence: paired-ends
     command: __BUGBUILDER_BIN__/run_sspace --tmpdir __TMPDIR__ --scaff_dir __SCAFFDIR__ --contigs __CONTIGS__ --insert_size __INSSIZE__ --insert_sd __INSSD__
     scaffold_output: BugBuilder.scaffolds.fasta
     create_dir: 1
     priority: 3

merge_tools:
   - name: gfinisher
     command: __BUGBUILDER_BIN__/run_gfinisher --tmpdir __TMPDIR__  --assembler __ASSEMB1__ --assembler __ASSEMB2__ --reference __REFERENCE__
     contig_output: renamed.fasta
     create_dir: 1
     priority: 1
     allow_scaffolding: 1
   - name: minimus
     command: __BUGBUILDER_BIN__/run_minimus --tmpdir __TMPDIR__  --assembler __ASSEMB1__ --assembler __ASSEMB2__
     contig_output: renumbered.fasta
     create_dir: 1
     priority: 2
     allow_scaffolding: 1

finishers:
   - name: gapfiller
     command: __BUGBUILDER_BIN__/run_gapfiller --tmpdir __TMPDIR__ --insert_size __INSSIZE__ --insert_sd __INSSD__ --threads __THREADS__
     create_dir: 1
     ref_required: 0
     paired_reads: 1
     priority: 2
   - name: abyss-sealer
     command: __BUGBUILDER_BIN__/run_abyss-sealer --tmpdir __TMPDIR__ --encoding __ENCODING__ --threads __THREADS__
     create_dir: 1
     ref_required: 0
     paired_reads: 1
     priority: 3
   - name: pilon
     command: __BUGBUILDER_BIN__/run_pilon --tmpdir __TMPDIR__ --threads __THREADS__
     create_dir: 1
     ref_required: 0
     paired_reads: 1
     priority: 1

varcallers:
   - name: pilon
     command: __BUGBUILDER_BIN__/run_pilon_var --tmpdir __TMPDIR__ --threads __THREADS__
     ref_required: 1
     create_dir: 1
     priority: 1
"""


"""
=item B<scaffolder>: Scaffolder to run

=item B<scaffolder-args>: Any additional arguments to pass to the scaffolder. Overides the setting
of the 'default_args' setting in the scaffolder configuration

=item B<merge-method>: Assembly merging tool to use

=item B<finisher>: Method for assembly finishing/gap closure

=item B<varcall>: Method for variant calling

=item B<insert-size>: Size of insert in paired-read library. This will be determined empircally if a reference
genome sequence is provided, so only needs specifying when assembling paired-read sequences for
which no reference genome is available.

=item B<insert-stddev>: standard deviation of insert in paired-read library. This will be
determined empircally if a reference genome sequence is provided, so only needs specifying when
assembling paired-read sequences for which no reference genome is available.

=item B<genome-size>: Approximate genome size. Required for PacBio/MinION assemblies using PBcR

=item B<genus>: Genus of organism sequenced (i.e. Streptococcus). Included in resutling EMBL file,
and passed to Prokka during annotation stage.

=item B<species>: Specific name of species if known (i.e. pyogenes). Included in resutling EMBL
file, and passed to Prokka during annotation stage.

=item B<strain>: Name of strain used for inclusion in annotation results.

=item B<mode>: Mode to run in - valid modes are 'submission' (default) or 'draft'

=item B<locustag>: Locustag argument to pass to Prokka. Used to customise
locus_tag in generated EMBL records.

=item B<centre>: Sequence centre argument to pass to Prokka. Used to customise
locus_tag in generated EMBL records.

=item B<[no]-fastqc>: Determine whether to run fastqc: Default: on

=item B<[no]-trim>: Determine wheter to quality trim reads. Default: on

=item B<trim-qv>: Quality threshold for trimming reads. Default: 20

=item B<trim-length>: Min. length of read to retain following trimming. Default: 50 (25 for reads <50bp)

=item B<[no-]split-origin>: Determine wether to attempt to split assembly around
the origin or not. Should the assembly not be scaffolded using a reference
sequence, or where the reference sequence is in a significant number of
contigs, then starting the sequence at the origin makes little sense. Default: on

=item B<[no-]gap-fill>: Determe wether to run GapFiller to close scaffold gaps.
Assemblies with a large number of scaffold gaps can result in the gap filling
stage taking a significant amount of time. Default: on

=item B<keepall>: Return full working directory with intermediate files, rather the just returning
the results (default off)

=item B<[no-]gap-seal>: Determine whether to run Abyss sealer to close scaffold gaps. Default: on

=item B<help>: display short help text

=item B<man>: display full help text

=item B<threads>: Number of threads to use for multi-threaded applications. Default: 1

=item B<out-dir>: Path for the output directory. If not provided, will be generated under the current
working directory.

=item B<tmp-dir>: Path for the temporary working directory.
If not provided, unique temporary directory will be created under the temp_dir path set in
a config file.

=back

=head1 REPORTING BUGS

Please report any bugs/issues via github:
https://github.com/jamesabbott/BugBuilder/issues/new.

All bug reports should include the output of the 'check_config.pl' script,
which reports on the installed software packages and versions

=head1 AUTHOR - James Abbott

Email j.abbott@imperial.ac.uk

=cut

use warnings;
use strict;

"""

import argparse
import os
import yaml
import sys
import arrow
import shutil
import logging
import statistics
import pkg_resources

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# def parse_config():
#     pass

def parse_available_platforms(config):
    platform_list = ["short_illumina", "long_illumina", "de_fere"]
    return platform_list


def parse_available_assemblers():
    assembler_list = ["abyss", "spades", "mascura"]
    return assembler_list


def parse_available_scaffolders():
    scaffolder_list = ["sis", "mauve"]
    return scaffolder_list


def parse_available_mergers():
    mergers_list = ["gfinisher"]
    return mergers


def parse_available_finishing():
    finishers = ["gapfiller"]
    return finishers


def parse_available_varcaller():
    finishers = ["pilon"]
    return finishers


def configure(config_path):
    programs = [
        'fastqc', 'sickle', 'seqtk', 'samtools', 'picard' , 'sis', 'mauve',
        'R', 'barrnap', 'prokka', 'aragorn', 'prodigal', 'hmmer3', 'rnammer',
        'mummer', 'infernal', 'blast', 'bwa', 'tbl2asn', 'abyss', 'spades',
        'celera','gapfiller', 'sspace', 'asn2gb', 'amos' , 'masurca', 'gfinisher',
        'pilon', 'vcflib', 'cgview']
    program_dict = dict((k, None) for k in programs)
    for prog in programs:
        if shutil.which(prog):
            program_dict[prog] = shutil.which(prog)
    with open(config_path, 'w') as outfile:
        outfile.write("STATUS: COMPLETE\n")
    with open(config_path, 'a') as outfile:
        yaml.dump(program_dict, outfile, default_flow_style=False)


def return_config(config_path):
    config = parse_config(config_path)
    if config.STATUS != "COMPLETE":
        configure(config_path)
        config = parse_config(config_path)
    return config


def get_args():  # pragma: no cover
    """
    """
    parser = argparse.ArgumentParser(
        prog="BugBuilder",
        description="Automated pipeline for assembly of draft quality " +
        "bacterial genomes with reference guided scaffolding and annotation." +
        "Please see accompanying userguide for full documentation",
        add_help=False)  # to allow for custom help
    parser.add_argument("contigs", action="store",
                        help="either a (multi)fasta or a directory " +
                        "containing one or more chromosomal " +
                        "sequences in fasta format")

    # taking a hint from http://stackoverflow.com/questions/24180527
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("--platform", dest='platform', action="store",
                               help="Sequencing platform  used i.e. illumina, 454, iontorrent",
                               choices=parse_available_platforms,
                               type=str, required=True)
    requiredNamed.add_argument("-o", "--output", dest='output', action="store",
                               help="output directory; " +
                               "default: %(default)s", default=os.getcwd(),
                               type=str, required=True)

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--fastq1", dest='fastq1', action="store",
                               help="Path to first read of paired library, or fragment library",
                               type=str, required="--long-fastq" not in sys.argv)
    optional.add_argument("--fastq2", dest='fastq2', action="store",
                               help="Path to second read of paired library",
                               type=str)
    optional.add_argument("--de_fere_contigs", dest='de_fere_contigs',
                               action="store",
                               help="contigs to be used in de fere novo assembly!",
                               type=str)
    optional.add_argument("--long_fastq", dest='long_fastq', action="store",
                               help="Path to fastq file from long-read sequencer",
                               type=str, required="--fastq1" not in sys.argv)
    optional.add_argument("--reference", dest='reference',
                          action="store",
                          help="Path to fasta formatted reference sequence",
                          type=str)
    optional.add_argument("--prefix", dest='prefix',
                          action="store",
                          help="Prefix to use for output file naming",
                          type=str)
    optional.add_argument("--assemblers", dest='assemblers', action="store",
                          help="Assembler(s) to run - may be specified twice" +
                          ", in which case the two assemblers will be run " +
                          "in parallel and the results merged using minimus." +
                          " If no assembler is specified, BugBuilder will " +
                          "try to select an appropriate assembler " +
                          "automatically",
                          choices=parse_available_assemblers,
                          nargs="*",
                          type=str)
    optional.add_argument("--assembler-args", dest='assembler-args',
                          action="store",
                          help="Any additional arguments to pass to the " +
                          "assembler. Default values are set in the " +
                          "'default_args' attribute of the configuration " +
                          "file. If running multiple assemblers, " +
                          "assembler_args should be specified twice, once " +
                          "for each assemler, in the same order than the " +
                          "assemblers are specified.",
                          nargs="*",
                          type=str)
    optional.add_argument("--scaffolder", dest='scaffolder', action="store",
                          help="scaffolder to use",
                          choices=parse_available_scaffolders,
                          type=str)
    optional.add_argument("--assemlber-args", dest='scaffolder-args',
                          action="store",
                          help="args to pass to the scaffolder, in single quotes",
                          nargs="*",
                          type=str)
    optional.add_argument("--merge-method", dest='merge_method', action="store",
                          help="merge method to use",
                          choices=parse_available_mergers,
                          type=str)
    optional.add_argument("--finisher", dest='finisher', action="store",
                          help="finisher to use",
                          choices=parse_available_finisher,
                          nargs="*",
                          type=str)
    optional.add_argument("--varcaller", dest='varcaller', action="store",
                          help="varcaller to use",
                          choices=parse_available_varcaller,
                          nargs="*",
                          type=str)
    optional.add_argument("--insert-size", dest='insert_size', action="store",
                          help="insert size (bp)",
                          type=int)
    optional.add_argument("--insert-stddev", dest='insert_stddev', action="store",
                          help="insert size standard deviation ",
                          type=int)
    optional.add_argument("--genome-size", dest='genome_size', action="store",
                          help="size of the genome ",
                          type=int, default=0) # 0 is better than None for addition :)
    optional.add_argument("--downsample", dest='downsample', action="store_true",
                          help="size of the genome ", default=False)
    optional.add_argument("--species", dest='species', action="store",
                          help="species of your bug", default="unknown_species",
                          type=str)
    optional.add_argument("--genus", dest='genus', action="store",
                          help="genus of your bug", default="unknown_genus",
                          type=str)
    optional.add_argument("--strain", dest='strain', action="store",
                          help="strain of your bug", defualt="unknown_strain",
                          type=str)
    optional.add_argument("--locustag", dest='locustag', action="store",
                          help="locus tag prefix for prokka",
                          type=str)
    optional.add_argument("--centre", dest='centre', action="store",
                          help="Sequence centre argument to pass to Prokka. " +
                          "Used to customise locus_tag in generated EMBL records.",
                          type=str)
    optional.add_argument("--mode", dest='mode', action="store",
                          help="Mode to run in",
                          choices=["submission", "draft"], default="submission",
                          type=str)
    optional.add_argument("--skip-fastqc", dest='skip_fastqc', action="store_true",
                          help="size of the genome ", default=False)
    optional.add_argument("--skip-trim", dest='skip_time', action="store_true",
                          help="Quality threshold for trimming reads ", default=False)
    optional.add_argument("--trim-qv", dest='trim_qv', action="store",
                          help="quality threshold to trim  ", default=20,
                          type=int)
    optional.add_argument("--trim-length", dest='trim_length', action="store",
                          help="min read klength to retain", default=50,
                          type=int)
    optional.add_argument("--skip-split-origin", dest='skip_split_origin',
                          action="store_true",
                          help="split at origin ", default=False)
    optional.add_argument("--keepall", dest='keepall',
                          action="store_true",
                          help="keep all intermediate files ", default=False)
    optional.add_argument("--threads", dest='threads', action="store",
                          help="threads  ",
                          type=int)
    optional.add_argument("--out-dir", dest='out_dir', action="store",
                          help="dir for results",
                          type=str)
    optional.add_argument("--tmp-dir", dest='tmp_dir', action="store",
                          help="dir for results",
                          type=str)
    optional.add_argument("--already_assembled", dest='already_assembled',
                          action="store_true",
                          help="use existing assembly from --scratchdir ",
                          default=False)
    optional.add_argument("--scratchdir", dest='scratchdir', action="store",
                          help="dir for results",
                          type=str)
    optional.add_argument("-v", "--verbosity", dest='verbosity',
                          action="store",
                          default=2, type=int, choices=[1, 2, 3, 4, 5],
                          help="Logger writes debug to file in output dir; " +
                          "this sets verbosity level sent to stderr. " +
                          " 1 = debug(), 2 = info(), 3 = warning(), " +
                          "4 = error() and 5 = critical(); " +
                          "default: %(default)s")
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")
    # args = parser.parse_args(sys.argv[2:])
    args = parser.parse_args()
    return args


def match_assembler_args(assembler, assembler_args):
    """these will both be lists coming off of argparse
    """
    assert len(assembler) == len(assembler_args), \
        "length of assemblers must equal length of assembler args"
    assembler_list = []
    for i, v in enumerate(assembler):
        assembler_list.append([v, assembler_args[i]])
    return assembler_list


def check_file_paths(args):
    for f in [args.fastq1, args.fastq2, args.reference]:
        if f is not None and not os.path.exists(f):
            raise ValueError("file %s does not exist!" % f)
    if args.fastq1 is not None and args.fastq1 == args.fastq2:
        raise valueError("fastq1 and fastq2 are the same!")
    if args.mode == "draft" and args.reference is None:
        raise ValueError("Draft mode requiresa reference!")
    if args.reference is not None:
         if not os.path.exists(args.reference):
             raise ValueError(
                 "reference sequence %s does not exist!" % args.reference )
    if args.merge-method is None and len(args.assemblers) > 1:
        raise ValueError("Must provide merge method if using multiple assemblers")

def setup_tempdir(args, output_root):
    """  setup_tempdir creates a temporary directory and copies
    over the relevent files

    requried params: $ (path to fastq1)
                  $ (path to fastq2)
                  $ (path to longread fastq file)
                  $ (path to fastafile)
    optional params: $ (platform)
                  $ (preset value for tempdir)

    returns        : $ (path to tempdir)

    """
    if args.tmpdir is None:
        scratch_dir = os.path.join(output_root, "temp_dir")
    else:
        scratch_dir = args.tmpdir;

    print("Created working directory: %s" % scratch_dir)
    try:
        os.makedirs(scratch_dir, exist_ok=False)
    except OSError:
        print("Temp directory already exists; exiting...")
        sys.exit(1)
    # copy to temp dir
    for f in [args.fastq1, args.fastq2, args.long_fastq,
              args.de_fere_contigs, args.reference]:
        if os.path.exists(f):
            newpath = os.path.join(scratch_dir, os.path.basename(f))
            shutil.copyfile(f, newpath)
            # does this even work?
            f = newpath

    # implement this later
    # # we've seen some miseq fastq files have -1/-2 rather that /1 /2 pair ids which cause
    # # problems with older software which doesn't expect this
    # my @fastqs = ( basename($fastq1) ) if ($fastq1);
    # push @fastqs, basename($fastq2) if ($fastq2);
    # foreach my $file (@fastqs) {
    #     open FASTQ, "$tmpdir/$file" or die "Error opening $tmpdir/$file:$!";
    #     open NEW, ">$tmpdir/$file.new"
    #       or die "Error opening $tmpdir/$file.new: $!";
    #     while ( my $line = <FASTQ> ) {
    #         $line =~ s/^([@\+][0-9A-Z\:\-]+)-([12])$/$1\/$2/;
    #         print NEW $line;
    #     }
    #     close NEW;
    #     close FASTQ;
    #     move( "$tmpdir/$file.new", "$tmpdir/$file" )
    #       or die "Error copying $tmpdir/$file.new -> $tmpdir/$file: $!";
    # }


    # Different tools have differing interpretations of various fasta ids,
    # so, rewrite the reference fasta file to ensure these just contain the id.
    # accepted formats are for ENA and Genbank format fasta headers, as well as plain IDs
    if args.reference:
        new_reference = os.path.join(scratch_dir, os.path.basename(args.reference))
        with open(args.reference, "r") as inf:
            for rec in SeqIO.parse(inf, "fasta"):
                with open(new_reference, "a") as outd:
                    SeqIO.write(rec, writepath, 'fasta')
        args.reference = new_reference
    return scratch_dir;

def return_open_fun(f):
    if os.path.splitext(f)[-1] in ['.gz', '.gzip']:
        open_fun = gzip.open
    else:
        open_fun = open
    return open_fun


def get_read_len_from_fastq(fastq1, logger=None):
    """from LP; return total read length and count in fastq1 file from all reads
    """
    lengths = []
    open_fun = return_open_fun(fastq1)
    with open_fun(fastq1, "rt") as file_handle:
        data = SeqIO.parse(file_handle, "fastq")
        for read in data:
            lengths.append(len(read))
    if logger:
        logger.info(str("From the first {0} reads in {1}, " +
                        "mean length is {2}").format(N,
                                                     os.path.basename(fastq1),
                                                     float(tot / count)))
    file_handle.close()
    return (lengths)


def id_fastq_encoding(tempdir):
    """
    Identifies fastq quality encoding type, based upon the information
    at http://en.wikipedia.org/wiki/FASTQ_format:

    SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
    ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
    ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
    .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
    LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
    !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
    |                         |    |        |                              |                     |
    33                        59   64       73                            104                   126

    S - Sanger        Phred+33,  raw reads typically (0, 40)
    X - Solexa        Solexa+64, raw reads typically (-5, 40)
    I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
    J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
    with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
    L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)

    required params: $ (tmpdir)

    returns:         $ (encoding)
    """
    min_val = 104;
    max_val = 0;
    if os.path.exists(args.long_fastq):
        target = args.long_fastq
    else:
        target = args.fastq1
    open_fun = return_open_fun(target)
    with open_fun(target, "r") as inf:
        counter = 3
        for i, line in enumerate(inf):
            if i <= 1000:
                break
            if counter == 0:
                quals = [ord(x) for x in line.split()]
                tmp_min =  min(quals)
                tmp_max =  max(quals)
                min_val = tmp_min if tmp_min < min_val else min_val
                max_val = tmp_max if tmp_max > max_val else max_val
            else:
                counter = counter - 1
                continue
    if min_val <= 59 and max_val <= 74:
        encoding = 'sanger'
    elif min_val > 59 and min_val < 64 and max_val < 104:
        encoding = 'solexa'
    elif min_val > 64:
        encoding = 'illumina'
    else:
        print("assumina illumina encoding")
        encoding = 'illumina'
    return encoding


def assess_reads(args, conf, platform, tempdir, logger=None):
    """
    Determine some characteristics of the provided reads. Need to determine
    the mean read length, wheterh they are variable length or not, and the
    quality encoding.

    We'll also look at the predicted coverage if a reference genome is available

    required params: $ (tmpdir)

    optional params: $ (platform)

    returns: $ (assembler libtype)
         $ (quality encoding)
         $ (mean length)
         $ (stddev of read length)
         $ (projected genome coverage)
         $ (number of reads)
    my ( @mean, @stddev, $long_mean, $long_stddev, $tot_length, $long_tot_length );
    """
    logger.info("Checking reads...");
    means = []
    stddevs = []
    long_mean = 0
    long_stddev = 0
    tot_length = 0
    long_tot_length = 0
    # get the lengths of all the short reads
    type_list = ["short", "short", "long"]
    for i, read in enumerate([args.fastq1, args.fastq2, args.long_fastq]):
        if not os.path.exists(read):
            next
        lengths = get_read_len_from_fastq(read, logger=logger)
        if type_list[i] == "long":
            tot_long_length = tot_long_length + sum(lengths)
        else:
            tot_length = tot_length + sum(lengths)
        mean = statistics.mean(lengths)
        # gid rid of superlong 454 reads that mess with stats
        lengths = [x for x in lengths if x < (mean * 3)]
        if type_list[i] == "long":
            long_count = len(lengths)
            long_mean = statistics.mean(lengths)
            long_stddev = statistics.stddev(lengths)
        else:
            counts.append(len(lengths))
            means.append(statistics.mean(lengths))
            stddevs.append(statistics.stddev(lengths))

    # ensure paired reads have same numebr of reads
    if len(counts) > 1:
        if counts[0] != counts[1]:
            raise ValueError("Paired reads must have same number of reads!")
    #determin mean read length and stddev....
    mean = statistics.mean(means)
    stddev = statistics.mean(stddevs)
    # set type of library
    # TODO type is a bad variable name in python
    lib_type = None
    for category in config.assembler_categories:
        if platform is not None:
            if platform == category.platform:
                lib_type = category.name
        if lib_type == None:  # if that failed, try setting according to read length
            if mean > category.min_length and mean <= category.max_length:
                lib_type = category.name
            elif long_mean != 0:
                if long_mean > category.min_length and long_mean <= category.max_length:
                    lib_type = category.name
    if lib_type is None:
        raise ValueError(
            "Cound not find appropriate platform based on read length!")

    coverage = 0
    long_coverage = 0

    # not making a new one here
    if os.path.exists(args.reference):
        with open(args.reference, "r") as inf:
            for rec in SeqIO.parse(inf, "fasta"):
                args.genome_size = args.genome_size + len(rec.seq)
    try:
        coverage = tot_length / args.genome_size
        long_coverage = long_tot_length / args.genome_size
    except:
        logger.error("error calculating coverage!")
        sys.exit(1)
    encoding = id_fastq_encoding(tmpdir)
    return Namespace(lib_type=lib_type, encoding=encoding, mean_read_length=mean,
                     read_length_stddev=stddev,
                     mean_long_read_length=long_mean,
                     long_read_length_stddev=long_stddev,
                     coverage=coverage, long_read_coverage=long_coverage,
                     read_bases=tot_length)


def check_ref_needed(args, lib_type):
    if lib_type == "long" :
        if (args.genome_size is 0 and args.reference):
            raise ValueError("Please supply a genome size or reference " +
                             "sequence when running a long-read assembly")

def get_config_path():
    resource_package = pkg_resources.Requirement.parse("BugBuilder")
    config_path = '/'.join(('BugBuilder','config_data', 'BugBuilder.yaml'))
    config = pkg_resources.resource_filename(resource_package, config_path)
    return config


def parse_config(config_file):
    """ Read the config file, make a namespace object out of it
    """
    with open(config_file, 'r') as stream:
        try:
            yamldata = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            raise yaml.YAMLError("error parsing config file!")
    newns = argparse.Namespace()
    newns.__dict__ = yamldata
    return(newns)


def get_assemblers(args, config, paired):
    if len(args.assemblers) > 2:
        raise ValueError("A maximum of 2 assemblers can be requested")
    # sanity check assemblers where these have been manually requested...
    for assembler in args.assemblers:
        assembler_name = os.path.basename(assembler)
        if shutil.which(assembler) is None:
            raise ValueError("%s is nt in PATH! exiting" % assembler)
    # If no assembler is requested, we need to select one based on the
    # 'assembly_type' from the configuration file....
    if len(args.assemblers) > 1:
        for category in config.assembler_categories:
            if category.name == reads_ns.lib_type:
                args.assemblers.append(category.assemblers)
    # then check the requested assemblers are appropriate for
    # the sequence characteristics...
    for assembler in args.assemblers:
        assemblers = config.assemblers
        for conf_assembler in assemblers:
            if conf_assembler.name != assembler:
                continue
            if paired and conf_assembler.command_se:
                raise ValueError("%s requires paired reads, but you only " +
                                 "specified one fastq file..." % assembler)
            if conf_assembler.min_length and reads_ns.read_length < conf_assembler.min_length:
                raise ValueError("%s does not support reads less than %d" % \
                                 (assembler, conf_assembler.min_length))
            if conf_assembler.max_length and reads_ns.read_length > conf_assembler.max_length:
                raise ValueError("%s does not support reads greather than %d" % \
                                 (assembler, conf_assembler.max_length))
            if conf_assembler.max_length and reads_ns.read_length > conf_assembler.max_length:
                raise ValueError("%s requires the library insert size and " +
                                 "stddev to be provided. Please add the " +
                                 "--insert-size and --insert-stddev " +
                                 "parameters, or provide a reference " +
                                 "sequence" % assembler)


def get_scaffolder_and_linkage(args, config, paired):
    if args.scaffolder is None:
        return None
    if args.scaffolder not in [x.name for x in config.scaffolders]:
        raise ValueError("%s not an available scaffolder!" %args.scaffolder)
    linkage_evidence = None
    for conf_scaffolder in config.scaffolders:
        if conf_scaffolder.name == args.scaffolder:
            if conf_scaffolder.linkage_evidence == "paired-ends" and not paired:
                    raise ValueError("%s requires paired reads, but you " +
                                     "only specified one fastq file." % \
                                     args.scaffolder)
            elif "align" in conf_scaffolder.linkage_evidence and args.reference is None:
                raise ValueError("%s requires a reference for alignment, " +
                                 "but none is specified." % args.scaffolder)

            return conf_scaffolder.linkage_evidence


def get_merger_and_linkage(args, config, paired):
    if args.merge_method is None:
        return None
    if args.merge_method not in [x.name for x in config.merge_tools]:
        raise ValueError("%s not an available merge_method!" % args.merge_method)
    linkage_evidence = None
    for merger_methods in config.merge_methods:
        if merger_methods.name == args.merge_method:
            if merge_methods.allow_scaffolding:
                args.scaffolder = None

def get_finisher(args, config, paired):
    if args.finisher is None:
        return None
    if args.finisher not in [x.name for x in config.finishers]:
        raise ValueError("%s not an available finisher!" %args.finisher)
    for conf_finisher in config.finishers:
        if conf_finisher.name == args.finisher:
            if args.reference is None and conf_finisher.ref_required:
                raise ValueError("%s requires a reference." % \
                                     args.finisher)
            elif not paired and conf_finisher.paired_reads:
                raise ValueError("%s requires paired reads" % args.finisher)

def get_varcaller(args, config, paired):
    if args.varcaller is None:
        return None
    if args.varcaller not in [x.name for x in config.varcallers]:
        raise ValueError("%s not an available varcaller!" %args.varcaller)
    for conf_varcallers in config.varcallers:
        if conf_varcallers.name == args.varcaller:
            if args.reference is None and conf_varcallers.ref_required:
                raise ValueError("%s requires reference, but you " +
                                 "only specified one fastq file." % \
                                 args.varcaller)


def select_tools(args, paired, config):
    """
    select_tools
    options, while checking for validity of selections...

    required params: $ ($config hash ref)
                 $ (arrayref of requested assemblers)
                 $ (name of requested scaffolder)
                 $ (name of requested merge tool)
                 $ (name of requested finishing tool)
                 $ (name of requested variant caller)
                 $ (paired - flag to indicate we are using paired reads
                 $ (reference - indicates a reference has been provided)
                 $ (lib_type - category of sequence i.e. short_illumina)
                 $ (mean read length)

    returns          $ (arrayref of assemblers to use)
                 $ (name of scaffolder to use)

    """
    get_assemblers(args, config, paired)
    args.scaffolder, linkage_evidence = get_scaffolder_and_linkage(
        args=args, config=config, paired=paired)

    get_merger_and_linkage(args, config, paired)
    get_finisher(args, config, paired)
    get_varcaller(args, config, paired)

    print ("linkage: %s" % linkage_evidence)
    return Namespace(assemblers=args.assemblers,
                     scaffolder=args.scaffolder,
                     scaffold_type=linkage_evidence,
                     merge_method=args.merger,
                     finisher=args.finisher,
                     varcall=args.varcaller)


def get_downsampling(args, config):
    assembler_conf = config.assemblers
    for conf_assembler in config.assemblers:
        if conf_assembler.name in args.assemblers:
            return True
    return False


def make_fastqc_cmd(args, fastqc_dir):
    cmd = \
        "fastqc -t {4} --extract -o {0}{1}{2}{3} > {0}fastqc.log 2>&1".format(
            fastqc_dir,
            " " + args.fastq1 if args.fastq1 is not None else "",
            " " + args.fastq2 if args.fastq2 is not None else "",
            " " + args.long_fastq if args.long_fastq is not None else "",
            args.threads)
    return cmd

def run_fastqc(reads_ns, tmpdir, logger=None):
    """
    Carries out QC analysis using fastqc...

    required params: $ (tmpdir)
                 $ (type)

    returns: $ (0)
    """
    logger.info("Running FastQC...")
    fastqc_dir = os.path.join(tmpdir, "fastqc", "")
    os.makedirs(fastqc_dir)
    fastqc_cmd = make_fastqc_cmd(args, tmpdir)
    fastqc_res = subprocess.run(fastqc_cmd,
                                shell=sys.platform != "win32",
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                check=False)
    if result.returncode != 0:
        logger.warning("Error running fastqc with command %s", fastqc_cmd)
        sys.exit(1)
    report = []
    fails = 0
    for f in [args.fastq1, args.fastq2, args.long_fastq]:
        name = os.path.splitext(os.path.basename(fastqc))[0] + "_fastqc"
        report.append(name + "\n----------------------------------")
        with open(os.path.join(
                fastqc_dir, name + "_fastqc", "summary.txt"), "r") as s:
            for line in s:
                if "FAIL" in line:
                    fails = fails + 1
                report.append(s.split("\t")[:2])
        shutil.copyfile(os.path.join(fastqc_dir, name + "_fastqc.html"),
                        os.path.join(tmpdir, name + "_fastqc.html"))
    for line in report:
        logger.info(line)
    if fails > 0:
        if reads_ns.type in ["long", "hybrid", "de_fere"]:
            logger.info("NB: Reported quality issues from_fastq are normal " +
                        "when analysing PacBio sequence with FastQC...");
        else:
            logger.warning("Some QC tests indicate quality issues with this " +
                           "data.Please examine the fastqc outputs for " +
                           "these reads")

def make_sickle_cmd(args, tmpdir):
    if args.fastq2:
        cmd = str("sickle pe -f {0} -r {1} -t {2} -q {3} -l {4} -o " +
                  "{5}read1.fastq -p {5}read2.fastq -s {5}singles.fastq" +
                  "> {5}sickle.log").format(
                      args.fastq1,  # 0
                      args.fastq2,  # 1
                      reads_ns.encoding,  # 2
                      args.trim_qv,  #3
                      args.min_length,  #4
                      tmpdir)
    else:
        cmd = str("sickle se -f {0} -t {1} -q {2} -l {3} -o {4}read1.fastq " +
                  "> {4}sickle.log").format(
                      args.fastq1,  # 0
                      reads_ns.encoding,  # 1
                      args.trim_qv,  #2
                      args.min_length,  #3
                      tmpdir)
    return cmd

def quality_trim_reads(args, tmpdir, config, reads_ns):
    """
    Quality trims reads using sickle

    required params: $ (tmpdir)
                 $ (encoding)

    returns          $ (0)
    """
    if encoding is None:
        logger.warning("Sequence quality encoding could not be determined;" +
                       "Sequence quality trimming will be skipped...")
    trim_dir = os.path_join(tmpdir, "qc_trim", "")
    os.makedirs(trim_dir)
    trim_cmd = make_sickle_cmd(args, tmpdir)
    trim_res = subprocess.run(trim_cmd,
                                shell=sys.platform != "win32",
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                check=False)
    if trim_res.returncode != 0:
        logger.warning("Error running trimmer with command %s", trim_cmd)
        sys.exit(1)
    else:
        # reasign the reads with the trimmed treads
        args.fastq1 = os.path.join(trim_dir, "read1.fastq")
        args.fastq2 = os.path.join(trim_dir, "read2.fastq")
    kept = 0
    disc = 0
    logger.info("Quality Trimming result:s")
    with open(os.path.join(trim_dir, "sickle.log"), "r") as inf:
        for line in inf:
            logger.info(line)
            if "kept" in line:
                kept = kept + int(line.split(": ")[1].split("(")[0].strip())
            elif "discarded" in line:
                disc = desc + int(line.split(": ")[1].split("(")[0].strip())
            else:
                pass
    if float(disc / kept) * 100 > 10:
        logger.warning(">10\% of reads discarded during read trimming. "+
                       "Please examine the FastQC outputs...")


def make_seqtk_ds_cmd(args, reads_ns, new_coverage, outdir, logger):
    assert is_instance(new_coverage, int), "new_coverage must be an int")
    frac = float(new_coverage / reads_ns.coverage)
    cmd_list = []
    logger.info("downsampleing to %f X  to %f X (%f)",
                reads_ns.coverage, new_coverage , frac )
    for reads in [args.fastq1, args.fastq2]:
        if os.path.exists(reads):
            cmd_list.append("{0} sample -s 100 {1} {2} > {3}".format(
                config.seqtk,  reads, frac,
                os.path.join(outdir, os.path.basename(reads))))
    return cmd_list


def downsample_reads(args, reads_ns, config, new_cov=100):
    """
    Downsamples reads by default to a maximum of 100x if higher
    coverages are present, or to specified coverage

    required params: $ (tmpdir)
                 $ (config data)
                 $ (coverage required)
                 $ (original coverage)
                 $ (number of bases in original reds)
                 $ (mean read length)
    """
    ds_dir = os.path_join(tmpdir, "seqtk_dir", "")
    os.makedirs(ds_dir)
    ds_cmds = make_seqtk_ds_cmd(args, reads_ns, new_coverage=new_cov,
                                outdir=ds_dir, logger=logger,)
    for cmd in ds_cmds:
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    args.fastq1 = os.path.join(ds_dir, os.path.basename(args.fastq1))
    args.fastq2 = os.path.join(ds_dir, os.path.basename(args.fastq2))
    return new_cov


def make_bwa_cmds(args, config, outdir, ref, reads_ns, fastq1, fastq2):
    """ map with bwa
    returns a list of cmds and the path to the mapping file, as a tuple
    """
    cmd_list = []
    index_cmd = "{0} index {1}  > {2} bwa_index.log 2>&1".format(
        config.bwa, ref, outdir)
    cmd_list.append(index_cmd)
    # Use bwa-bwt for 'short' reads less than 100 bp, and bwa-mem for longer reads
    if reads_ns.read_length <= 100:
        cmdF = "{0} aln -t {1} {2} {3} > {4}read1.sai".format(
            config.bwa, args.threads, ref, fastq1, outdir)
        cmd_list.append(cmdF)
        if fastq2 is not None:
            cmdR = "{0} aln -t {1} {2} {3} > {4}read2.sai".format(
                config.bwa, args.threads, ref, fastq2, outdir)
            cmd_list.append(cmdR)
            sampe_cmd = str("{0} sampe {1} {2}read1.sai {2}read2.sai {3} {4}" +
                            " 2> {2}sampe.log > {2}mapping.sam").format(
                                config.bwa, ref, outdir, fastq1, fastq2)
        else:
            sampe_cmd = str("{0} samse {1} {2}read1.sai {3} " +
                            "2> {2}sampe.log > {2}mapping.sam").format(
                                config.bwa, ref, outdir, fastq1)
        cmd_list.append(sampe_cmd)

    else: # if long reads
        if fastq2 is None:
            mem_cmd = "{0} mem -t {1} -M {2} {3} >{4}mapping.sam 2>{4}bwa_mem.log".format(
                config.bwa, args.threads, ref, fastq1, outdir)
        else:
            mem_cmd = "{0} mem -t {1} -M {2} {3} {5} >{4}mapping.sam 2>{4}bwa_mem.log".format(
                config.bwa, args.threads, ref, fastq1, outdir, fastq2)
        cmd_list.append(mem_cmd)

    return (cmd_list, os.path.join(outdir, "mapping.sam"))


def make_samtools_cmds(config, mapping, outdir, sorted_bam):
    convert = str("{0} view -q 10 -Sb {1} 2>{2}samtoolsview.log | {0} sort - >" +
                  "{3}",).format(
                      config.samtools, mapping, outdir, sorted_bam)
    index = str("{0} index {1} 2>{2}samtools_index.log").format(
                      config.samtools, sorted_bam, outdir)
    return [convert, index]


def align_reads(dirname, downsample, logger):
    """
    Maps reads to reference using bwa

    required params: $ (tmp directory);
                 $ (reference);
                 $ (read length)
                 $ (flag to indicate downsampling required...)

               : $ (0)
    """
    align_dir = os.path_join(tmpdir, "align_" + dirname, "")
    bwa_dir = os.path_join(align_dir, "bwa", "")
    samtools_dir = os.path_join(tmpdir, "samtools" + dirname, "")
    seqtk_dir = os.path_join(tmpdir, "seqtk" + dirname, "")
    for d in [align_dir, bwa_dir, samtools_dir, seqtk_dir]:
        os.makedirs(d)
    bwa_reference = os.path.join(bwa_dir, os.path.basename(args.reference))
    shutil.copyfile(args.reference, bwa_reference)
    if downsample:
        logger.info("Downsampling reads for insert-size estimation...")
        ds_cmds = make_seqtk_ds_cmd(args, reads_ns, 10, outdir=seqtk_dir, logger=logger)
        for cmd in ds_cmds:
            subprocess.run(cmd,
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
        fastq1 = os.path.join(seqtk_dir, os.path.basename(args.fastq1))
        fastq2 = os.path.join(seqtk_dir, os.path.basename(args.fastq2))
    else:
        fastq1 = args.fastq1
        fastq2 = args.fastq2
    logger.info("BWA aligning reads vs $reference...")
    bwa_cmds, mapping = make_bwa_cmds(args, config, outdir, ref=bwa_reference,
                                      reads_ns=reads_ns, fastq1=fastq1,
                                      fastq2=fastq2)
    samtools_cmds =  make_samtools_cmds(config, mapping, outdir, sorted_bam)
    for cmd in bwa_cmds + samtools_cmds:
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    return sorted_bam


def main(args):
    if args is None:
        args = get_args()
    assembler_list = match_assembler_args(
        assembler=args.assemblers,
        assembler_args=args.assembler_args)
    check_file_paths(args)

    seq_ids = [x.id for x in  list(SeqIO.parse(args.reference, "fasta"))]

    organism = "{0}{1}{2}".format(args.genus, " " + args.species, " sp. " + args.strain)
    dt = arrow.utcnow()
    if not args.outdir:
        if "unknown" in organism:
            args.outdir = "BugBuilder_" + dt.local.format('YYYY-MM-DD_HHmmss')
        else:
            args.outdir = "BugBuilder_" + organism.replace(" ", "_")
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    logger = logging

    # my $log = $out_dir . ".log";
    log = "~/bugbuilder.log"
    PAIRED = False if not args.fastq2 else True

    logger.debug( "Logging output to %s", log)

    logger.info("Welcome to BugBuilder")
    logger.info("\nPreparing to build your bug...\n\n")

    tmpdir = setup_tempdir(args, output_root)

    reference = os.path.basename(args.reference)
    #  need to know a bit about the provided reads so we can act appropriately
    reads_ns  = assess_reads(args=args, conf=conf, platform=platform,
                                    tempdir=tmpdir, logger=logger)
    check_ref_needed(args=args, lib_type=reads_ns.lib_type)
    config = get_config_path()
    tools = select_tools(args, paired, config)

    ################################################################
    downsample_reads = get_downsampling(args, config)
    if not args.skip_fastqc and config.fastqc is not None:
       run_fastqc(reads_ns, tmpdir, logger=None)
    if reads_ns.mean_read_length is not None:
        if reads_ns.mean_read_length < 50 and reads_ns.mean_read_length < args.trim_length:
            logger.info("trim-length reset to 25 due to mean read length: %i",
                        reads_ns.mean_read_length)
            args.trim_length = 25
    if reads_ns.lib_type == "long" or args.skip_trim_reads:
        quality_trim_reads(args, tmpdir, config, reads_ns)
    if (reads_ns.coverage is not None and reads_ns.coverage > 100) and \
       (downsample_reads or args.downsample):
         downsampled_coverage = downsample_reads(args, reads_ns, config, new_cov=100)

    if args.fastq2 and args.reference:
        sorted_bam = align_reads(dirname, downsample=True, logger=logger)
# 	( $insert_size, $stddev ) = get_insert_stats( "$tmpdir", $reference );
#     }
#     my $paired_str;
#     ($paired) ? ( $paired_str = "Paired" ) : ( $paired_str = "Fragment" );
#     my $tb = Text::ASCIITable->new( { 'headingText' => 'Read Information' } );
#     $tb->setCols( 'Aaaaaaaaargh', 'Phwwwwwwweb' );
#     $tb->setOptions( { 'hide_HeadRow' => 1 } );
#     $tb->alignCol( 'Aaaaaaaaargh', 'left' );
#     $tb->alignCol( 'Phwwwwwwweb',  'left' );
#     $tb->addRow( "Mean Read Length",                  $mean_read_length )          if ($mean_read_length);
#     $tb->addRow( "Read Length Standard Deviation",    $read_length_stddev )        if ($read_length_stddev);
#     $tb->addRow( "Insert size",                       $insert_size )               if ($insert_size);
#     $tb->addRow( "Insert size Standard Deviation",    $stddev )                    if ($stddev);
#     $tb->addRow( "Mean Long Read Length",             $mean_long_read_length )     if ($mean_long_read_length);
#     $tb->addRow( "Mean Long Read Standard Deviation", $long_read_length_stddev )   if ($long_read_length_stddev);
#     $tb->addRow( "Library type",                      $paired_str );
#     $tb->addRow( "Platform",                          $platform )                  if ($platform);
#     $tb->addRow( "Quality Encoding",                  $encoding );
#     $tb->addRow( "Projected Coverage",                "${coverage}x" )             if ($coverage);
#     $tb->addRow( "Projected Coverage (Downsampled)",  "${downsampled_coverage}x" ) if ($downsampled_coverage);
#     $tb->addRow( "Projected Long Read Coverage",      "${long_read_coverage}x" )   if ($long_read_coverage);
#     print "$tb\n";

#     $tb = Text::ASCIITable->new( { 'headingText' => 'Assembly Information' } );
#     $tb->setCols( 'Aaaaaaaaargh', 'Phwwwwwwweb' );
#     $tb->setOptions( { 'hide_HeadRow' => 1 } );
#     $tb->alignCol( 'Aaaaaaaaargh', 'left' );
#     $tb->alignCol( 'Phwwwwwwweb',  'left' );
#     $tb->addRow( "Assembly category",        $type );
#     $tb->addRow( "Selected assemblers",      join( ",", @assemblers ) );
#     $tb->addRow( "Selected assembly merger", $merge_method ) if ($merge_method);
#     $tb->addRow( "Selected scaffolder",      $scaffolder ) if ($scaffolder);
#     $tb->addRow( "Selected finisher",        $finisher ) if ($finisher);
#     $tb->addRow( "Selected variant caller",  $varcall ) if ($varcall);
#     $tb->addRow( "Trim QV",                  $trim_qv ) if ($trim_reads);
#     $tb->addRow( "Trim length",              $trim_length ) if ($trim_reads);
#     $tb->addRow( "Split Origin",             $split_origin );
#     print "$tb\n";

#     my $pm = new Parallel::ForkManager(2);
#     my $id_ok = 1;
#     if ( !$already_assembled ){
# 	print "Running in $mode mode\n";
# 	my $ret;
# 	my $link = 0;

# 	#if we are running a single assembler we need to symlink its contig output
# 	# into $tmpdir - pass a $link flag to run_assembler
# 	$link = 1 if ( $#assemblers == 0 );
# 	for ( my $i = 0 ; $i <= $#assemblers ; $i++ ) {

# 	    my $assembler = $assemblers[$i];
# 	    my $args      = $assembler_args[$i];

# 	    run_assembler(
# 		$tmpdir,      $assembler, $reference,   $args,     $link,
# 		$type,        $encoding,  $genome_size, $platform, $mean_read_length,
# 		$insert_size, $stddev,    $threads
# 		);
# 	}
# 	merge_assemblies( $tmpdir, \@assemblers, $merge_method ) if ( $#assemblers == 1 );
# 	if ( $scaffolder && $reference ) {
# 	    $id_ok = check_id($tmpdir);
# 	}
#     } else { # end if already assembled
#     	print "using assemblies from another directory\n";
#     	$id_ok = 1;  # it better be, that is
# 	$tmpdir = $scratchdir;
# 	my $oridir = $tmpdir . "/origin";
# 	rmtree "$tmpdir/origin" or print "Unable to remove origin dir: $!\n";
# 	rmtree "$tmpdir/agp" or print "Unable to remove agp dir: $!\n";
# 	unlink "$tmpdir/scaffolds.agp" or print "unable to delete scaffolds.agp\n";
# 	unlink "$tmpdir/contigs.embl" or print "unable to delete contigs.embl\n";
# 	unlink "$tmpdir/scaffolds.embl" or print "unable to delete scaffolds.embl\n";
# 	# unlink "$tmpdir/scaffolds.fasta" or print "unable to delete scaffolds.fasta\n";
# 	# unlink "$tmpdir/contigs.fasta" or print "unable to delete scaffolds.contigs\n";
# 	rmtree "$tmpdir/annotation_merge" or print "unable to delete merge\n";
# 	rmtree "$tmpdir/cgview" or print "unable to delete cgview\n";
# 	rmtree "$tmpdir/comparisons" or print "unable to delete comparisons\n";
# 	unlink "$tmpdir/contigs_cgview.png" or print "unable to delete contigs_cgview.png\n";
# 	unlink "$tmpdir/scaffolds_cgview.png" or print "unable to delete scaffolds_cgview.png\n";
# 	unlink glob "$tmpdir/comparison_vs*" or print "unable to delete comparisons";
#     	print "idok: $id_ok\n";
#     }
#     run_scaffolder( $tmpdir, "$tmpdir/reference_parsed_ids.fasta",
#                     $scaffolder, $scaffolder_args, $insert_size, $stddev, 1, "$tmpdir/contigs.fasta", $mean_read_length,
#                     $threads )
#       if ( $scaffolder && ( $id_ok && $scaffold_type eq 'align_genus' )
#            || ( $id_ok && $scaffold_type eq 'paired-ends' )
#            # why is that id_ok negated?
#            || ( !$id_ok && $scaffold_type eq 'paired-ends' ) );

#     if ( -e "$tmpdir/scaffolds.fasta" && $reference ) {

#         #we may have been given a reference for a long read assembly but no scaffolder is used for these by default
#         $scaffolder = "mauve" if ( !$scaffolder );

#         find_origin( $tmpdir, $scaffolder, $scaffolder_args, "$tmpdir/reference_parsed_ids.fasta",
#                      $insert_size, $stddev, $mean_read_length, $threads )
#           if ($split_origin);
#         order_scaffolds( $tmpdir, basename($reference) ) if ($id_ok);

#         finish_assembly( $tmpdir, $finisher, $insert_size, $stddev, $encoding, $threads ) if ( $finisher && $id_ok );

#     }

#     my $gaps;

#     if ( -e "$tmpdir/scaffolds.fasta" ) {
#         $gaps = build_agp( $tmpdir, $organism, $mode, $scaffold_type );
#     }

#     # sequence stable from this point, only annotations altered
#     for my $i (qw(1 2)) {
# 	print "pm: $i \n";
#         $pm->start and next();
#         if ( $i == 1 ) {

#             ##amosvalidate fails if we don't have mate-pairs
#             if ( -e "$tmpdir/read2.fastq" && ( $mode eq 'draft' ) ) {
#                 my $seq_file;
#                 ( -e "$tmpdir/scaffolds.fasta" ) ? ( $seq_file = 'scaffolds.fasta' ) : ( $seq_file = 'contigs.fasta' );

#                 align_reads( $tmpdir, $seq_file, $mean_read_length );
#                 amosvalidate( $tmpdir, $insert_size, $stddev );
#             }
#         }
#         elsif ( $i == 2 ) {

#             run_prokka( $tmpdir, $genus, $species, $strain, $locustag, $centre );
#         }

#         $pm->finish();
#     }

#     $pm->wait_all_children();

#     my $amosvalidate_results;
#     if ( -e "$tmpdir/read2.fastq" && ( $mode eq 'draft' ) ) {
#         get_contig_to_iid_mapping($tmpdir);
#         $amosvalidate_results = summarise_amosvalidate($tmpdir);
#     }

#     run_varcaller( $tmpdir, $varcall, $threads, $mean_read_length ) if ($varcall);
#     merge_annotations( $tmpdir, $amosvalidate_results, $gaps, $genus, $species, $strain );
#     #  This kept throwing an error about Bio::SeqIO
#     # run_cgview($tmpdir);

#     build_comparisons( $tmpdir, basename($reference), $organism ) if ($reference);

#     message("Final Assembly Statistics...");
#     get_contig_stats( "$tmpdir/contigs.fasta", 'contigs' );
#     get_contig_stats( "$tmpdir/scaffolds.fasta", 'scaffolds' ) if ( -e "$tmpdir/scaffolds.fasta" );
#     my $emblfile;
#     ( -e "$tmpdir/scaffolds.embl" ) ? ( $emblfile = "$tmpdir/scaffolds.embl" ) : ( $emblfile = "$tmpdir/contigs.embl" );

#     my $io = Bio::SeqIO->new( -format => 'embl', -file => "$emblfile" );
#     my ( $cds, $tRNA, $rRNA );
#     while ( my $contig = $io->next_seq() ) {
#         foreach my $feature ( $contig->get_SeqFeatures() ) {
#             $cds++  if ( $feature->primary_tag eq 'CDS' );
#             $tRNA++ if ( $feature->primary_tag eq 'tRNA' );
#             $rRNA++ if ( $feature->primary_tag eq 'rRNA' );
#         }
#     }

#     $tb = Text::ASCIITable->new();
#     $tb->setCols( "Feature Type", "Number" );
#     $tb->addRow( "CDS",  $cds );
#     $tb->addRow( "tRNA", $tRNA );
#     $tb->addRow( "rRNA", $rRNA );
#     print "\nAnnotated features\n==================\n\n";
#     print $tb, "\n";

#     return_results( $tmpdir, $out_dir, $prefix, $keepall, $mode );
#     chdir $orig_dir or warn "Failed to chdir:$ !";

#     message("All done...");

# }

# ######################################################################
# #
# # set_paths
# #
# # Setups the PATH, PYTHONPATH and PERL5LIB environmental variables
# # according to the specifications in the config...
# #
# # required arguments: $ (config hashref)
# #
# # returns:            0
# #
# ######################################################################

# sub set_paths {

#     my $config = shift;

#     # Although we execute tools using fully qualified paths,
#     # some of these expect certain things to be on path...
#     $ENV{'PATH'} = $ENV{'PATH'} . ':' . $FindBin::Bin;
#     $ENV{'PATH'} = $ENV{'PATH'} . ":" . "$FindBin::Bin/../packages/bin";
#     $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'R_dir'} . '/bin'
#       if ( $config->{'R_dir'} );
#     $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'blast_dir'}
#       if ( $config->{'blast_dir'} );
#     $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'mummer_dir'}
#       if ( $config->{'mummer_dir'} );
#     $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'ncbi_utils_dir'}
#       if ( $config->{'ncbi_utils_dir'} );
#     $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'tbl2asn_dir'}
#       if ( $config->{'tbl2asn_dir'} );
#     $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'asn2gb_dir'}
#       if ( $config->{'asn2gb_dir'} );
#     $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'aragorn_dir'} . '/bin'
#       if ( $config->{'aragorn_dir'} );
#     $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'prodigal_dir'}
#       if ( $config->{'prodigal_dir'} );
#     $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'hmmer3_dir'} . '/bin'
#       if ( $config->{'hmmer3_dir'} );

#     #$ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'rnammer_dir'}
#     #  if ( $config->{'rnammer_dir'} );
#     $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'infernal_dir'} . '/bin'
#       if ( $config->{'infernal_dir'} );
#     $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'prokka_dir'}
#       if ( $config->{'prokka_dir'} );
#     $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'barrnap_dir'}
#       if ( $config->{'barrnap_dir'} );
#     $ENV{'PATH'} = $ENV{'PATH'} . ":" . $config->{'mauve_dir'} . '/linux-x64'
#       if ( $config->{'mauve_dir'} );
#     $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'abyss_sealer_dir'}
#       if ( $config->{'abyss_sealer_dir'} );

#     if ( $config->{'python_lib_path'} ) {
#         if ( defined( $ENV{'PYTHONPATH'} ) ) {
#             $ENV{'PYTHONPATH'} = "$ENV{'PYTHONPATH'}:" . $config->{'python_lib_path'};
#         }
#         else {
#             $ENV{'PYTHONPATH'} = $config->{'python_lib_path'};
#         }
#     }

#     if ( $config->{'perl_lib_path'} ) {
#         if ( $ENV{'PERL5LIB'} ) {
#             $ENV{'PERL5LIB'} = "$ENV{'PERL5LIB'}:" . $config->{'perl_lib_path'};
#         }
#         else {
#             $ENV{'PERL5LIB'} = $config->{'perl_lib_path'};
#         }
#     }
#     return;
# }


# #######################################################################
# #
# ######################################################################
# ######################################################################
# #
# # return_results
# #
# # Copies results back from tmpdir, returning full working directory
# # if dircopy argument specified
# #
# # required params: $ (tmpdir)
# #                  $ (strain)
# #                  $ (dircopy - flag to indicate entire directory
# #                     should be returned)
# #                  $ (mode - draft mode also needs amos bank copying)
# #
# # returns        : $ (0)
# #
# ######################################################################

# sub return_results {

#     my $tmpdir  = shift;
#     my $dir     = shift;
#     my $prefix  = shift;
#     my $dircopy = shift;
#     my $mode    = shift;

#     if ($dircopy) {
#         dircopy( $tmpdir, "$dir" )
#           or die "Error copying $tmpdir: $!";
#     }
#     else {
#         my @files = qw(annotated.embl contigs.fasta scaffolds.fasta scaffolds.embl scaffolds.agp
#           unplaced_contigs.fasta BugBuilder.log read1_fastqc.html read2_fastqc.html
#           scaffolds_cgview.png contigs_cgview.png circleator.png circleator.svg reference.variants.vcf
#           );

#         opendir TMP, "$tmpdir" or die "Error opening $tmpdir: $!";
#         my @all_files = readdir TMP;
#         close TMP;

#         foreach my $pattern (qw(blastout png)) {
#             my @found = grep /$pattern/, @all_files;
#             push @files, @found;
#         }

#         mkdir "$dir"
#           or die "Error creating $dir: $!";

#         foreach my $file (@files) {
#             my $outfile;
#             if ($prefix) {
#                 $outfile = "$dir/${prefix}_${file}";
#             }
#             else {
#                 $outfile = "$dir/$file";
#             }
#             my $target;
#             if ( -l "$tmpdir/$file" ) {
#                 $target = readlink("$tmpdir/$file");
#             }
#             else {
#                 $target = "$tmpdir/$file";
#             }
#             copy( "$target", "$outfile" )
#               or die "Error copying $file: $!"
#               if ( -e "$tmpdir/$file" );
#         }
#         if ( $mode eq 'draft' ) {
#             my $outfile;
#             if ($prefix) {
#                 $outfile = "$dir/${prefix}_assembly.bnk";
#             }
#             else {
#                 $outfile = "$dir/assembly.bnk";
#             }
#             dircopy( "$tmpdir/amos/assembly.bnk", "$outfile" )
#               or die "Error copying $tmpdir/amos/assembly.bnk: $!";
#         }

#     }

# }










# #######################################################################
# #
# # run_assembler
# #
# # Runs specified assembler on fastq files
# #
# # required params: $ (tmpdir)
# #                  $ (assembler name)
# #                  $ (reference fasta sequece - probably not needed)
# #                  $ (arguments to pass to assembler)
# #                  $ (link - flag to indicate contigs should be symlinked into tmpdir)
# #                  $ (category - type of assembly to run, since an assembler may fall into more
# #                  than one category)
# #                  $ (encoding - some assemblers need explicitly telling)
# #                  $ (genome size)
# #                  $ (average read length)
# #	           $ (insert size)
# #                  $ (insert size stddev)
# #
# # returns          $ (0)
# #
# #######################################################################

# sub run_assembler {

#     my $tmpdir      = shift;
#     my $assembler   = shift;
#     my $reference   = shift;
#     my $args        = shift;
#     my $link        = shift;
#     my $category    = shift;
#     my $encoding    = shift;
#     my $genome_size = shift;
#     my $platform    = shift;
#     my $read_length = shift;
#     my $insert_size = shift;
#     my $stddev      = shift;
#     my $threads     = shift;

#     my ( $cmd, $contig_output, $scaffold_output, $create );

#     my $assemblers = $config->{'assemblers'};
#     foreach my $conf_assembler (@$assemblers) {
#         if ( $conf_assembler->{'name'} eq $assembler ) {
#             if ( $category eq 'hybrid' ) {
#                 $cmd = $conf_assembler->{'command_hybrid'};
#             }
#             elsif ( $category eq 'de_fere' ) {
#                 $cmd = $conf_assembler->{'command_de_fere'};
#             }
#             elsif ( -e "$tmpdir/read2.fastq" ) {
#                 $cmd = $conf_assembler->{'command_pe'};
#             }
#             else {
#                 $cmd = $conf_assembler->{'command_se'};
#             }
#             if ($args) {
#                 $cmd .= " $args";
#             }
#             elsif ( $conf_assembler->{'default_args'} ) {
#                 $cmd .= " " . $conf_assembler->{'default_args'};
#             }
#             $contig_output   = $conf_assembler->{'contig_output'};
#             $scaffold_output = $conf_assembler->{'scaffold_output'}
#               if ( $conf_assembler->{'scaffold_output'} );
#             $create = $conf_assembler->{'create_dir'};
#         }
#     }
#     die "Assembler $assembler is not defined" unless ($cmd);

#     if ($create) {
#         mkdir "$tmpdir/$assembler"
#           or die "could not create $tmpdir/$assembler: $! ";
#         chdir "$tmpdir/$assembler"
#           or die "could not chdir to $tmpdir/$assembler: $! ";
#     }

#     if ( $reference && !$genome_size ) {
#         my $io = Bio::SeqIO->new( -file => "$tmpdir/$reference", -format => 'fasta' );
#         while ( my $ref = $io->next_seq() ) {
#             $genome_size += $ref->length();
#         }
#     }
#     my $asm_dir = $config->{"${assembler}_dir"};
#     foreach ( $cmd, $contig_output, $scaffold_output ) {
#         s/__BUGBUILDER_BIN__/$FindBin::Bin/g;
#         s/__ASMDIR__/$asm_dir/;
#         s/__TMPDIR__/$tmpdir/g;
#         s/__FASTQ1__/$tmpdir\/read1.fastq/;
#         s/__FASTQ2__/$tmpdir\/read2.fastq/;
#         s/__ORIG_FASTQ1__/$tmpdir\/orig_read1.fastq/;
#         s/__ORIG_FASTQ2__/$tmpdir\/orig_read2.fastq/;
#         s/__DE_FERE_CONTIGS__/$tmpdir\/de_fere_contigs.fasta/;
#         s/__LONGFASTQ__/$tmpdir\/long.fastq/;
#         s/__REFERENCE__/$reference/ if ($reference); #
#         s/__CATEGORY__/$category/;
#         s/__ENCODING__/$encoding/;
#         s/__GENOME_SIZE__/$genome_size/ if defined($genome_size);
#         s/__PLATFORM__/$platform/       if ($platform);
#         s/__READ_LENGTH__/$read_length/ if ($read_length);
#         s/__INSSIZE__/$insert_size/     if ($insert_size);
#         s/__INSSD__/$stddev/            if ($stddev);
#         s/__THREADS__/$threads/;
#     }

#     $cmd .= " > $tmpdir/$assembler.log 2>&1 ";

#     message(" Starting $assembler assembly ... ");
#     system($cmd) == 0 or die "Error running $cmd: $! ";

#     my $header = "$assembler assembly statistics";
#     print $header . "\n" . '=' x length($header) . "\n\n";
#     get_contig_stats( "$contig_output", 'contigs' );

#     # Only generate scaffold stats for paired read alignments...
#     get_contig_stats( "$scaffold_output", 'scaffolds' )
#       if ( $scaffold_output && -e $scaffold_output && -e "$tmpdir/read2.fastq" );

#     # rename contigs/scaffolds for consistent naming, since we need to retrieve by
#     # id later, so it helps if we know what the ids look like...
#     if ( !$create ) {
#         my $outdir = dirname($contig_output);
#         chdir $outdir or die " Error changing to dir $outdir: $! ";
#     }

#     my $inIO = Bio::SeqIO->new( -format => 'fasta', -file => $contig_output );
#     my $outIO =
#       Bio::SeqIO->new( -format => 'fasta',
#                        -file   => ">BugBuilder.contigs.fasta" );
#     my $contig_count = 0;
#     while ( my $seq = $inIO->next_seq() ) {
#         my $contig_id = sprintf( "contig_%06s", ++$contig_count );
#         $seq->id($contig_id);
#         $outIO->write_seq($seq);
#     }

#     if ( $scaffold_output && ( -e $scaffold_output ) ) {    #&& ( -e "$tmpdir/read2.fastq" ) ) {
#         my $inIO = Bio::SeqIO->new( -format => 'fasta', -file => $scaffold_output );
#         my $outIO =
#           Bio::SeqIO->new( -format => 'fasta',
#                            -file   => ">BugBuilder.scaffolds.fasta" );
#         my $scaffold_count = 0;
#         while ( my $seq = $inIO->next_seq() ) {
#             my $scaffold_id = sprintf( "scaffold_%06s", ++$scaffold_count );
#             $seq->id($scaffold_id);
#             $outIO->write_seq($seq);
#         }
#     }

#     if ($link) {

#         # need to use File::Find::Rule to get the right path since some outputs are more nested
#         # than others...
#         my $contigs = ( File::Find::Rule->file()->name("BugBuilder.contigs.fasta")->in("$tmpdir/$assembler") )[0];
#         symlink( $contigs, "../contigs.fasta" )
#           or die "Error creating symlink: $!";
#         if ( $scaffold_output && ( -e $scaffold_output ) ) {
#             my $scaffolds =
#               ( File::Find::Rule->file()->name("BugBuilder.scaffolds.fasta")->in("$tmpdir/$assembler") )[0];
#             symlink( $scaffolds, "../scaffolds.fasta" )
#               or die "Error creating symlink: $!"
#               if ($scaffolds);
#         }
#     }

#     chdir $tmpdir or die "Could not chdir to $tmpdir: $!" if ($create);

#     return ();
# }

# ######################################################################
# #
# # merge_assemblies
# #
# # combines two assemblies using the selected merge-method
# #
# # required params: $ (tmpdir)
# #                  $ (arrayref of assemblers used)
# #                  $ (method)
# #                  $ (reference)
# #
# # returns        : $ (0)
# #
# ######################################################################

# sub merge_assemblies {

#     my $tmpdir     = shift;
#     my $assemblers = shift;
#     my $method     = shift;
#     my $reference  = shift;

#     message("Merging assemblies ($method)...");

#     my ( $cmd, $create_dir, $contig_output );
#     my $merge_tools = $config->{'merge_tools'};
#     foreach my $tool (@$merge_tools) {
#         my $name = $tool->{'name'};
#         if ( $name eq $method ) {
#             $cmd           = $tool->{'command'};
#             $create_dir    = $tool->{'create_dir'};
#             $contig_output = $tool->{'contig_output'};
#         }
#         if ($create_dir) {
#             mkdir "$tmpdir/$method" or die "Error creating $tmpdir/$method: $! ";
#             chdir "$tmpdir/$method" or die "Error chdiring to $tmpdir/$method: $!";
#         }

#         $cmd =~ s/__BUGBUILDER_BIN__/$FindBin::Bin/;
#         $cmd =~ s/__TMPDIR__/$tmpdir/;
#         $cmd =~ s/__ASSEMB1__/$assemblers->[0]/;
#         $cmd =~ s/__ASSEMB2__/$assemblers->[1]/;
#         $cmd =~ s/__REFERENCE__/${tmpdir}\/reference.fasta/;

#         system($cmd) == 0 or die "Error running $cmd: $!";
#         chdir $tmpdir     or die " Error chdiring to $tmpdir: $! ";
#         symlink( "$method/$contig_output", "contigs.fasta" )
#           or die " Error creating symlink : $! ";

#         print "Merged assembly statistics : \n============================\n";
#         get_contig_stats( "$tmpdir/$method/$contig_output", 'contigs' );

#         return (0);

#     }
# }

# ######################################################################
# #
# # finish_assembly
# #
# # Carried out assembly finishing using selected method
# #
# # required params: $ (tmpdir)
# #                  $ (finisher)
# #                  $ (insert size)
# #                  $ (insert stddev)
# #                  $ (base quality encoding)
# #                  $ (no. threads)
# #
# # returns        : $ (0)
# #
# ######################################################################

# sub finish_assembly {

#     my $tmpdir        = shift;
#     my $finisher      = shift;
#     my $insert_size   = shift;
#     my $insert_stddev = shift;
#     my $encoding      = shift;
#     my $threads       = shift;

#     message("Finishing assembly ($finisher)...");

#     my ( $cmd, $create_dir );
#     my $finishers = $config->{'finishers'};
#     foreach my $tool (@$finishers) {
#         my $name = $tool->{'name'};
#         if ( $name eq $finisher ) {
#             $cmd        = $tool->{'command'};
#             $create_dir = $tool->{'create_dir'};
#         }
#     }
#     if ($create_dir) {
#         mkdir "$tmpdir/$finisher" or die "Error creating $tmpdir/$finisher: $! ";
#         chdir "$tmpdir/$finisher" or die "Error chdiring to $tmpdir/$finisher: $!";
#     }

#     $cmd =~ s/__BUGBUILDER_BIN__/$FindBin::Bin/;
#     $cmd =~ s/__TMPDIR__/$tmpdir/;
#     $cmd =~ s/__REFERENCE__/${tmpdir}\/reference.fasta/;
#     $cmd =~ s/__INSSIZE__/$insert_size/;
#     $cmd =~ s/__INSSD__/$insert_stddev/;
#     $cmd =~ s/__ENCODING__/$encoding/;
#     $cmd =~ s/__THREADS__/$threads/;

#     system("perl $cmd") == 0 or die "Error running $cmd: $!";
#     chdir $tmpdir     or die " Error chdiring to $tmpdir: $! ";

#     print "Finished assembly statistics : \n============================\n";
#     get_contig_stats( "$tmpdir/scaffolds.fasta", 'scaffolds' );

#     return (0);
# }




# ######################################################################
# #
# # get_insert_stats
# #
# # converts contigs.sam -> bam, sorts, indexes and generates
# # insert stats with Picard
# #
# # required params: $ (tmp directory);
# #                  $ (reference);
# #
# # returns        : $ (insert size)
# #                : $ (stddev)
# #
# ######################################################################

# sub get_insert_stats {

#     my $tmpdir    = shift;
#     my $reference = shift;

#     message(" Collecting insert size statistics ");

#     mkdir "$tmpdir/insert_stats"
#       or die "Error creating $tmpdir/insert_stats: $! "
#       if ( !-d "$tmpdir/insert_stats" );
#     chdir "$tmpdir/insert_stats"
#       or die "Error chdiring to $tmpdir/insert_stats : $! ";

#     # picard command for earlier versions....
#     #my $cmd = "java -jar " . $config->{'picard_dir'} . "CollectInsertSizeMetrics.jar ";
#     my $cmd = "java -jar " . $config->{'picard_dir'} . "picard.jar CollectInsertSizeMetrics ";
#     $cmd .= "INPUT=$tmpdir/bwa/$reference.bam HISTOGRAM_FILE=insert_histogram.pdf OUTPUT=insert_stats.txt ";
#     $cmd .=
#       " QUIET=true VERBOSITY=ERROR ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT > CollectInsertMetrics.log 2>&1";
#     system($cmd) == 0 or die " Error executing $cmd: $! ";

#     open STATS, "insert_stats.txt"
#       or die "Error opening insert_stats.txt: $! ";
#     my ( $min_insert, $max_insert, $insert, $stddev );
#   LINE: while ( my $line = <STATS> ) {
#         if ( $line =~ /^MEDIAN/ ) {
#             my @stats = split( /\t/, <STATS> );
#             $min_insert = $stats[2];
#             $max_insert = $stats[3];
#             $insert     = sprintf( "%d", $stats[4] );
#             $stddev     = sprintf( "%d", $stats[5] );
#             last LINE;
#         }
#     }
#     close STATS;

#     chdir $tmpdir or die "Error chdiring to $tmpdir: $!";

#     return ( $insert, $stddev );

# }

# ######################################################################
# #
# # check_id
# #
# # Checks wether the assembled contigs have sufficient identity to the
# # reference for reference-based scaffolding
# #
# # required params: $ (tmp directory);
# #
# # returns        : $ id_ok (boolean)
# #
# ######################################################################

# sub check_id {

#     my $tmpdir = shift;

#     message("Checking identity of assembly with reference...");

#     my $id_ok = 0;

#     mkdir "$tmpdir/id_check"
#       or die "Error creating $tmpdir/id_check: $! "
#       if ( !-d "$tmpdir/id_check" );
#     chdir "$tmpdir/id_check"
#       or die "Error chdiring to $tmpdir/id_check: $! ";

#     my $cmd = $config->{'blast_dir'}
#       . "/blastn -query $tmpdir/contigs.fasta -subject $tmpdir/reference_parsed_ids.fasta -outfmt 5 -evalue 0.01 -out blastout.xml 2>&1 > blastn.log";
#     system($cmd) == 0 or die "Error running $cmd: $!";

#     my $blio = Bio::SearchIO->new( -format => 'blastxml',
#                                    -file   => 'blastout.xml' );

#     my ( $aligned, $unaligned );
#     while ( my $result = $blio->next_result() ) {
#         foreach my $hit ( $result->hits() ) {
#             my $tiling = Bio::Search::Tiling::MapTiling->new($hit);
#             $aligned   += $tiling->num_aligned();
#             $unaligned += $tiling->num_unaligned();
#         }
#     }
#     my $percent_id = sprintf( '%.2f', ( $aligned / ( $aligned + $unaligned ) * 100 ) );
#     print "ID=$percent_id %\n";

#     if ( $percent_id > 80 ) {
#         $id_ok = 1;
#     }
#     else {
#         $id_ok = 0;
#         print RED, "WARNING", RESET, ": Percentage identitiy of assembly with reference looks to low (${percent_id}%\n";
#         print "Reference will not be used for scaffolding or ordering contigs\n";
#     }

#     chdir $tmpdir or die "Error chdiring to $tmpdir: $!";

#     return ($id_ok);

# }

# ######################################################################
# #
# # run_scaffolder
# #
# # Runs specified scaffolder....
# #
# # Scaffolders which don't separate unscaffolded contigs can be wrapped
# # in a script which created a "$scaffolder.contig_ids" file listing the
# # IDs of contigs scaffolded. These will then be used following scaffolding
# # to create our own file of unscaffolded contigs
# #
# # required params: $ (tmpdir)
# #                  $ (reference)
# #                  $ (scaffolder)
# #                  $ (scaffolder args)
# #                  $ (library insert size)
# #	           $ (library insert sd)
# #                  $ (run_id - appended to tmpdir to allow multiple runs)
# #                  $ (path to contigs to scaffold)
# #                  $ (mean read length)
# #
# # returns        : $ (linkage evidence type)
# #
# ######################################################################

# sub run_scaffolder {

#     my $tmpdir           = shift;
#     my $reference        = shift;
#     my $scaffolder       = shift;
#     my $scaffolder_args  = shift;
#     my $insert_size      = shift;
#     my $insert_sd        = shift;
#     my $run_id           = shift;
#     my $contigs          = shift;
#     my $mean_read_length = shift;
#     my $threads          = shift;

#     message(" Starting $scaffolder");

#     my ( $cmd, $scaffold_output, $unscaffolded_output, $create, $linkage_evidence, $default_args );
#     my $blast_dir = $config->{'blast_dir'};

#     my $scaffolders = $config->{'scaffolders'};
#     foreach my $conf_scaffolder (@$scaffolders) {
#         if ( lc( $conf_scaffolder->{'name'} ) eq lc($scaffolder) ) {
#             $cmd                 = $conf_scaffolder->{'command'};
#             $scaffold_output     = $conf_scaffolder->{'scaffold_output'};
#             $unscaffolded_output = $conf_scaffolder->{'unscaffolded_output'}
#               if ( $conf_scaffolder->{'unscaffolded_output'} );
#             $create           = $conf_scaffolder->{'create_dir'};
#             $linkage_evidence = $conf_scaffolder->{'linkage_evidence'};
#             $default_args     = $conf_scaffolder->{'default_args'};
#         }
#     }
#     die "Scaffolder $scaffolder is not defined" unless ($cmd);
#     my $run_dir .= "${tmpdir}/${scaffolder}_${run_id}";
#     mkdir("$run_dir") or die "Error creating $run_dir: $! ";
#     chdir("$run_dir") or die "Error in chdir $run_dir: $! ";

#     # Treat reference-based scaffolder separately from paired-read scaffolders,
#     # since we need to scaffold per-reference, which doesn't work if your not using one...
#     if ( $linkage_evidence eq 'align_genus' ) {

#         # If the reference contains multiple contigs, we first neeed to align out contigs to these
#         # to identify which contigs to scaffold against which reference, since
#         # some scaffolders targeted at bacteria don't handle multiple reference
#         # sequences
#         my $io = Bio::SeqIO->new( -format => 'fasta', -file => "$reference" );
#         my @ref_ids;
#         while ( my $seq = $io->next_seq() ) {
#             my $ref_id = $seq->display_id();

#             #$ref_id =~ s/lcl\|//;
#             $ref_id = parse_ref_id($ref_id);

#             #push @ref_ids, $seq->id();
#             push @ref_ids, $ref_id;
#             my $outIO = Bio::SeqIO->new( -format => 'fasta', -file => ">$run_dir/reference_${ref_id}" );
#             $outIO->write_seq($seq);
#         }

#         if ( $#ref_ids > 0 ) {

#             # blast indexing doesn't produce a workable index from a symlink, so need to copy the reference sequences
#             copy( $reference, "$run_dir/reference.fasta" )
#               or die "Error copying reference.fasta:$! ";
#             symlink( "$contigs", "$run_dir/contigs.fasta" );
#             my $fasta_db = Bio::DB::Fasta->new("$run_dir/contigs.fasta");

#             my $blast_cmd =
#               "$blast_dir/makeblastdb -in reference.fasta -dbtype nucl -parse_seqids 2>&1 > makeblastdb.log";
#             system($blast_cmd) == 0 or die "Error building reference.fasta blast database: $!";
#             $blast_cmd = "$blast_dir/blastn -query $contigs -task blastn -db reference.fasta "
#               . "-out clusters.blast  2>&1 > blastn.log";
#             system($blast_cmd) == 0 or die "Error executing $blast_cmd: $!";

#             # create hash of contigs per reference sequence, or unaligned
#             # only need to worry about the top hit for each...
#             my %ref_seqs;
#             $ref_seqs{'unaligned'} = [];
#             my $io = Bio::SearchIO->new( -format => 'blast', -file => 'clusters.blast' );
#             while ( my $result = $io->next_result() ) {
#                 my $top_hit = $result->next_hit();
#                 my $query   = $result->query_name();
#                 if ($top_hit) {
#                     my $ref_seq = $top_hit->name();
#                     if ( $ref_seqs{$ref_seq} ) {
#                         push( @{ $ref_seqs{$ref_seq} }, $query );
#                     }
#                     else {
#                         $ref_seqs{$ref_seq} = [$query];
#                     }
#                 }
#                 else {
#                     push( @{ $ref_seqs{'unaligned'} }, $query );
#                 }
#             }

#             # generate a fasta file of contigs which align to each reference
#             foreach my $ref ( keys(%ref_seqs) ) {
#                 my $out_name = "$run_dir/reference_${ref}_contigs";
#                 $out_name =~ s/lcl\|//;
#                 my $outIO = Bio::SeqIO->new( -format => 'fasta', -file => ">$out_name" );
#                 foreach my $id ( @{ $ref_seqs{$ref} } ) {
#                     my $contig = $fasta_db->get_Seq_by_id($id);
#                     $outIO->write_seq($contig);
#                 }
#             }

#         }
#         else {

#             # if we don't have mulitple references, we just need to make the
#             # reference and contigs available under consistent names
#             symlink( $reference, "$run_dir/reference_" . $ref_ids[0] )
#               or die "Error creating $run_dir/reference_${ref_ids[0]} symlink: $!"
#               if ( !-e "$run_dir/reference_" . $ref_ids[0] );

#             #my $src_contigs = "$run_dir/" . basename($contigs);
#             symlink( $contigs, "$run_dir/reference_${ref_ids[0]}_contigs" )
#               or die "Error creating $run_dir/reference_${ref_ids[0]}_contigs symlink";
#         }

#         my $mergedIO = Bio::SeqIO->new( -format => 'fasta', -file => ">$run_dir/scaffolds.fasta" );
#         my $merged_scaffolds;    #for final merging of per-reference scaffolds
#         my $merged_scaff_count = 0;

#         # Now run the selecting scaffolder on each set of reference and contigs...
#         foreach my $ref_id (@ref_ids) {

#             my $scaff_contigs = "reference_${ref_id}_contigs";

#             # Only if anything mapped to this reference
#             next unless ( -e "${run_dir}/${scaff_contigs}" );

#             print "\nScaffolding vs. $ref_id\n";

#             my $exec_cmd  = $cmd;
#             my $reference = "reference_${ref_id}";

#             my @replace = ( $cmd, $scaffold_output );
#             foreach ( $exec_cmd, $scaffold_output ) {
#                 s/__BUGBUILDER_BIN__/$FindBin::Bin/;
#                 s/__TMPDIR__/$tmpdir/g;
#                 s/__SCAFFDIR__/$run_dir\/${scaffolder}_${ref_id}/g;
#                 s/__RUN__/$run_id/;
#                 s/__FASTQ1__/$tmpdir\/read1.fastq/;
#                 s/__FASTQ2__/$tmpdir\/read2.fastq/;
#                 s/__REFERENCE__/${run_dir}\/${reference}/;
#                 s/__CONTIGS__/${run_dir}\/${scaff_contigs}/;
#                 s/__INSSIZE__/$insert_size/;
#                 s/__INSSD__/$insert_sd/;
#                 s/__THREADS__/$threads/;
#             }
#             if ( defined($unscaffolded_output) ) {
#                 $unscaffolded_output =~ s/__BUGBUILDER_BIN__/$FindBin::Bin/;
#                 $unscaffolded_output =~ s/__TMPDIR__/$tmpdir/g;
#                 $unscaffolded_output =~ s/__SCAFFDIR__/$run_dir\/$ref_id/g;
#                 $unscaffolded_output =~ s/__RUN__/$run_id/g;
#                 $unscaffolded_output =~ s/__FASTQ1__/$tmpdir\/read1.fastq/;
#                 $unscaffolded_output =~ s/__FASTQ2__/$tmpdir\/read2.fastq/;
#                 $unscaffolded_output =~ s/__REFERENCE__/${run_dir}_${reference}/;
#                 $unscaffolded_output =~ s/__CONTIGS__/$run_dir\/$scaff_contigs/;
#                 $unscaffolded_output =~ s/__INSSIZE__/$insert_size/;
#                 $unscaffolded_output =~ s/__INSSD__/$insert_sd/;
#             }
#             my $run_scaffold_output = "$run_dir/${scaffolder}_${ref_id}/$scaffold_output";

#             if ($scaffolder_args) {
#                 $exec_cmd .= "$scaffolder_args";
#             }
#             elsif ($default_args) {
#                 $exec_cmd .= "$default_args" if ($default_args);
#             }

#             $exec_cmd .= " 2>&1 >$tmpdir/${scaffolder}_${run_id}_${ref_id}.log";

#             mkdir("$run_dir/${scaffolder}_${ref_id}") or die "Error creating $run_dir/${scaffolder}_${ref_id}: $! ";
#             chdir("$run_dir/${scaffolder}_${ref_id}") or die "Error in chdir $run_dir/${scaffolder}_${ref_id}: $! ";
#             symlink( "$run_dir/${reference}", "$run_dir/${scaffolder}_${ref_id}/$reference" )
#               or die "Error linking $reference:$! ";

#             #symlink( "$run_dir/reference_${ref_ids[0]}_contigs", "$run_dir/${scaffolder}_${ref_id}/$contigs" )
#             #or die "Error linking $contigs:$! ";

#             system("perl $exec_cmd") == 0 or die "Error executing $exec_cmd: ";
#             my $scaffIO =
#               Bio::SeqIO->new( -format => 'fasta', -file => "$run_dir/${scaffolder}_${ref_id}/scaffolds.fasta" )
#               or die "Error opening $run_dir/${scaffolder}_${ref_id}/scaffolds.fasta: $!";
#             while ( my $scaff = $scaffIO->next_seq() ) {
#                 $scaff->display_id( 'scaffold_' . ++$merged_scaff_count );
#                 $mergedIO->write_seq($scaff);
#             }
#         }
#     }
#     else {    #non reference-guided scaffolding....

#         my $exec_cmd = $cmd;

#         # A kludge to work when tmpdir is not the top run dir - needed
#         # because align_reads concatenates path from tmpdir and contigs.fasta
#         if ( !-e "$tmpdir/contigs.fasta" ) {
#             symlink( $contigs, "$tmpdir/contigs.fasta" );
#         }

#         # if no reference provided we won't have an estimate of insert size,
#         # so need to get this by read alignment vs the assembly.
#         if ( !$insert_size ) {
#             align_reads( $tmpdir, "contigs.fasta", $mean_read_length, 1 );
#             ( $insert_size, $insert_sd ) = get_insert_stats( "$tmpdir", "contigs.fasta" );
#         }

#         my @replace = ( $cmd, $scaffold_output );
#         foreach ( $exec_cmd, $scaffold_output ) {
#             s/__BUGBUILDER_BIN__/$FindBin::Bin/;
#             s/__TMPDIR__/$tmpdir/g;
#             s/__SCAFFDIR__/$run_dir/g;
#             s/__RUN__/$run_id/;
#             s/__FASTQ1__/$tmpdir\/read1.fastq/;
#             s/__FASTQ2__/$tmpdir\/read2.fastq/;
#             s/__CONTIGS__/$contigs/;
#             s/__INSSIZE__/$insert_size/;
#             s/__INSSD__/$insert_sd/;
#             s/__THREADS__/$threads/;
#         }
#         if ( defined($unscaffolded_output) ) {
#             $unscaffolded_output =~ s/__BUGBUILDER_BIN__/$FindBin::Bin/;
#             $unscaffolded_output =~ s/__TMPDIR__/$tmpdir/g;
#             $unscaffolded_output =~ s/__SCAFFDIR__/$run_dir/g;
#             $unscaffolded_output =~ s/__RUN__/$run_id/g;
#             $unscaffolded_output =~ s/__FASTQ1__/$tmpdir\/read1.fastq/;
#             $unscaffolded_output =~ s/__FASTQ2__/$tmpdir\/read2.fastq/;
#             $unscaffolded_output =~ s/__CONTIGS__/$contigs/;
#             $unscaffolded_output =~ s/__INSSIZE__/$insert_size/;
#             $unscaffolded_output =~ s/__INSSD__/$insert_sd/;
#         }
#         my $run_scaffold_output = "$run_dir/${scaffolder}/$scaffold_output";

#         if ($scaffolder_args) {
#             $exec_cmd .= "$scaffolder_args";
#         }
#         elsif ($default_args) {
#             $exec_cmd .= "$default_args" if ($default_args);
#         }

#         $exec_cmd .= " 2>&1 >$tmpdir/${scaffolder}_${run_id}.log";

#         system("perl $exec_cmd") == 0 or die "Error executing $exec_cmd:$!";

#         #mkdir("$run_dir/${scaffolder}") or die "Error creating $run_dir/${scaffolder}: $! ";
#         #chdir("$run_dir/${scaffolder}") or die "Error in chdir $run_dir/{scaffolder}: $! ";
#         #symlink( "$run_dir/${reference}", "$run_dir/${scaffolder}/$reference" )
#         #  or die "Error linking $reference:$! ";

#         #symlink( "$run_dir/reference_${ref_ids[0]}_contigs", "$run_dir/${scaffolder}_${ref_id}/$contigs" )

#     }

#     # Create a fasta file of unplaced contigs if files of contig_ids are generated by the scaffolder wrapper
#     my @id_files = File::Find::Rule->file()->name("${scaffolder}.contig_ids")->in($run_dir);
#     if ( scalar(@id_files) ) {
#         my %used_contigs;
#         foreach my $id_file (@id_files) {
#             open CONTIG_IDS, $id_file or die "Error opening $id_file: $!";
#             while (<CONTIG_IDS>) {
#                 chomp;
#                 $used_contigs{$_}++;
#             }
#             close CONTIG_IDS;
#         }
#         my $inIO  = Bio::SeqIO->new( -format => 'fasta', -file => "$contigs" );
#         my $outIO = Bio::SeqIO->new( -format => 'fasta', -file => ">$run_dir/unplaced_contigs.fasta" );
#         while ( my $contig = $inIO->next_seq() ) {
#             $outIO->write_seq($contig) unless ( $used_contigs{ $contig->display_id() } && $contig->length() > 200 );
#         }
#     }

#     #renumber scaffolds to ensure they are unique...
#     my $count = 0;
#     my $inIO  = Bio::SeqIO->new( -file => "$run_dir/scaffolds.fasta", -format => "fasta" );
#     my $outIO = Bio::SeqIO->new( -file => ">$run_dir/scaffolds_renumbered.fasta", -format => "fasta" );

#     while ( my $seq = $inIO->next_seq() ) {
#         $seq->display_id( "scaffold_" . ++$count );
#         $outIO->write_seq($seq);
#     }

#     print "\nScaffolded assembly stats\n=========================\n\n";
#     get_contig_stats( "$run_dir/scaffolds.fasta", 'scaffolds' );

#     if ( $run_id == 1 ) {
#         chdir $tmpdir or die "Error chdiring to $tmpdir: $! ";
#         unlink "scaffolds.fasta"
#           or die "Error removing scaffolds.fasta: $! "
#           if ( -l "scaffolds.fasta" );
#         symlink( "$run_dir/scaffolds_renumbered.fasta", "scaffolds.fasta" )
#           or die "Error creating symlink: $! ";
#         if ( defined($unscaffolded_output) ) {
#             symlink( "$run_dir/$unscaffolded_output", "unplaced_contigs.fasta" )
#               or die "Error creating symlink: $!";
#         }
#         elsif ( -e "$tmpdir/${scaffolder}/unplaced_contigs.fasta" ) {
#             symlink( "${scaffolder}/unplaced_contigs.fasta", "unplaced_contigs.fasta" )
#               or die "Error creating symlink: $!";
#         }
#     }
#     return ($linkage_evidence);

# }

# ######################################################################
# #
# # build_agp
# #
# # Creates an AGP file from the scaffolds, while generating new
# # contig/scaffold outputs meeting ENA requirements (no consecutive runs
# # of >=10 N, minimum contig size of 200 bp) if running in 'submission'
# # mode, otherwise leaves short contigs and gaps<100bp intact. The
# # scaffold_type argument is used to determine the linkage evidence type
# # for scaffold gaps
# #
# # required parameters: $ (tmpdir)
# #		     : $ (organism description)
# #                    : $ (mode - submission or draft)
# #                    : $ (scaffold_type: align or mate_pair)
# #
# # returns            : $ (0)
# #
# ######################################################################

# sub build_agp {

#     my $tmpdir   = shift;
#     my $organism = shift;
#     my $mode     = shift;
#     my $evidence = shift;

#     die "Unknown evidence type: $evidence"
#       unless (    $evidence eq 'paired-ends'
#                || $evidence eq 'align_genus'
#                || $evidence eq 'align_xgenus' );

#     message("Creating AGP file...");

#     mkdir "$tmpdir/agp" or die "Error creating agp dir: $!";
#     chdir "$tmpdir/agp" or die "Error chdiring to $tmpdir/agp: $!";

#     open AGP, ">scaffolds.agp"
#       or die "Error opening scaffolds.agp for writing: $! ";
#     print AGP "##agp-version 2.0\n";
#     print AGP "#$organism\n";

#     my $scaffold_inIO  = Bio::SeqIO->new( -file => "$tmpdir/scaffolds.fasta", -format => 'fasta' );
#     my $scaffold_outIO = Bio::SeqIO->new( -file => ">scaffolds.fasta",        -format => "fasta" );
#     my $contig_outIO   = Bio::SeqIO->new( -file => ">contigs.fasta",          -format => 'fasta' );

#     my $contig_count = 0;
#     my %gaps;    #overall per-scaffold gaps to return....

#     while ( my $scaffold = $scaffold_inIO->next_seq() ) {

#         my $contig_start = 0;
#         my $scaffold_loc = 1;
#         my $scaffold_id  = $scaffold->display_id();
#         my $scaff_count  = $1 if ( $scaffold_id =~ /([0-9]+)$/ );
#         $scaffold_id = sprintf( "scaffold_%06s", $scaff_count );
#         my $scaffold_end = $scaffold->length();
#         my ( $contig_end, $gap_start, $gap_end );

#         my @bases = split( //, $scaffold->seq() );
#         my ( %contigs, @gaps );    #location tracking for included contigs/gaps

#         # this is kind of crude, but seems to work...
#       BASE: for ( my $b_count = 0 ; $b_count <= $#bases ; $b_count++ ) {
#             if ( ( $bases[$b_count] ne 'N' ) && ( !$gap_start ) ) {
#                 $contig_end = $b_count;
#                 next BASE;
#             }
#             elsif ( ( $bases[$b_count] eq 'N' ) & !($gap_start) ) {
#                 $gap_start = $b_count if ( !$gap_start );
#                 next BASE;
#             }
#             elsif ( ( $bases[$b_count] ne 'N' ) && ($gap_start) ) {

#                 #we have left the gap
#                 $gap_end = $b_count;
#                 my $gap_length = $gap_end - $gap_start;
#                 if ( ( $gap_length < 10 ) && ( $mode eq 'submission' ) ) {

#                     # skip gaps <10 bp which are acceptable by EMBL
#                     $gap_start = undef;
#                     next BASE;
#                 }
#                 elsif ( ( $mode eq 'draft' ) && ( $gap_length <= 1 ) ) {

#                     # we can leave single ambiguous bases alone
#                     $gap_start = undef;
#                     next BASE;
#                 }
#                 else {
#                     $gap_end = $b_count;

#                     # Output contigs only > 200 bp
#                     if (    ( ( $contig_end - $contig_start ) < 200 )
#                          && ( $mode eq 'submission' ) )
#                     {

#                         #need to extend previous gap to new position including short contig
#                         if ( scalar(@gaps) ) {
#                             my $last_gap = pop @gaps;
#                             my ( $last_start, $last_end ) = split( /-/, $last_gap );
#                             my $new_gap = $last_start . '-' . $gap_end;
#                             push @gaps, $new_gap;
#                             $gap_start    = undef;
#                             $contig_start = $b_count;
#                             next BASE;
#                         }
#                     }
#                     else {
#                         my $contig_id = sprintf( "contig_%06s", ++$contig_count );
#                         my $contig_seq =
#                           Bio::Seq->new( -display_id => $contig_id,
#                                          -seq        => join( '', @bases[ $contig_start .. $contig_end ] ) );
#                         $contig_outIO->write_seq($contig_seq);
#                         $contigs{$contig_id} = {
#                                                  'coords' => $contig_start . '-' . $contig_end,
#                                                  length   => $contig_seq->length()
#                                                };
#                         my $gap = $gap_start . '-' . $gap_end;
#                         push @gaps, $gap;

#                         $gap_start    = undef;
#                         $contig_start = $b_count;
#                     }
#                 }
#             }
#         }

#         # Output last contig from last contig_start position to scaffold end if it is longer than
#         # 200 bp and we are running in submission mode, otherwise remove the last gap to truncate the scaffold...
#         if (    ( ( $scaffold_end - $contig_start ) > 200 )
#              || ( $mode eq 'draft' )
#              || $contig_count == 0 )
#         {
#             my $contig_id = sprintf( "contig_%06s", ++$contig_count );
#             my $contig_seq =
#               Bio::Seq->new( -display_id => $contig_id,
#                              -seq        => join( '', @bases[ $contig_start .. $#bases ] ) );
#             $contig_outIO->write_seq($contig_seq);
#             $contigs{$contig_id} = {
#                                      'coords' => $contig_start . '-' . $scaffold_end,
#                                      length   => $contig_seq->length()
#                                    };
#         }
#         else {
#             pop @gaps;
#         }

#         $gaps{$scaffold_id} = \@gaps;

#         my $scaffold_part = 0;

#         #write AGP output and new scaffolds fasta file
#         unlink("contigs.fasta.index") if ( -e "contigs.fasta.index" );
#         my $contig_db = Bio::DB::Fasta->new("contigs.fasta");
#         my $scaffold_seq;

#         my @contig_ids = map { $_->[0] }
#           sort { $a->[1] <=> $b->[1] }
#           map { [ $_, /(\d+)$/ ] } keys(%contigs);

#         if ( $#contig_ids > -1 ) {
#             for ( my $i = 0 ; $i <= $#contig_ids ; $i++ ) {
#                 my $contig_id   = $contig_ids[$i];
#                 my $contig_data = $contigs{$contig_id};
#                 my $coords      = $contig_data->{'coords'};
#                 my $length      = $contig_data->{'length'};

#                 if ( $contig_id ne "" ) {
#                     my $contig = $contig_db->get_Seq_by_id($contig_id);

#                     my ( $contig_start, $contig_end ) = split( /-/, $coords );

#                     print AGP "$scaffold_id\t"
#                       . ( $contig_start + 1 ) . "\t"
#                       . ( $contig_end + 1 ) . "\t"
#                       . ++$scaffold_part
#                       . "\tW\t$contig_id\t1\t$length\t+\n";
#                     $scaffold_seq .= $contig->seq();
#                     if ( $i < $#contig_ids ) {

#                         my $gap = $gaps[$i];
#                         my ( $gap_start, $gap_end ) = split( /-/, $gap );
#                         my $gap_size = $gap_end - $gap_start;

#                         print AGP "$scaffold_id\t"
#                           . ( $gap_start + 1 )
#                           . "\t$gap_end\t"
#                           . ++$scaffold_part
#                           . "\tN\t$gap_size\tscaffold\tyes\t$evidence\n";
#                         $scaffold_seq .= 'N' x $gap_size;
#                     }
#                 }
#             }
#         }
#         if ($scaffold_seq) {
#             my $scaffold_seqobj = Bio::Seq->new( -display_id => $scaffold_id, -seq => $scaffold_seq );
#             $scaffold_outIO->write_seq($scaffold_seqobj);
#         }
#     }

#     close AGP     or warn "Error closring $tmpdir/scaffolds.agp: $!";
#     chdir $tmpdir or die "Error chdiring to $tmpdir: $!";

#     unlink "scaffolds.fasta" or die "Error removing scaffolds.fasta: $!";
#     unlink "contigs.fasta"   or die "Error removing contigs.fasta; $!";

#     symlink( "agp/scaffolds.agp", "scaffolds.agp" )
#       or die "Error creating scaffolds.agp symlink: $!";
#     symlink( "agp/scaffolds.fasta", "scaffolds.fasta" )
#       or die "Error creating scaffolds.fasta symlink: $!";
#     symlink( "agp/contigs.fasta", "contigs.fasta" )
#       or die "Error creating contigs.fasta symlink: $!";

#     return ( \%gaps );

# }

# ######################################################################
# #
# # run_prokka
# #
# # generates annotation on the assembly using prokka
# #
# # required parameters: $ (tmpdir)
# #		       $ (genus)
# #                      $ (species)
# #                      $ (strain)
# #                      $ (locustag)
# #                      $ (centre)
# #
# # returns            : $ (0)
# #
# ######################################################################

# sub run_prokka {

#     my $tmpdir   = shift;
#     my $genus    = shift;
#     my $species  = shift;
#     my $strain   = shift;
#     my $locustag = shift;
#     my $centre   = shift;

#     message("Starting PROKKA...");

#     #use scaffolds if we have them, otherwise contigs....
#     my $seqs;
#     if ( -e "$tmpdir/scaffolds.fasta" ) {
#         $seqs = "$tmpdir/scaffolds.fasta";
#     }
#     else {
#         $seqs = "$tmpdir/contigs.fasta";
#     }

#     my $type = $1 if ( $seqs =~ /\/(contigs|scaffolds).fasta/ );

#     my $cmd = "prokka --addgenes --outdir $tmpdir/prokka --prefix prokka ";
#     $cmd .= "--genus $genus " if ( $genus && ( $genus ne "unknown_genus" ) );
#     $cmd .= "--species $species "
#       if ( $species && ( $species ne "unknown_species" ) );
#     $cmd .= "--strain $strain "
#       if ( $strain && ( $strain ne "unknown_strain" ) );
#     $cmd .= "--locustag $locustag " if ($locustag);
#     $cmd .= "--centre $centre "     if ($centre);
#     $cmd .= " $seqs ";
#     $cmd .= " >prokka.log 2>&1";

#     system($cmd);
#     if ( $? != 0 ) { print "prokka exited with $?...\n" }

#     # my $inIO        = Bio::SeqIO->new( -file => "$tmpdir/prokka/prokka.gbf",     -format => 'genbank' );
#     my $inIO        = Bio::SeqIO->new( -file => "$tmpdir/prokka/prokka.gbk",     -format => 'genbank' );
#     my $embl_outIO  = Bio::SeqIO->new( -file => ">$tmpdir/prokka/prokka.embl",   -format => 'embl' );
#     my $fasta_outIO = Bio::SeqIO->new( -file => ">$tmpdir/prokka/${type}.fasta", -format => 'fasta' );

#     # something is losing the scaffold/contig naming so needs to be regenerated...
#     # at least ordering should be conserved
#     my $count = 0;
#     while ( my $seq = $inIO->next_seq() ) {
#         my $label;
#         ( $type eq 'scaffolds' ) ? ( $label = 'scaffold' ) : ( $label = 'contig' );
#         my $id = sprintf( "${label}_%06s", ++$count );
#         $seq->display_id($id);
#         $seq->accession($id);
#         $embl_outIO->write_seq($seq);
#         $fasta_outIO->write_seq($seq);
#     }

#     chdir $tmpdir or die "Error chdiring to $tmpdir: $!";
#     unlink "scaffolds.fasta" or die "Error removing scaffolds.fasta: $!" if ( $type eq 'scaffolds' );

#     symlink( "prokka/prokka.embl", "$type.embl" )
#       or die "Error creating $type.embl symlink: $!";
#     symlink( "prokka/scaffolds.fasta", "scaffolds.fasta" )
#       or die "Error creating scaffolds.fasta symlink: $!"
#       if ( $type eq 'scaffolds' );

#     return (0);
# }

# ######################################################################
# #
# # amosvalidate
# #
# # Uses abyss's abyss-samtoafg script to convert spades sam and contigs
# # to an amos bank
# #
# # required params: $ (tmpdir)
# #                  $ (insert size)
# #                  $ (insert size stddev)
# #
# # returns        : $ (0)
# #
# ######################################################################

# sub amosvalidate {

#     my $tmpdir        = shift;
#     my $insert_size   = shift;
#     my $insert_stddev = shift;

#     mkdir "$tmpdir/amos" or die "Error creating $tmpdir/amos: $!";
#     chdir "$tmpdir/amos" or die "Error running chdir $tmpdir/amos: $!";

#     my $seq_file;
#     ( -e "$tmpdir/scaffolds.fasta" ) ? ( $seq_file = "scaffolds.fasta" ) : ( $seq_file = "contigs.fasta" );

#     open SEQS, "$tmpdir/$seq_file"
#       or die "Error opening $tmpdir/$seq_file: $!";
#     open OUT, ">$tmpdir/amos/$seq_file"
#       or die "Error opening $tmpdir/amos/contigs.fasta: $!";

#     my $first_line = 0;    #kludgetastic...
#     while ( my $line = <SEQS> ) {
#         if ( $line =~ /^>/ ) {
#             if ( $first_line > 0 ) {
#                 print OUT "\n$line";
#             }
#             else {
#                 print OUT "$line";
#                 $first_line++;
#             }
#         }
#         else {
#             chomp $line;
#             print OUT $line;
#         }
#     }
#     close SEQS;
#     close OUT;

#     message("Converting to amos bank...");
#     my $cmd = $config->{'sam2afg'};
#     $cmd .= " -m $insert_size -s $insert_stddev " if ( $insert_size && $insert_stddev );
#     $cmd .= " $tmpdir/amos/${seq_file} $tmpdir/bwa/${seq_file}.sam";
#     $cmd .= " > amos.afg";

#     system($cmd) == 0 or die "Error executing $cmd: $!";

#     $cmd = $config->{'amos_dir'}
#       . "bank-transact -cb $tmpdir/amos/assembly.bnk -m $tmpdir/amos/amos.afg  > $tmpdir/amos/bank-transact.log 2>&1";
#     system($cmd) == 0 or die "Error executing $cmd: $!";

#     message("Running amosvalidate");

#     # read amosvalidate script and comment out '4xx' lines, which run SNP checks and are extreeeeemly slow....
#     open AMOSVALIDATE, $config->{'amos_dir'} . "/amosvalidate"
#       or die "Error opening " . $config->{'amos_dir'} . "/amosvalidate: $!";
#     open SCRIPT, ">$tmpdir/amos/amosvalidate" or die "Error opening $tmpdir/amos/amosvalidate: $!";
#     while ( my $line = <AMOSVALIDATE> ) {
#         $line =~ s/^4/#4/;
#         print SCRIPT $line;
#     }
#     close AMOSVALIDATE;
#     close SCRIPT;
#     chmod 0755, "$tmpdir/amos/amosvalidate" or die "Error running chmod $tmpdir/amos/amosvalidate: $!";

#     $cmd = "$tmpdir/amos/amosvalidate $tmpdir/amos/assembly.bnk > $tmpdir/amos/amosvalidate.log 2>&1";
#     system($cmd) == 0 or die "Error executing $cmd: $!";

#     chdir $tmpdir or die "Error chainging to $tmpdir: $!";
#     return (0);

# }

# ######################################################################
# #
# # find_origin
# #
# # Attempts to identify location of origin based upon contig overlapping
# # base 1 of the reference sequence. This assumes  each reference sequence
# # is a complete circular molecular i.e.a chromosome or a plasmid
# #
# # required params: $ (tmpdir)
# #                  $ (scaffolder)
# #                  $ (scaffolder_args)
# #                  $ (reference)
# #                  $ (insert_size)
# #                  $ (insert_stddev)
# #                  $ (mean_read_length)
# #
# # returns: $ (0)
# #
# ######################################################################

# sub find_origin {

#     my $tmpdir           = shift;
#     my $scaffolder       = shift;
#     my $scaffolder_args  = shift;
#     my $reference        = shift;
#     my $insert_size      = shift;
#     my $stddev           = shift;
#     my $mean_read_length = shift;
#     my $threads          = shift;

#     my $ori_dir = $tmpdir . "/origin";
#     mkdir "$ori_dir" or die "Error creating $ori_dir: $!";
#     chdir "$ori_dir" or die "Error changing to $ori_dir: $!";
#     symlink( "$tmpdir/read1.fastq", "read1.fastq" ) or die "Error creating read1.fastq symlink: $!";
#     symlink( "$tmpdir/read2.fastq", "read2.fastq" ) or die "Error creating read1.fastq symlink: $!";

#     message("Attempting to identify origin...");

#     my $cmd = $config->{'mummer_dir'}
#       . "/nucmer $tmpdir/reference_parsed_ids.fasta $tmpdir/scaffolds.fasta -p $ori_dir/ori > $ori_dir/nucmer.log 2>&1";
#     system($cmd) == 0 or die "Error executing $cmd: $!";
#     $cmd =
#       $config->{'mummer_dir'} . "/delta-filter -1 $ori_dir/ori.delta 2>$ori_dir/delta-filter.log > $ori_dir/ori.filter";
#     system($cmd) == 0 or die "Error executing $cmd: $!";
#     $cmd =
#       $config->{'mummer_dir'} . "/show-coords -H $ori_dir/ori.filter 2>$ori_dir/show-coords.log > $ori_dir/ori.coords";
#     system($cmd) == 0 or die "Error executing $cmd: $!";

#     open COORDS, "$ori_dir/ori.coords"
#       or die "Error opening ori.coords: $!";

#     my $origin;
#   LINE: while ( my $line = <COORDS> ) {
#         chomp $line;
#         $line =~ s/\|//g;
#         $line =~ s/^ *//;
#         my @fields = split( /\s+/, $line );
#         if ( !$origin && $fields[0] == 1 ) {
#             $origin = "$fields[8]:$fields[2]";
#             print "Potential origin found at $origin...\n\n";
#         }
#     }
#     close COORDS;

#     if ($origin) {
#         my $io      = Bio::SeqIO->new( -file => "$tmpdir/scaffolds.fasta",    -format => 'fasta' );
#         my $outIO   = Bio::SeqIO->new( -file => '>scaffolds_ori_split.fasta', -format => 'fasta' );
#         my $splitIO = Bio::SeqIO->new( -file => ">split_ori_scaff.fasta",     -format => 'fasta' );
#       SCAFFOLD: while ( my $scaffold = $io->next_seq() ) {
#             my ( $ori_scaffold, $pos ) = split( /:/, $origin );
#             if ( $ori_scaffold eq $scaffold->display_id() ) {
#                 my $part_a = $scaffold->subseq( 1, $pos );
#                 my $part_b = $scaffold->subseq( $pos + 1, $scaffold->length() );
#                 my $ori_a = Bio::Seq->new( -display_id => $scaffold->display_id . "_A", -seq => $part_a );
#                 my $ori_b = Bio::Seq->new( -display_id => $scaffold->display_id . "_B", -seq => $part_b );

#                 $splitIO->write_seq($ori_a);
#                 $splitIO->write_seq($ori_b);
#                 undef($splitIO);

#                 # Now rerun scaffolder on split sequence containing origin.
#                 run_scaffolder(
#                                 $ori_dir,          "$tmpdir/reference_parsed_ids.fasta",
#                                 $scaffolder,       $scaffolder_args,
#                                 $insert_size,      $stddev,
#                                 2,                 "$tmpdir/origin/split_ori_scaff.fasta",
#                                 $mean_read_length, $threads
#                               );
#                 my $scaffIO =
#                   Bio::SeqIO->new( -format => 'fasta', -file => "$ori_dir/${scaffolder}_2/scaffolds.fasta" );
#                 while ( my $scaffold = $scaffIO->next_seq() ) {
#                     $outIO->write_seq($scaffold);
#                 }
#                 next SCAFFOLD;
#             }
#             else {
#                 $outIO->write_seq($scaffold);
#             }
#         }

#         #renumber scaffolds to ensure they are unique...
#         my $count = 0;
#         my $inIO = Bio::SeqIO->new( -file => "$tmpdir/origin/scaffolds_ori_split.fasta", -format => "fasta" );
#         $outIO = Bio::SeqIO->new( -file => ">$tmpdir/origin/scaffolds_renumbered.fasta", -format => "fasta" );

#         while ( my $seq = $inIO->next_seq() ) {
#             $seq->display_id( "scaffold_" . ++$count );
#             $outIO->write_seq($seq);
#         }

#         chdir $tmpdir                     or die "Error changing to $tmpdir: $!";
#         unlink("$tmpdir/scaffolds.fasta") or die "Error unlinking $tmpdir/scaffolds.fasta:$!";
#         symlink( "$tmpdir/origin/scaffolds_renumbered.fasta", "$tmpdir/scaffolds.fasta" )
#           or die "Error symlinking scaffolds_ori_split.fasta:$!";
#     }

#     return (0);
# }

# ######################################################################
# #
# # order_scaffolds
# #
# # Identifies origin based on homology with reference.
# # Resulting scaffolds are then ordered and oriented relative to the reference...
# #
# # required params: $ (tmpdir)
# #                  $ (fasta reference)
# #
# # returns        : $ (0)
# #
# ######################################################################

# sub order_scaffolds {

#     my $tmpdir    = shift;
#     my $reference = shift;

#     mkdir "$tmpdir/orientating" or die "Error creating $tmpdir/orientating: $!";
#     chdir "$tmpdir/orientating"
#       or die "Error changing to $tmpdir/orientating: $!";

#     message("Orienting scaffolds vs. reference...");

#     my $cmd = $config->{'mummer_dir'}
#       . "/nucmer $tmpdir/$reference $tmpdir/scaffolds.fasta -p $tmpdir/orientating/ori2 > $tmpdir/orientating/nucmer2.log 2>&1";
#     system($cmd) == 0 or die "Error executing $cmd: $!";
#     $cmd = $config->{'mummer_dir'}
#       . "/delta-filter -1 $tmpdir/orientating/ori2.delta 2>$tmpdir/orientating/delta-filter2.log > $tmpdir/orientating/ori2.filter";
#     system($cmd) == 0 or die "Error executing $cmd: $!";
#     $cmd = $config->{'mummer_dir'}
#       . "/show-coords -H $tmpdir/orientating/ori2.filter 2>$tmpdir/orientating/show-coords2.log > $tmpdir/orientating/ori2.coords";
#     system($cmd) == 0 or die "Error executing $cmd: $!";

#     open COORDS, "$tmpdir/orientating/ori2.coords"
#       or die "Error opening ori.coords: $!";
#     my ( %orientations, $start, $end, $orient );
#     my %orient_count = ( '+' => 0, '-' => 0 );
#     my $contig = '';
#   LINE: while ( my $line = <COORDS> ) {
#         chomp $line;
#         $line =~ s/\|//g;
#         $line =~ s/^ *//;
#         my @fields = split( /\s+/, $line );

#         if ( ( $contig ne $fields[8] ) ) {
#             if ( $contig ne '' ) {    #end of previous contig
#                 store_orientation( \%orientations, $contig, \%orient_count );
#             }
#         }
#         $contig = $fields[8];
#         $start  = $fields[2];
#         $end    = $fields[3];
#         my $length;
#         if ( $start < $end ) {
#             $orient = '+';
#             $length = $end - $start;
#         }
#         else {
#             $orient = '-';
#             $length = $start - $end;
#         }
#         if ( $orient_count{$orient} ) {
#             $orient_count{$orient} = $orient_count{$orient} + $length;
#         }
#         else {
#             $orient_count{$orient} = $length;
#         }
#     }

#     close COORDS;

#     # record data for last contig...
#     store_orientation( \%orientations, $contig, \%orient_count );

#     my $orig_length     = 0;
#     my $oriented_length = 0;    #track how much sequence we align ok...
#     my $unplaced_length = 0;
#     my @unplaced;

#     my $io    = Bio::SeqIO->new( -file => '../scaffolds.fasta', -format => 'fasta' );
#     my $outIO = Bio::SeqIO->new( -file => '>scaffolds.fasta',   -format => 'fasta' );

#     # Reorientate contigs and break origin, rewriting to a new file...
#     my $i = 0;                  #for renumbering scaffolds...
#     while ( my $scaffold = $io->next_seq() ) {
#         $orig_length += $scaffold->length();

#         my $id = "scaffold_" . ++$i;
#         if ( $orientations{ $scaffold->display_id() } ) {
#             if ( $orientations{ $scaffold->display_id() } eq '+' ) {
#                 $scaffold->display_id($id);
#                 $outIO->write_seq($scaffold);
#                 $oriented_length += $scaffold->length();
#             }
#             else {
#                 my $rev = $scaffold->revcom();
#                 $rev->display_id($id);
#                 $outIO->write_seq($rev);
#                 $oriented_length += $scaffold->length();
#             }
#         }
#         else {
#             push @unplaced, $scaffold;
#         }

#     }

#     foreach my $scaffold (@unplaced) {
#         $outIO->write_seq($scaffold);
#         $unplaced_length += $scaffold->length();
#     }

#     my $ref_length;
#     my $refIO = Bio::SeqIO->new( -format => 'fasta', -file => "$tmpdir/$reference" );
#     while ( my $seq = $refIO->next_seq() ) {
#         $ref_length += $seq->length();
#     }

#     my $tb = Text::ASCIITable->new();
#     $tb->setCols( "", "Length (bp)" );
#     $tb->addRow( "Reference Sequence", $ref_length );
#     $tb->addRow( "Assembly",           $orig_length );
#     $tb->addRow( "Orientated contigs", $oriented_length );
#     $tb->addRow( "Unaligned contigs",  $unplaced_length );

#     print $tb, "\n";

#     chdir $tmpdir or die "Error changing to $tmpdir: $!";
#     unlink("scaffolds.fasta")
#       or die "Error removing scaffolds.fasta symlink";
#     symlink( "orientating/scaffolds.fasta", "scaffolds.fasta" )
#       or die "Error linking orientating/scaffolds.fasta: $!";

#     return (0);

# }

# ######################################################################
# #
# # store_orientation
# #
# # save correct contig orientation in hash passed as first parameter
# #
# # since alignment fof scaffold can contain blocks in reverse orientation
# # need to use the 'prevailing' orientation based on number of +/- blocks
# #
# # required parameters: $ (orientations hash)
# #                      $ (contig)
# #                      $ (hash of +/- base counts)
# #
# # returns           : none
# #
# ######################################################################

# sub store_orientation {

#     my $orientations = shift;
#     my $contig       = shift;
#     my $orient_count = shift;

#     if (
#          (
#               ( defined( $orient_count->{'+'} ) && defined( $orient_count->{'-'} ) )
#            && ( $orient_count->{'+'} > $orient_count->{'-'} )
#          )
#          || ( defined( $orient_count->{'+'} ) && !defined( $orient_count->{'-'} ) )
#        )
#     {
#         $orientations->{$contig} = '+';
#     }
#     else {
#         $orientations->{$contig} = '-';
#     }

# }

# ######################################################################
# #
# # run_varcaller
# #
# # Carries out variant calling using requested variant caller
# #
# # required params: $ (tmpdir)
# #                  $ (varcall)
# #                  $ (no. threads)
# #                  $ (read length)
# #
# # returns        : $ (0)
# #
# ######################################################################

# sub run_varcaller {

#     my $tmpdir      = shift;
#     my $varcall     = shift;
#     my $threads     = shift;
#     my $read_length = shift;

#     message("Running variant calling ($varcall)...");

#     my ( $cmd, $caller_cmd, $create_dir );
#     my $varcallers = $config->{'varcallers'};
#     foreach my $caller (@$varcallers) {
#         my $name = $caller->{'name'};
#         if ( $name eq $varcall ) {
#             $caller_cmd = $caller->{'command'};
#             $create_dir = $caller->{'create_dir'};
#         }
#     }
#     my $vardir = "$tmpdir/var_${varcall}/";
#     if ($create_dir) {
#         mkdir "$tmpdir/var_${varcall}" or die "Error creating $tmpdir/var_${varcall}: $! ";
#         chdir "$tmpdir/var_${varcall}" or die "Error chdiring to $tmpdir/var_${varcall}: $!";
#     }

#     symlink( "$tmpdir/reference_parsed_ids.fasta", "$vardir/reference.fasta" )
#       or die "Error creating $vardir/reference.fasta symlink: $!";
#     $cmd = $config->{'bwa_dir'} . symlink( "$tmpdir/read1.fastq", "$vardir/read1.fastq" )
#       or die "Error creating symlink: $! ";
#     symlink( "$tmpdir/read2.fastq", "$vardir/read2.fastq" )
#       if ( -e "$tmpdir/read2.fastq" )
#       or die "Error creating symlink: $! ";

#     print "BWA aligning reads to assembly...\n";
#     my $samtools_dir = $config->{'samtools_dir'};

#     $cmd = $config->{'bwa_dir'} . "/bwa index $vardir/reference.fasta >$vardir/bwa_index.log 2>&1";
#     system($cmd) == 0 or die " Error running $cmd";

#     # Use bwa-bwt for 'short' reads less than 100 bp, and bwa-mem for longer reads
#     if ( $read_length <= 100 ) {
#         $cmd =
#             $config->{'bwa_dir'}
#           . "/bwa aln -t $threads $vardir/reference.fasta $vardir/read1.fastq > $vardir/read1.sai"
#           . " 2> $vardir/bwa_sai1.log";
#         system($cmd) == 0 or die "Error running $cmd";
#         if ( -e "$vardir/read2.fastq" ) {
#             $cmd =
#                 $config->{'bwa_dir'}
#               . "/bwa aln -t $threads $vardir/reference.fasta $vardir/read2.fastq > $vardir/read2.sai"
#               . " 2> $vardir/bwa_sai2.log";
#             system($cmd) == 0 or die "Error running $cmd";
#             $cmd =
#                 $config->{'bwa_dir'}
#               . "/bwa sampe $vardir/reference.fasta $vardir/read1.sai $vardir/read2.sai "
#               . "$vardir/read1.fastq $vardir/read2.fastq";
#             $cmd .= " 2> $vardir/sampe.log > $vardir/scaffolds.sam";
#             system($cmd) == 0 or die "Error running $cmd";
#         }
#         else {
#             $cmd = $config->{'bwa_dir'} . "/bwa samse $vardir/reference.fasta $vardir/read1.sai $vardir/read1.fastq";
#             $cmd .= "2> $vardir/samse.log > $vardir/scaffolds.sam";
#             system($cmd) == 0 or die "Error running $cmd";
#         }
#     }
#     else {
#         if ( !-e "$vardir/read2.fastq" ) {

#             # single-ended long reads
#             $cmd =
#                 $config->{'bwa_dir'}
#               . "/bwa mem -t $threads -M $vardir/reference.fasta $vardir/read1.fastq > reference.sam "
#               . "2>$vardir/bwa_mem.log";
#             system($cmd) == 0 or die "Error running $cmd: $!";
#         }
#         else {

#             # paired-end long reads
#             $cmd =
#                 $config->{'bwa_dir'}
#               . "/bwa mem -t $threads -M $vardir/reference.fasta $vardir/read1.fastq $vardir/read2.fastq >reference.sam "
#               . "2>$vardir/bwa_mem.log";
#             system($cmd) == 0 or die "Error running $cmd: $!";
#         }
#     }

#     $cmd = "$samtools_dir/samtools view -q 10 -Sb $vardir/reference.sam 2>$vardir/samtoolsview.log"
#       . "|$samtools_dir/samtools sort - $vardir/reference";
#     system($cmd) == 0 or die "Error running $cmd";

#     $cmd = "$samtools_dir/samtools index $vardir/reference.bam 2>$vardir/samtools_index.log";
#     system($cmd) == 0 or die "Error running $cmd";

#     $caller_cmd =~ s/__BUGBUILDER_BIN__/$FindBin::Bin/;
#     $caller_cmd =~ s/__TMPDIR__/$tmpdir/;
#     $caller_cmd =~ s/__THREADS__/$threads/;

#     system($caller_cmd) == 0 or die "Error running $cmd: $!";
#     chdir $tmpdir            or die " Error chdiring to $tmpdir: $! ";

#     symlink( "var_${varcall}/var.filtered.vcf", "reference.variants.vcf" )
#       or die "Error creating $tmpdir/reference.variants.vcf: $!";
#     my $varcount = `grep -vc ^# $tmpdir/reference.variants.vcf`;
#     chomp $varcount;

#     print "\nIdentified $varcount variants...\n";

#     return (0);
# }

# ######################################################################
# #
# # build_comparisons
# #
# # generates a comparison appropriate for viewing with ACT and a
# # MUMmerplot to provide a quick overview
# #
# # required params: $ (tmpdir)
# #                  $ (reference)
# #                  $ (organism)
# #
# # returns        : $ (0)
# #
# ######################################################################

# sub build_comparisons {

#     my $tmpdir    = shift;
#     my $reference = shift;

#     mkdir "$tmpdir/comparisons" or die "Error creating $tmpdir/comparisons: $!";
#     chdir "$tmpdir/comparisons"
#       or die "Error changing to $tmpdir/comparisons: $!";
#     my $io     = Bio::SeqIO->new( -format => 'fasta', -file => "$tmpdir/$reference" );
#     my $ref    = $io->next_seq();
#     my $ref_id = $ref->display_id();
#     $ref_id = ( split( /\|/, $ref_id ) )[0] if ( $ref_id =~ /\|/ );

#     my $query;
#     if ( -l "$tmpdir/scaffolds.fasta" ) {
#         $query = "$tmpdir/scaffolds.fasta";
#     }
#     else {
#         $query = "$tmpdir/contigs.fasta";
#     }
#     my $cmd =
#         $config->{'blast_dir'}
#       . "/blastn  -query $query -subject $tmpdir/$reference "
#       . "-out $tmpdir/comparisons/comparison_vs_$ref_id.blastout -outfmt 6 > $tmpdir/comparisons/blast.log";
#     system($cmd) == 0 or die "Error running $cmd";

#     # also build a mummerplot in png format...
#     $cmd = $config->{'mummer_dir'} . "nucmer --prefix $ref_id $tmpdir/$reference $query >nucmer.log 2>&1";
#     system($cmd) == 0 or die "Error running $cmd";
#     $cmd = $config->{'mummer_dir'}
#       . "mummerplot -large --filter --layout -p $ref_id -t png -R $tmpdir/$reference -Q $query $ref_id.delta  >mummerplot.log 2>&1";
#     system($cmd) == 0 or die "Error running $cmd";

#     $cmd = $config->{'mummer_dir'}
#       . "mummerplot -large --filter --layout -c -p ${ref_id}_pip -t png -R $tmpdir/$reference -Q $query $ref_id.delta  >mummerplot_pip.log 2>&1";
#     system($cmd) == 0 or die "Error running $cmd";

#     chdir $tmpdir or die "Error chdiring to $tmpdir: $!";
#     symlink( "comparisons/comparison_vs_$ref_id.blastout", "comparison_vs_$ref_id.blastout" )
#       or die "Error creating symlink: $!";
#     symlink( "comparisons/$ref_id.png", "comparison_vs_$ref_id.png" )
#       or die "Error creating symlink: $!";
#     symlink( "comparisons/${ref_id}_pip.png", "comparison_vs_${ref_id}_pip.png" )
#       or die "Error creating symlink: $!";

#     return (0);

# }

# ######################################################################
# #
# # get_contig_to_iid_mapping
# #
# # generates a mapping of contig ids to amos IID
# #
# # required params: $ (tmpdir)
# #
# # returns        : $ (0)
# #
# ######################################################################

# sub get_contig_to_iid_mapping {

#     my $tmpdir = shift;
#     message("Extracting contig -> AMOS iid mapping...");

#     my $amos_dir = $config->{'amos_dir'};
#     my $mapping =
# `$amos_dir/bank-report -i -p -b $tmpdir/amos/assembly.bnk CTG 2> /dev/null|cut -f2,3 > $tmpdir/amos/ctg_to_iid.txt`;

#     return (0);
# }

# ######################################################################
# #
# # summarise_amosvalidate
# #
# # postprocesses amosvalidate outputs to make more readily digestible
# #
# # required params: $ (tmpdir)
# #
# # returns        : $ ($ - hashref to parsed results)
# #
# ######################################################################

# sub summarise_amosvalidate {

#     my $tmpdir = shift;

#     message("processing amosvalidate results...");

#     my $amosvalidate = $tmpdir . "/amos";
#     my $iid_mapping  = $tmpdir . "/amos/ctg_to_iid.txt";
#     my %results;

#     # Read Contig->iid mapping into a hash
#     my %contig_to_iid;
#     open IID, $iid_mapping or die "Error opening $iid_mapping: $!";
#     while (<IID>) {
#         my ( $contig, $iid ) = split(/\s+/);
#         chomp $iid;
#         $contig_to_iid{$iid} = $contig;
#         $results{$contig}    = [];
#     }
#     close IID;

#     opendir AMOSVALIDATE, $amosvalidate
#       or die "Error opening $amosvalidate: $!";
#     my @outputs = grep /feat$/, readdir AMOSVALIDATE;
#     close AMOSVALIDATE;

#     open ALL, "$amosvalidate/assembly.all.feat" or die "Could not open $amosvalidate/assembly.all.feat: $!";
#     while (<ALL>) {
#         my @fields    = split(/\s+/);
#         my $contig_id = $contig_to_iid{ $fields[0] };
#         my $start     = $fields[2];
#         my $end       = $fields[3];
#         my $type      = $fields[4];
#         my $res_arr   = $results{$contig_id};
#         $start = 0 if ( $start < 0 );    #no idea how it ends up with a negative start...but it does...
#         push @$res_arr, { 'start' => $start, 'end' => $end, type => $type, };
#         $results{$contig_id} = $res_arr;
#     }
#     close ALL;
#     return ( \%results );
# }
# ######################################################################
# #
# # merge_annotations
# #
# # updates annotated generated embl file with amosvalidate results
# # and gap locations
# #
# # required parameters: $ (tmpdir)
# #                      $ (hashref of amosvalidate results, keyed on contigid)
# #                      $ (hashref of scaffold gaps)
# #                      $ (genus)
# #                      $ (species)
# #                      $ (strain)
# #
# # returns            : $ (none)
# #
# ######################################################################

# sub merge_annotations {

#     my $tmpdir               = shift;
#     my $amosvalidate_results = shift;
#     my $gaps                 = shift;
#     my $genus                = shift;
#     my $species              = shift;
#     my $strain               = shift;

#     message("Merging annotations");

#     mkdir "$tmpdir/annotation_merge"
#       or die "Error creating $tmpdir/annotation_merge: $!";
#     chdir "$tmpdir/annotation_merge"
#       or die "Error chdiring to $tmpdir/annotation_merge: $!";
#     my $filename;

#     ( -e "$tmpdir/scaffolds.embl" ) ? ( $filename = "scaffolds.embl" ) : ( $filename = "contigs.embl" );

#     my $IO = Bio::SeqIO->new( -format => 'embl',
#                               -file   => "$tmpdir/$filename" );

#     my $outIO =
#       Bio::SeqIO->new( -format => 'embl',
#                        -file   => ">$filename" );

#     my %amos_colours = (
#                          'CE_STRETCH'      => '0 128 128',
#                          'CE_COMPRESS'     => '0 128 128',
#                          'HIGH_SNP'        => '128 128 0',
#                          'HIGH_READ_CVG'   => '255 0 0',
#                          'HIGH_KMER'       => '255 0 0',
#                          'KMER_COV'        => '255 0 0',
#                          'HIGH_OUTIE_CVG'  => '255 0 0',
#                          'HIGH_NORMAL_CVG' => '255 0 0',
#                          'LOW_GOOD_CVG'    => '0 0 255',
#                        );

#     my %amos_notes = (
#                        'CE_STRETCH' => 'Stretched mate-pairs: Possible repeat copy number expansion or other insertion',
#                        'CE_COMPRESS'     => 'Compressed mate-pairs; Possible collapsed repeat',
#                        'HIGH_SNP'        => 'High SNP frequency',
#                        'HIGH_READ_CVG'   => 'High read coverage; Possible collapsed repeat',
#                        'HIGH_KMER'       => 'High frequency of normalized kmers: Possible collapsed repeat',
#                        'KMER_COV'        => 'High frequency of normalized kmers: Possible collapsed repeat',
#                        'LOW_GOOD_CVG'    => 'Low coverage',
#                        'HIGH_NORMAL_CVG' => 'High coverage',
#                        'HIGH_OUTIE_CVG'  => 'High outie coverage',
#                      );

#     while ( my $embl_record = $IO->next_seq() ) {
#         my $orig_id = $embl_record->display_id();
#         $embl_record->display_id("$orig_id");
#         $embl_record->accession_number("$orig_id");
#         $embl_record->division('PRO');
#         $embl_record->molecule('genomic DNA');
#         $embl_record->is_circular(1);

#         $embl_record->add_date( get_embl_date() );
#         $embl_record->description("$genus $species $strain genome scaffold");

#         my @comments = $embl_record->annotation->get_Annotations('comment');
#         my $annot    = new Bio::Annotation::Collection;
#         my $comment  = Bio::Annotation::Comment->new;
#         $comment->text("Assembled using BugBuilder from http://github.com/jamesabbott/BugBuilder");
#         push @comments, $comment;
#         foreach my $c (@comments) {
#             $annot->add_Annotation( 'comment', $c );
#         }
#         $embl_record->annotation($annot);

#         # remove source entries from feature table
#         my @features = $embl_record->get_SeqFeatures();
#         $embl_record->flush_SeqFeatures();

#         # retrieve amosvalidate results for this contig and sort by start co-ordinate...
#         if ($amosvalidate_results) {
#             my @amos_features =
#               map  { $_->[0] }
#               sort { $a->[1] <=> $b->[1] }
#               map  { [ $_, $_->{'start'} ] } @{ $amosvalidate_results->{$orig_id} };

#             foreach my $feature (@amos_features) {
#                 if ( $feature->{'start'} != $feature->{'end'} ) {
#                     my $colour = $amos_colours{ $feature->{'type'} };
#                     my $note   = $amos_notes{ $feature->{'type'} };
#                     my $feature =
#                       new Bio::SeqFeature::Generic(
#                                                     -start       => $feature->{'start'} + 1,
#                                                     -end         => $feature->{'end'},
#                                                     -primary_tag => 'misc_feature',
#                                                     -tag         => {
#                                                               'note'   => $note,
#                                                               'colour' => $colour,
#                                                             }
#                                                   );
#                     push @features, $feature;
#                 }
#             }
#         }

#         if ($gaps) {
#             foreach my $scaffold ( keys(%$gaps) ) {
#                 if ( $orig_id eq $scaffold ) {
#                     my $scaffold_gaps = $gaps->{$scaffold};
#                     foreach my $gap (@$scaffold_gaps) {
#                         my ( $gap_start, $gap_end ) = split( /-/, $gap );
#                         my $est_length;
#                         if ( $gap_end - $gap_start == 100 ) {
#                             $est_length = 'unknown';
#                         }
#                         else {
#                             $est_length = $gap_end - $gap_start;
#                         }
#                         my $feature =
#                           new Bio::SeqFeature::Generic(
#                                           -start       => $gap_start,
#                                           -end         => $gap_end,
#                                           -primary_tag => 'assembly_gap',
#                                           -tag => { 'estimated_length' => $est_length, 'gap_type' => 'within_scaffold' }
#                           );
#                         push( @features, $feature );
#                     }
#                 }
#             }
#         }

#         my $source =
#           new Bio::SeqFeature::Generic(
#                                         -start       => 1,
#                                         -end         => $embl_record->length(),
#                                         -primary_tag => 'source',
#                                         -tag         => {
#                                                   'organism' => "$genus $species $strain",
#                                                   'strain'   => $strain
#                                                 }
#                                       );

#         $embl_record->add_SeqFeature($source);
#         @features = map { $_->[0] }
#           sort { $a->[1] <=> $b->[1] }
#           map { [ $_, $_->start() ] } @features;

#         foreach my $feat (@features) {
#             if ( $feat->primary_tag() ne 'source' ) {
#                 $embl_record->add_SeqFeature($feat);
#             }
#         }

#         $outIO->write_seq($embl_record);
#     }

#     chdir $tmpdir      or die "Error chdiring to $tmpdir: $!";
#     unlink "$filename" or die "Error unlinking $filename: $!";
#     symlink( "annotation_merge/$filename", "$filename" )
#       or die "Error creating $filename symlink: $!";

# }

# ######################################################################
# #
# # run_cgview
# #
# # Runs cgview to generate a genome map from the annotations
# #
# # Required parameters: $ (tmpdir)
# #
# # Returns: $ (0)
# #
# ######################################################################

# sub run_cgview {

#     my $tmpdir = shift;

#     my $cgview_dir  = $config->{'cgview_dir'};
#     my $xml_creator = $cgview_dir . "cgview_xml_builder/cgview_xml_builder.pl";
#     my $java        = $config->{'java'};

#     message("Creating genome visualisaton...");

#     chdir $tmpdir  or die "Error chdiring to $tmpdir: $!";
#     mkdir "cgview" or die "Error creating cgview directory: $!";
#     chdir "cgview" or die "Error chdiring to cgview $!";

#     my ( $embl, $outfile );
#     if ( -e "../scaffolds.embl" ) {
#         $embl    = "../scaffolds.embl";
#         $outfile = "scaffolds_cgview.png";
#     }
#     else {
#         $embl    = "../contigs.embl";
#         $outfile = "contigs_cgview.png";
#     }

#     my $cmd = "$xml_creator -sequence $embl -output scaffolds_cgview.xml -gc_skew T >xml_creator.log 2>&1";
#     system("$cmd") == 0 or die "Error running $cmd: $!";

#     $cmd = "${java} -jar ${cgview_dir}/cgview.jar -f png -i scaffolds_cgview.xml -o $outfile >cgview.log 2>&1";
#     system("$cmd") == 0 or die "Error running $cmd: $!";

#     chdir $tmpdir or die "Error chdiring to $tmpdir: $!";
#     symlink( "cgview/$outfile", "$outfile" ) or die "Error creating symlink: $!";

#     return (0);
# }

# ######################################################################
# #
# # get_contig_stats
# #
# # Reports contig statistics on assembly. Reports on scaffolds or
# # contigs depending upon 2nd argument passed - contigs gives values
# # for all contigs and those >200bp
# #
# # required params: $ (path to contigs)
# #                  $ ('scaffolds'|'contigs')
# #
# # returns          $ (0)
# #
# ######################################################################

# sub get_contig_stats {

#     my $file = shift;
#     my $type = shift;

#     my $IO = Bio::SeqIO->new( -format => 'fasta', -file => $file );
#     my ( %lengths, @all_lengths, $count, $tot_length, $max, $progress, $n50, $l50, $n50_tot );
#     my (
#          %lengths_200,  @all_lengths_200, $count_200, $tot_length_200, $max_200,
#          $progress_200, $n50_200,         $l50_200,   $n50_tot_200
#        );

#     while ( my $seq = $IO->next_seq() ) {
#         $lengths{ $seq->length() }++;
#         $count++;
#         $tot_length += $seq->length();
#         push @all_lengths, $seq->length();
#         if ( $seq->length() > 200 ) {
#             $lengths_200{ $seq->length() }++;
#             $count_200++;
#             $tot_length_200 += $seq->length();
#             push @all_lengths_200, $seq->length();
#         }
#     }

#     my $fifty     = $tot_length / 2;
#     my $fifty_200 = $tot_length_200 / 2;

#     my @sorted_lengths     = sort { $b <=> $a } @all_lengths;
#     my @sorted_lengths_200 = sort { $b <=> $a } @all_lengths_200;

#     # l50
#     foreach ( sort { $b <=> $a } keys(%lengths) ) {
#         $max = $_ if ( !$max );
#         $progress += $_ * ( $lengths{$_} );
#         if ( $progress >= $fifty ) {
#             $l50 = $_;
#             last;
#         }
#     }
#     foreach ( sort { $b <=> $a } keys(%lengths_200) ) {
#         $max_200 = $_ if ( !$max_200 );
#         $progress_200 += $_ * ( $lengths_200{$_} );
#         if ( $progress_200 >= $fifty_200 ) {
#             $l50_200 = $_;
#             last;
#         }
#     }

#     # n50
#     foreach my $length (@sorted_lengths) {
#         $n50_tot += $length;
#         $n50++;

#         # $l50 = $length;
#         last if ( $n50_tot >= $fifty );
#     }
#     foreach my $length (@sorted_lengths_200) {
#         $n50_tot_200 += $length;
#         $n50_200++;

#         # $l50_200 = $length;
#         last if ( $n50_tot_200 >= $fifty_200 );
#     }

#     my $tb = Text::ASCIITable->new();
#     if ( $type eq 'contigs' ) {
#         $type = ucfirst($type);
#         $tb->setCols( "", "All $type", "$type >200bp" );
#         $tb->addRow( "$type count",   $count,      $count_200 );
#         $tb->addRow( "Max Length",    $max,        $max_200 );
#         $tb->addRow( "Assembly size", $tot_length, $tot_length_200 );
#         $tb->addRow( "L50",           $l50,        $l50_200 );
#         $tb->addRow( "N50",           $n50,        $n50_200 );
#     }
#     else {
#         $type = ucfirst($type);
#         $tb->setCols( "", "All $type" );
#         $tb->addRow( "$type count",   $count );
#         $tb->addRow( "Max Length",    $max );
#         $tb->addRow( "Assembly size", $tot_length );
#         $tb->addRow( "L50",           $l50 );
#         $tb->addRow( "N50",           $n50 );
#     }
#     print $tb . "\n";

#     return (0);
# }

# ######################################################################
# #
# #  Pretty formats a status message
# #
# #  requried params: $ (message to display)
# #
# #  returns        : $ (0)
# #
# ######################################################################

# sub message {

#     my $message = shift;
#     print "\n\n", "*" x 80, "\n";
#     print "*",    " " x 78, "*\n";
#     print "* $message", " " x ( 77 - length($message) ), "*\n";
#     print "*", " " x 78, "*\n";
#     print "*" x 80, "\n\n";

# }

# ######################################################################
# #
# # get_embl_date
# #
# # returns the current date in EMBL style...
# #
# # required params: none
# #
# # returns: $ (formatted date)
# #
# ######################################################################

# sub get_embl_date {

#     my @months = qw( JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC );
#     my @vals   = ( localtime() )[ 3 .. 5 ];
#     my $date   = $vals[0] . '-' . $months[ $vals[1] ] . '-' . ( 1900 + $vals[2] );

#     return ($date);

# }

# ######################################################################
# #
# # parse_fasta_id
# #
# # attempt to parse an ID from the reference fasta file. Different fasta formats
# # make this tricky, so we will specifically parse plain IDs, ENA and NCBI
# # formats
# #
# # required params: $ (fasta file)
# #
# # returns: $ (hashref of parsed sequence IDs)=
# #
# ######################################################################

# sub parse_fasta_id {

#     my $file = shift;
#     my $IO = Bio::SeqIO->new( -format => 'fasta', -file => $file );
#     my @seq_ids;
#     while ( my $seq = $IO->next_seq() ) {
#         my $id = $seq->display_id();
#         if ( $id =~ /\|/ ) {
#             my @fields = split( /\|/, $id );
#             if ( $fields[0] eq 'ENA' ) {
#                 push( @seq_ids, $fields[1] );
#             }
#             elsif ( $fields[0] eq 'gi' ) {
#                 push( @seq_ids, $fields[3] );
#             }
#         }
#         else {
#             if ( $id =~ />([^ ])/ ) {
#                 push( @seq_ids, $1 );
#             }
#         }
#     }
#     return (@seq_ids);
# }

# ######################################################################
# #
# # parse_ref_id
# #
# # Similar to parse_fasta_id above, but works directly on a passed ID
# #
# # required params: $ (id)
# #
# # returns: $ (id)
# #
# ######################################################################

# sub parse_ref_id {

#     my $id = shift;
#     my $parsed_id;
#     if ( $id =~ /\|/ ) {
#         my @fields = split( /\|/, $id );
#         if ( $fields[0] eq 'ENA' ) {
#             $parsed_id = $fields[1];
#         }
#         elsif ( $fields[0] eq 'gi' ) {
#             $parsed_id = $fields[3];
#         }
#     }
#     else {
#         if ( $id =~ />([^ ])/ ) {
#             $parsed_id = $1;
#         }
#         else {
#             $parsed_id = $id;
#         }
#     }
#     return ($parsed_id);
# }

# ######################################################################
# #
# # show_tools
# #
# # Reports on configured assemblers, scaffolders and platforms
# #
# # Required params: $ (config hash)
# #
# # Returns:         $ ()
# #
# ######################################################################

# sub show_tools {

#     my $config = shift;

#     my @available_assemblers  = map { $_->{'name'} } @{ $config->{'assemblers'} };
#     my @available_scaffolders = map { $_->{'name'} } @{ $config->{'scaffolders'} };
#     my @available_mergers     = map { $_->{'name'} } @{ $config->{'merge_tools'} };
#     my @available_finishers   = map { $_->{'name'} } @{ $config->{'finishers'} };
#     my @platforms             = map { $_->{'platforms'} } @{ $config->{'assembler_categories'} };

#     my %available_platforms;
#     foreach my $platform (@platforms) {
#         foreach my $p (@$platform) {
#             $available_platforms{$p}++;
#         }
#     }
#     my @available_platforms = sort( keys(%available_platforms) );

#     print "\nWelcome to BugBuilder\n\n";
#     print "Available assemblers: " . join( ", ",                @available_assemblers ),
#       "\n" . "Available scaffolders: " . join( ", ",            @available_scaffolders ),
#       "\n" . "Available assembly merging tools: " . join( ", ", @available_mergers ),
#       "\n" . "Available finishing tools: " . join( ", ",        @available_finishers ),
#       "\n" . "Configured platforms: " . join( ", ", @available_platforms ), "\n\n";

#     return ();
# }
