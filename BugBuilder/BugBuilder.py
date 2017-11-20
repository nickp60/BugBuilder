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
     min_length: 200
     max_length: null
     command_pe: run_abyss --tmpdir __TMPDIR__ --fastq1 __FASTQ1__ --fastq2 __FASTQ2__ --read_length __READ_LENGTH__
     contig_output: __TMPDIR__/abyss/abyss-contigs.fa
     scaffold_output: __TMPDIR__/abyss/abyss-scaffolds.fa
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
     ref_required: 1
     linkage_evidence: align_genus
     command: run_sis --reference __REFERENCE__ --contigs __CONTIGS__ --tmpdir __TMPDIR__ --scaff_dir __SCAFFDIR__
     scaffold_output: scaffolds.fasta
     unscaffolded_output: unplaced_contigs.fasta
     default_args: null
     create_dir: 1
     priority: 2
   - name: mauve
     ref_required: 1
     linkage_evidence: align_genus
     command: run_mauve --reference __REFERENCE__ --run __RUN__ --contigs __CONTIGS__ --tmpdir __TMPDIR__ --scaff_dir __SCAFFDIR__
     default_args: null
     create_dir: 1
     priority: 1
     scaffold_output: scaffolds.fasta
   - name: sspace
     ref_required: 0
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
     priority: 2
   - name: abyss-sealer
     command: run_abyss-sealer --tmpdir __TMPDIR__ --encoding __ENCODING__ --threads __THREADS__
     create_dir: 1
     ref_required: 0
     paired_reads: 1
     priority: 3
   - name: pilon
     command: run_pilon --tmpdir __TMPDIR__ --threads __THREADS__
     create_dir: 1
     ref_required: 0
     paired_reads: 1
     priority: 1

varcallers:
   - name: pilon
     command: run_pilon_var --tmpdir __TMPDIR__ --threads __THREADS__
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
import hashlib
import re
import yaml
import sys
import arrow
import shutil
import logging
import statistics
import subprocess
import pkg_resources
import tabulate

from BugBuilder import __version__
from Bio import SeqIO # , SearchIO
from Bio.SeqRecord import SeqRecord
from argparse import Namespace
from .shared_methods import make_nucmer_delta_show_cmds
from .run_sis import run as run_sis
from .run_spades import run as run_spades
from .run_abyss import run as run_abyss
sub_mains = {
    # "canu": run_canu,
    "abyss": run_abyss,
    "spades": run_spades,
    "sis": run_sis}
def parse_available(thing, path=None):
    if path is None: # path var is to allow easier testing
        try:
            config = return_config(get_config_path(), force=False, hardfail=False, logger=None)
        except Exception as e:
            return []
    else:
        config = return_config(path, force=False, hardfail=False, logger=None)

    # get all the potential names from the config
    thing_list = [x['name'].lower() for x in getattr(config, thing)]
    # return all the names that have an executable available
    return [x for x in thing_list if getattr(config, x.replace("-", "_")) is not None]

def configure(config_path, hardfail=True):
    with open(config_path, 'w') as outfile:  # write header and config params
        for line in __config_data__:
            outfile.write(line)
    all_programs_dict = {
        "mandatory_programs": [
            'seqtk', 'samtools', 'picard', 'blastn', 'makeblastdb', "nucmer",
            'R', 'mummer', "delta-filter",  "show-coords",  'bwa'],
        "mandatory_python_programs": [],
        "opt_programs": [
            'fastqc', 'sickle', 'mauve', 'prokka', "minimus",
            'aragorn', 'prodigal', 'hmmer3', 'infernal',
            'rnammer', 'tbl2asn', 'abyss', 'celera', "abyss", "abyss-sealer",
            'gapfiller', 'sspace', 'asn2gb', 'amos' , 'masurca',
            'gfinisher', 'pilon', 'vcflib', 'cgview'],
        "opt_python_programs": [
            'spades', 'sis', 'multifasta', 'barrnap']
        }
    output_dict = {}
    for category, programs in all_programs_dict.items():
        if len(programs) == 0: # pyyaml will write out curly brackets otherwise
            continue
        program_dict = dict((k, None) for k in programs)
        for prog in programs:
            # append .py if we need to (like for spades, sis, etc)
            # cant modify prog cause we need to be able to access name without exention
            extension = ".py" if "python" in category else ""
            if shutil.which(prog + extension):
                # replace dashes with underscores in config names but not in prog
                output_dict[prog.replace("-", "_")] = shutil.which(prog + extension)
            else:
                output_dict[prog.replace("-", "_")] = None
                if "mandatory" in category:
                    if hardfail:
                        raise OSError("%s is a mandatory program; please install" % prog)
    with open(config_path, 'a') as outfile:  # write paths to exes
        yaml.dump(output_dict, outfile, default_flow_style=False)


def return_config(config_path, force=False, hardfail=True, logger=None):
    """returns config namespace object, configuring if needed
    """
    try: # if the config file exists, parse it
        config = parse_config(config_path)
    except yaml.YAMLError:
        force = True
    # if config file is missing or broken, recreate the config file
    if force or config.STATUS != "COMPLETE":
        if logger is not None:
            logger.info("(Re-)Configuring BugBuilder config file: %s", config_path)
        else:
            print("(Re-)Configuring BugBuilder config file: %s" % config_path)
        configure(config_path, hardfail=hardfail)
        config = parse_config(config_path)
    return config

class JustConfigure(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
        try:
            config_path = get_config_path()
        except:  # not sure what this exception is:  its not an IOError or an OSError
            resource_package = pkg_resources.Requirement.parse("BugBuilder")
            config_path = '/'.join(('BugBuilder','config_data', 'BugBuilder.yaml'))
        configure(config_path, hardfail=True)
        print("(Re-)Configured BugBuilder config file: %s" % config_path)
        sys.exit(1)
        # setattr(args, self.dest, values)

def get_args():  # pragma: no cover
    """
    """
    parser = argparse.ArgumentParser(
        prog="BugBuilder",
        description="Bugbuilder %s: " % __version__ +
        "Automated pipeline for assembly of draft quality " +
        "bacterial genomes with reference guided scaffolding and annotation." +
        "Please see accompanying userguide for full documentation",
        add_help=False)  # to allow for custom help
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("--platform", dest='platform', action="store",
                               help="Sequencing platform  used i.e. illumina, 454, iontorrent",
                               choices=["illumina", "454", "iontorrent"],
                               type=str, required=True)
    requiredNamed.add_argument("-o", "--output", dest='outdir', action="store",
                               help="output directory; " +
                               "default: %(default)s", default=os.getcwd(),
                               type=str, required=True)
    requiredNamed.add_argument("--assemblers", dest='assemblers', action="store",
                               help="Assembler(s) to run - may be specified twice" +
                               ", in which case the two assemblers will be run " +
                               "in parallel and the results merged using minimus." +
                               " If no assembler is specified, BugBuilder will " +
                               "try to select an appropriate assembler " +
                               "automatically",
                               choices=parse_available("assemblers"),
                               nargs="*", default=[],
                               type=str, required=True)

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--fastq1", dest='fastq1', action="store",
                          help="Path to first read of paired library, or fragment library",
                          type=str, required="--long-fastq" not in sys.argv)
    optional.add_argument("--fastq2", dest='fastq2', action="store",
                          help="Path to second read of paired library",
                          type=str)
    optional.add_argument("--untrimmed_fastq1", dest='untrimmed_fastq1', action="store",
                          # help="used to hold path of raw F reads for use " +
                          # "with mascura; dont set from command line",
                          type=str, help=argparse.SUPPRESS)
    optional.add_argument("--untrimmed_fastq2", dest='untrimmed_fastq2', action="store",
                          # help="used to hold path of raw R reads for use " +
                          # "with mascura; dont use from command line",
                          type=str, help=argparse.SUPPRESS)
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
    optional.add_argument("--assembler-args", dest='assembler_args',
                          action="store",
                          help="Any additional arguments to pass to the " +
                          "assembler. Default values are set in the " +
                          "'default_args' attribute of the configuration " +
                          "file. If running multiple assemblers, " +
                          "assembler_args should be specified twice, once " +
                          "for each assemler, in the same order than the " +
                          "assemblers are specified.",
                          nargs="*", default=[],
                          type=str)
    optional.add_argument("--scaffolder", dest='scaffolder', action="store",
                          help="scaffolder to use",
                          choices=parse_available("scaffolders"),
                          type=str)
    optional.add_argument("--scaffolder-args", dest='scaffolder_args',
                          action="store",
                          help="args to pass to the scaffolder, in single quotes",
                          type=str)
    optional.add_argument("--merge-method", dest='merge_method', action="store",
                          help="merge method to use",
                          choices=parse_available("merge_tools"),
                          type=str)
    optional.add_argument("--finisher", dest='finisher', action="store",
                          help="finisher to use",
                          choices=parse_available("finishers"),
                          nargs="*",
                          type=str)
    optional.add_argument("--varcaller", dest='varcaller', action="store",
                          help="varcaller to use",
                          choices=parse_available("varcallers"),
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
    optional.add_argument("--downsample", dest='downsample', action="store",
                          help="Downsample depth; set to 0 to skip " +
                          "downsampling. default is 100x ",
                          type=int, default=None)
    optional.add_argument("--species", dest='species', action="store",
                          help="species of your bug", default="unknown_species",
                          type=str)
    optional.add_argument("--genus", dest='genus', action="store",
                          help="genus of your bug", default="unknown_genus",
                          type=str)
    optional.add_argument("--strain", dest='strain', action="store",
                          help="strain of your bug", default="unknown_strain",
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
    optional.add_argument("--skip-trim", dest='skip_trim', action="store_true",
                          help="Quality threshold for trimming reads ", default=False)
    optional.add_argument("--trim-qv", dest='trim_qv', action="store",
                          help="quality threshold to trim  ", default=20,
                          type=int)
    optional.add_argument("--trim-length", dest='trim_length', action="store",
                          help="min read klength to retain;" +
                          "50 (25 for reads <50bp)", default=50,
                          type=int)
    optional.add_argument("--skip-split-origin", dest='skip_split_origin',
                          action="store_true",
                          help="split at origin ", default=False)
    optional.add_argument("--keepall", dest='keepall',
                          action="store_true",
                          help="keep all intermediate files ", default=False)
    optional.add_argument("--threads", dest='threads', action="store",
                          help="threads  ",
                          type=int, default=4)
    optional.add_argument("--out-dir", dest='out_dir', action="store",
                          help="dir for results",
                          type=str)
    optional.add_argument("--tmp-dir", dest='tmp_dir', action="store",
                          help="dir for results",
                          type=str)
    optional.add_argument("--already_assembled_dirs", dest='already_assembled_dirs',
                          action="store",
                          help="dir(s) with existing assembler results; must" +
                          " specify assembler(s), and lists much match " +
                          "length", default=[],
                          nargs="*", type=str)
     # These can be used instead of the --already-assembled_dirs arg
    optional.add_argument("--contigs",
                          dest='already_assembled_contigs',
                          action="store",
                          help="path to contig file(s) if assembler(s) " +
                          "have already been run; this helps save time when " +
                          "rerunning analyses. If running multiple assemble" +
                          "rs, assemblies should be specified twice, once " +
                          "for each assemler, in the same order than the " +
                          "assemblers are specified.",
                          nargs="*",  default=[],
                          type=str)
    optional.add_argument("--scaffolds",
                          dest='already_assembled_scaffolds',
                          action="store",
                          help="path to scafoold file(s) if  assembler(s) " +
                          "have already been run; this helps save time when " +
                          "rerunning analyses. If running multiple assemble" +
                          "rs, assemblies should be specified twice, once " +
                          "for each assemler, in the same order than the " +
                          "assemblers are specified.",
                          nargs="*",  default=[],
                          type=str)
    # optional.add_argument("--scratchdir", dest='scratchdir', action="store",
    #                       help="dir for results",
    #                       type=str)
    optional.add_argument("--configure", action=JustConfigure, type=None,
                          help="Reconfigure and exit ", nargs=0)
    optional.add_argument("--memory", dest='memory', action="store",
                          help="how much memory to use",
                          type=int, default=4)
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
    optional.add_argument('--version', action='version',
                          version='BugBuilder {}'.format(__version__))
    # args = parser.parse_args(sys.argv[2:])
    args = parser.parse_args()
    return args

def set_up_logging(verbosity, outfile, name):
    """
    Set up logging a la pyani, with
    a little help from:
    https://aykutakin.wordpress.com/2013/08/06/
        logging-to-console-and-file-in-python/
    requires logging, os, sys, time
    logs debug level to file, and [verbosity] level to stderr
    return a logger object
    """
    logger = logging.getLogger('root')
    if (verbosity * 10) not in range(10, 60, 10):
        raise ValueError('Invalid log level: %s' % verbosity)
    # logging.basicConfig(level=logging.DEBUG)
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to given verbosity
    console_err = logging.StreamHandler(sys.stderr)
    console_err.setLevel(level=(verbosity * 10))
    console_err_format = logging.Formatter(
        str("%(asctime)s " + "%(levelname)s" +" %(message)s"),
        "%H:%M:%S")
    console_err.setFormatter(console_err_format)
    # set some pretty colors, shorten names of loggers to keep lines aligned
    # logging.addLevelName(logging.DEBUG, "\u001b[30m%s\033[1;0m" % "..")
    # logging.addLevelName(logging.INFO,  "\u001b[32m%s\033[1;0m" % "--")
    # logging.addLevelName(logging.WARNING, "\u001b[33m%s\033[1;0m" % "!!")
    # logging.addLevelName(logging.ERROR, "\u001b[31m%s\033[1;0m" % "xx")
    # logging.addLevelName(logging.CRITICAL, "\u001b[31m%s\033[1;0m" % "XX")
    logging.addLevelName(logging.DEBUG, "..")
    logging.addLevelName(logging.INFO,  "--")
    logging.addLevelName(logging.WARNING, "!!")
    logging.addLevelName(logging.ERROR, "xx")
    logging.addLevelName(logging.CRITICAL,  "XX")
    logger.addHandler(console_err)
    # create debug file handler and set level to debug
    try:
        logfile_handler = logging.FileHandler(outfile, "w")
        logfile_handler.setLevel(logging.DEBUG)
        logfile_handler_formatter = \
            logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        logfile_handler.setFormatter(logfile_handler_formatter)
        logger.addHandler(logfile_handler)
    except:
        logger.error("Could not open {0} for logging".format(outfile))
        sys.exit(1)
    logger.debug("Initializing logger")
    logger.debug("logging at level {0}".format(verbosity))
    return logger


def md5(fname):
    """ return md5 of file
    from https://stackoverflow.com/questions/3431825
    """
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def match_assembler_args(args):
    """Ensure that both assembler and assemb_args are same length

    these will both be lists coming off of argparse; we need to ensure that
    they have the same number of items in list.

    Args:
       args (Namespace): The argparse namespace
    Raises:
        ValueError: if the length of the args is unequal
    Returns:
        list: [assembler_name, assembler_argument_string]
    """
    assembler_list = []
    if len(args.assemblers) == 0:
        return [None, None]
    if len(args.assembler_args) != 0:
        if len(args.assemblers) != len(args.assembler_args):
            raise ValueError("length of assemblers must equal " +
                             "length of assembler args")
        for i, v in enumerate(args.assemblers):
            assembler_list.append([v, args.assembler_args[i]])
    else:
        for i, v in enumerate(args.assemblers):
            assembler_list.append([v, None])
    return assembler_list


def check_files_present(args):
    """ ensure our files are here if they are needed

    Args:
        argparse namespace
    Return:
        None
    Raises:
        ValueError: files duplicated
        ValueError: file not found
        ValueError: reference missing in draft mode
    """
    for f in [args.fastq1, args.fastq2, args.reference]:
        if f is not None and not os.path.exists(f):
            raise ValueError("file %s does not exist!" % f)
    if args.fastq2 is not None and args.fastq1 == args.fastq2:
        raise ValueError("fastq1 and fastq2 are the same!")
    if args.mode == "draft" and args.reference is None:
        raise ValueError("Draft mode requiresa reference!")
    # if args.reference is not None:
    #      if not os.path.exists(args.reference):
    #          raise ValueError(
    #              "reference sequence %s does not exist!" % args.reference )


def fastq_needs_newname(args):
    """ read the first read of each file;  if ends with -1 or -2, return true
    """
    for fastq in [args.fastq1, args.fastq2]:
        with open(fastq, "r") as inf:
            for read in SeqIO.parse(inf, "fastq"):
                if read.id.endswith("-1") or read.id.endswith("-2"):
                    return True
                break
    return False


def rename_fastq_seqids(args):
    """ for each fastq of pair, replace the -[1|2] with /[1|2]
    """
    flist = ["fastq1", "fastq2"]
    for f in flist:
        # changed to allow perios
        pattern = re.compile("^([@\+][0-9A-Z\.\:\-]+)-([12])$")
        # this is almost exclusively for testing; should be checked prior
        fastq = getattr(args, f)
        if fastq is None:
            break
        tmp = os.path.splitext(fastq)[0] + "_renamed.fq"
        with open(tmp, "w") as outf:
            with open(fastq, "r") as inf:
                for line in inf:
                    line = pattern.sub('\\1/\\2', line)
                    outf.write(line)
        setattr(args, f, tmp)


def setup_tmp_dir(args, output_root, logger):
    """  setup_tmp_dir creates a temporary directory and copies
    over the relevent files
    Args:
        args (Namespace):
        output_root (str):
    Raises:
        Sys.Exit: temp dir cannot be created
    Returns:
        None
    """
    if args.tmp_dir is None:
        args.tmp_dir = os.path.join(output_root, "temp_dir")
    logger.info("Created working directory: %s" % args.tmp_dir)
    try:
        os.makedirs(args.tmp_dir, exist_ok=False)
    except FileExistsError: #OSError:
        logger.error("Temp directory already exists; exiting...")
        sys.exit(1)
    # copy to tmp dir
    for f in [args.fastq1, args.fastq2, args.long_fastq,
              args.de_fere_contigs]:
        if f is not None and os.path.exists(f):
            newpath = os.path.join(args.tmp_dir, os.path.basename(f))
            shutil.copyfile(f, newpath)
            # does this even work?
            f = newpath
     # we've seen some miseq fastq files have -1/-2 rather that /1 /2 pair
     # ids which cause problems with older software which doesn't expect this.
     # So, if you have PE data, and those read names are seen, we fix it
    if args.fastq1 is None and args.fastq2 is None:
        if fastq_needs_newname(args):
            rename_fastq_seqids(args)

    # Different tools have differing interpretations of various fasta ids,
    # so, rewrite the reference fasta file to ensure these just contain the id.
    # accepted formats are for ENA and Genbank format fasta headers, as well as plain IDs
    #  Or in this case, let BioPython complain
    if args.reference:
        new_reference = os.path.join(args.tmp_dir, os.path.basename(args.reference))
        with open(args.reference, "r") as inf:
            for rec in SeqIO.parse(inf, "fasta"):
                with open(new_reference, "a") as outf:
                    rec.desc = None
                    SeqIO.write(rec, outf, 'fasta')
        args.reference = new_reference


def return_open_fun(f):
    if os.path.splitext(f)[-1] in ['.gz', '.gzip']:
        open_fun = gzip.open
    else:
        open_fun = open
    return open_fun


def get_read_lens_from_fastq(fastq1, logger=None):
    """from LP; return total read length and count in fastq1 file from all reads
    """
    lengths = []
    open_fun = return_open_fun(fastq1)
    with open_fun(fastq1, "rt") as file_handle:
        data = SeqIO.parse(file_handle, "fastq")
        for read in data:
            lengths.append(len(read))
    tot = sum(lengths)
    if logger:
        logger.info("mean read length is {0} in {1}".format(
            os.path.basename(fastq1),
            float(tot / len(lengths))))
    file_handle.close()
    return (lengths)


def id_fastq_encoding(args, logger):
    """
    Identifies fastq quality encoding type, based upon the information
    at http://en.wikipedia.org/wiki/FASTQ_format:

    SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
    ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
    ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
    .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
    LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
    !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHoIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
    |                         |    |        |                              |                     |
    33                        59   64       73                            104                   126

    S - Sanger        Phred+33,  raw reads typically (0, 40)
    X - Solexa        Solexa+64, raw reads typically (-5, 40)
    I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
    J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
    with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
    L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)

    required params: $ (tmp_dir)

    returns:         $ (encoding)
    """
    min_val = 104;
    max_val = 0;
    if args.long_fastq is not None and os.path.exists(args.long_fastq):
        target = args.long_fastq
    else:
        target = args.fastq1
    open_fun = return_open_fun(target)
    with open_fun(target, "r") as inf:
        counter = 3
        for i, line in enumerate(inf):
            line = line.strip() #gets rid of eol, qual 10
            if i >= 1000:
                break
            if counter == 0:
                counter = 3
                quals = [ord(x) for x in list(line)]
                tmp_min =  min(quals)
                tmp_max =  max(quals)
                min_val = tmp_min if tmp_min < min_val else min_val
                max_val = tmp_max if tmp_max > max_val else max_val
            else:
                counter = counter - 1
                continue
    logger.info("Quality score Min: %d;  Max:  %d", min_val, max_val)
    if min_val <= 59 and max_val <= 74:
        encoding = 'sanger'
    elif min_val > 59 and min_val < 64 and max_val < 104:
        encoding = 'solexa'
    elif min_val > 64:
        encoding = 'illumina'
    else:
        logger.error("Error detecting fastq encoding!")
        raise ValueError
        encoding = 'illumina'
    logger.info("Read encoding detection summary:\n" + tabulate.tabulate(
        [["min:", min_val],["max", max_val], ["encoding" , encoding]]))
    return encoding


def assess_reads(args, config, platform, logger=None):
    """
    Determine some characteristics of the provided reads. Need to determine
    the mean read length, wheterh they are variable length or not, and the
    quality encoding.

    We'll also look at the predicted coverage if a reference genome is available

    required params: $ (tmp_dir)

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
    counts = []
    means = []
    stddevs = []
    long_mean = 0
    long_stddev = 0
    tot_length = 0
    long_tot_length = 0
    # get the lengths of all the short reads
    type_list = ["short", "short", "long"]
    for i, read in enumerate([args.fastq1, args.fastq2, args.long_fastq]):
        if  read is None or not os.path.exists(read):
            continue
        lengths = get_read_lens_from_fastq(read, logger=logger)
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
            long_stddev = statistics.stdev(lengths)
        else:
            counts.append(len(lengths))
            means.append(statistics.mean(lengths))
            stddevs.append(statistics.stdev(lengths))

    # ensure paired reads have same numebr of reads
    if len(counts) > 1:
        if counts[0] != counts[1]:
            raise ValueError("Paired reads must have same number of reads!")
    #determin mean read length and stddev....
    if len(means) > 1:
        mean = statistics.mean(means)
        stddev = statistics.stdev(stddevs)
    else:
        mean = means[0]
        stddev = stddevs[0]
    # set type of library
    lib_type = None
    for category in config.assembler_categories:
        if platform is not None:
            if platform in category['platforms']:
                lib_type = category['name']
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

    if args.genome_size == 0:
        logger.info("Infering genome size from reference")
        if args.reference is not None and os.path.exists(args.reference):
            with open(args.reference, "r") as inf:
                for rec in SeqIO.parse(inf, "fasta"):
                    args.genome_size = args.genome_size + len(rec.seq)
    logger.info("Genome Size: %d", args.genome_size)
    logger.info("Genome Length: %d", tot_length)
    coverage, long_coverage = None, None
    try:
        coverage = round(float(tot_length / args.genome_size), 3)
        long_coverage = round(float(long_tot_length / args.genome_size), 3)
    except:
        logger.warning("Cannot calculate coverage without a reference!")
    encoding = id_fastq_encoding(args, logger)
    return Namespace(lib_type=lib_type, encoding=encoding,
                     read_length_mean=mean,
                     read_length_stddev=stddev,
                     long_read_length_mean=long_mean,
                     long_read_length_stddev=long_stddev,
                     read_bases=tot_length,
                     # set these later
                     paired=True if args.fastq2 is not None else False,
                     insert_mean=None, insert_stddev=None,
                     coverage=coverage, long_read_coverage=long_coverage,
                     downsampled_coverage=None)


def check_ref_needed(args, lib_type):
    """ ensure either a reference or genome size is provided for long read assembly
    """
    if lib_type == "long" :
        if (args.genome_size is 0 and args.reference is None):
            raise ValueError("Please supply a genome size or reference " +
                             "sequence when running a long-read assembly")

def get_config_path():
    """ return path to where BugBuilder.yaml got installed
    """
    resource_package = pkg_resources.Requirement.parse("BugBuilder")
    config = '/'.join(('BugBuilder','config_data', 'BugBuilder.yaml'))
    return pkg_resources.resource_filename(resource_package, config)


def parse_config(config_file):
    """ Read the config file, make a namespace object out of it
    """
    with open(config_file, 'r') as stream:
        try:
            yamldata = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            with open(config_file, 'r') as stream2:
                for line in stream2:
                    print(line.strip())
            raise yaml.YAMLError("error parsing config file!")
    newns = argparse.Namespace()
    newns.__dict__ = yamldata
    return(newns)


def check_and_get_assemblers(args, config, reads_ns, logger):
    """
    """
    logger.info("Requested assemblers: %s", " ".join(args.assemblers))
    if len(args.assemblers) > 2:
        raise ValueError("A maximum of 2 assemblers can be requested")
    # sanity check assemblers where these have been manually requested...
    for assembler in args.assemblers:
        if shutil.which(assembler) is None and \
           shutil.which(assembler + ".py") is None:
            raise ValueError("%s is nt in PATH! exiting" % assembler)
    # If no assembler is requested, we need to select one based on the
    # 'assembly_type' from the configuration file....
    if len(args.assemblers) == 0:
        logger.info("No assemblers selected; we will choose based " +
                    " on library type: %s", reads_ns.lib_type)
        for category in config.assembler_categories:
            if category['name'] == reads_ns.lib_type:
                args.assemblers = category['assemblers']
    # then check the requested assemblers are appropriate for
    # the sequence characteristics...
    tools = []
    for assembler in args.assemblers:
        for conf_assembler in config.assemblers:
            if conf_assembler['name'] != assembler:
                continue
            if not reads_ns.paired and 'command_se' not in conf_assembler.keys():
                logger.warning(str("%s requires paired reads, but you only " +
                               "specified one fastq file. Removing from list" +
                                   " of assemblers..."), assembler)
                args.assemblers = [x for x in args.assemblers if x != assembler]
                if len(args.assemblers) == 0:
                    raise ValueError("No valid assemblers chosen!")
            if conf_assembler['min_length'] and reads_ns.read_length_mean < conf_assembler['min_length']:
                raise ValueError("%s does not support reads less than %d" % \
                                 (assembler, conf_assembler.min_length))
            if conf_assembler['max_length'] and reads_ns.read_length_mean > conf_assembler['max_length']:
                raise ValueError("%s does not support reads greather than %d" % \
                                 (assembler, conf_assembler['max_length']))
            if conf_assembler['insert_size_required'] and not reads_ns.insert_mean and not args.reference:
                raise ValueError(str("%s requires the library insert size and " +
                                 "stddev to be provided. Please add the " +
                                 "--insert-size and --insert-stddev " +
                                 "parameters, or provide a reference " +
                                 "sequence") % assembler)
            tools.append(conf_assembler)
    logger.info(str("Approved assembler(s) based on library type and " +
                    "availability: %s"), " ".join(args.assemblers))
    return tools


def get_scaffolder_and_linkage(args, config, paired, logger):
    if args.scaffolder is None:
        return None, None
    logger.debug(config.scaffolders)
    if args.scaffolder not in [k['name'].lower() for k in config.scaffolders]:
        raise ValueError("%s not an available scaffolder!" %args.scaffolder)
    linkage_evidence = None
    for conf_scaffolder in config.scaffolders:
        if conf_scaffolder['name'].lower() == args.scaffolder:
            if conf_scaffolder['linkage_evidence'] == "paired-ends" and not paired:
                    raise ValueError(str(
                        "%s requires paired reads, but you only specified " +
                        "one fastq file.") % args.scaffolder)
            elif "align" in conf_scaffolder['linkage_evidence'] and args.reference is None:
                raise ValueError(str("%s requires a reference for alignment, " +
                                     "but none is specified.") % args.scaffolder)
            else:
                pass
            return (conf_scaffolder, conf_scaffolder['linkage_evidence'])


def get_merger_tool(args, config, paired):
    """must be used before scaffolder getter
    """
    if args.merge_method is None:
        return None
    if args.merge_method not in [x.name for x in config.merge_tools]:
        raise ValueError("%s not an available merge_method!" % args.merge_method)
    linkage_evidence = None
    for merger_method in config.merge_methods:
        if merger_method['name'] == args.merge_method:
            if merge_methods['allow_scaffolding']:
                args.scaffolder = None
            return merge_method

def get_finisher(args, config, paired):
    if args.finisher is None:
        return None
    if args.finisher not in [x['name'] for x in config.finishers]:
        raise ValueError("%s not an available finisher!" %args.finisher)
    for conf_finisher in config.finishers:
        if conf_finisher['name'] == args.finisher:
            if args.reference is None and conf_finisher['ref_required']:
                raise ValueError("%s requires a reference." % \
                                     args.finisher)
            elif not paired and conf_finisher['paired_reads']:
                raise ValueError("%s requires paired reads" % args.finisher)
            else:
                pass
            return conf_finisher

def get_varcaller(args, config, paired):
    if args.varcaller is None:
        return None
    if args.varcaller not in [x['name'] for x in config.varcallers]:
        raise ValueError("%s not an available varcaller!" %args.varcaller)
    for conf_varcallers in config.varcallers:
        if conf_varcallers['name'] == args.varcaller:
            if args.reference is None and conf_varcallers['ref_required']:
                raise ValueError("%s requires reference, but you " +
                                 "only specified one fastq file." % \
                                 args.varcaller)
            return conf_varcallers


def select_tools(args, config, reads_ns, logger):
    """
    returns a namespace wit the tools from the config!
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
    assembler_tools = check_and_get_assemblers(
        args=args, config=config, reads_ns=reads_ns, logger=logger)
    logger.debug("get scaffolder")
    # merge tool can disqalify scaffolder, so it must be checked first
    merge_tool = get_merger_tool(args, config, paired=reads_ns.paired)
    scaffolder_tool, linkage_evidence = get_scaffolder_and_linkage(
        args=args, config=config, paired=reads_ns.paired, logger=logger)
    logger.debug("getting merger")
    finisher_tool = get_finisher(args, config, paired=reads_ns.paired)
    varcall_tool = get_varcaller(args, config, paired=reads_ns.paired)

    logger.info("linkage: %s" % linkage_evidence)
    return Namespace(assemblers=assembler_tools,
                     scaffolder=scaffolder_tool,
                     scaffold_type=linkage_evidence,
                     merge_method=args.merge_method,
                     finisher=finisher_tool,
                     varcall=varcall_tool)


def assembler_needs_downsampling(tools):
    for assembler in tools.assemblers:
        if assembler['downsample_reads']:
            return True
    return False


def make_fastqc_cmd(args, outdir):
    cmd = \
        "fastqc -t {4} --extract -o {0}{1}{2}{3} > {0}fastqc.log 2>&1".format(
            outdir,
            " " + args.fastq1 if args.fastq1 is not None else "",
            " " + args.fastq2 if args.fastq2 is not None else "",
            " " + args.long_fastq if args.long_fastq is not None else "",
            args.threads)
    return cmd

def run_fastqc(reads_ns, args, logger=None):
    """
    Carries out QC analysis using fastqc...

    required params: $ (tmpdir)
                 $ (type)

    returns: $ (0)
    """
    logger.info("Running FastQC...")
    fastqc_dir = os.path.join(args.tmp_dir, "fastqc", "")
    os.makedirs(fastqc_dir)
    fastqc_cmd = make_fastqc_cmd(args, outdir=fastqc_dir)
    logger.debug(fastqc_cmd)
    fastqc_res = subprocess.run(fastqc_cmd,
                                shell=sys.platform != "win32",
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                check=False)
    if fastqc_res.returncode != 0:
        logger.warning("Error running fastqc with command %s", fastqc_cmd)
        sys.exit(1)
    report = []
    report.append(["Status", "Metric", "File"])
    fails = 0
    for f in [args.fastq1, args.fastq2, args.long_fastq]:
        if f is None:
            continue
        name = os.path.splitext(os.path.basename(f))[0] + "_fastqc"
        report.append(["--", "--", "--"])
        with open(os.path.join(
                fastqc_dir, name, "summary.txt"), "r") as s:
            for line in s:
                if "FAIL" in line:
                    fails = fails + 1
                report.append(line.split("\t"))
        shutil.copyfile(os.path.join(fastqc_dir, name + ".html"),
                        os.path.join(args.tmp_dir, name + ".html"))
    logger.info("FastQC STATISTICS:\n" + tabulate.tabulate(report))
    if fails > 0:
        if reads_ns.lib_type in ["long", "hybrid", "de_fere"]:
            logger.info("NB: Reported quality issues from_fastq are normal " +
                        "when analysing PacBio sequence with FastQC...");
        else:
            logger.warning("Some QC tests indicate quality issues with this " +
                           "data.Please examine the fastqc outputs for " +
                           "these reads")

def make_sickle_cmd(args, reads_ns, paired, out_dir):
    if paired:
        cmd = str("sickle pe -f {0} -r {1} -t {2} -q {3} -l {4} -o " +
                  "{5}read1.fastq -p {5}read2.fastq -s {5}singles.fastq" +
                  "> {5}sickle.log").format(
                      args.fastq1,  # 0
                      args.fastq2,  # 1
                      reads_ns.encoding,  # 2
                      args.trim_qv,  #3
                      args.trim_length,  #4
                      out_dir)
    else:
        cmd = str("sickle se -f {0} -t {1} -q {2} -l {3} -o {4}read1.fastq " +
                  "> {4}sickle.log").format(
                      args.fastq1,  # 0
                      reads_ns.encoding,  # 1
                      args.trim_qv,  #2
                      args.trim_length,  #3
                      out_dir)
    return cmd

def quality_trim_reads(args, config, reads_ns, logger):
    """
    Quality trims reads using sickle

    required params: $ (tmpdir)
                 $ (encoding)

    returns          $ (0)
    """
    if reads_ns.encoding is None:
        logger.warning("Sequence quality encoding could not be determined;" +
                       "Sequence quality trimming will be skipped...")
    trim_dir = os.path.join(args.tmp_dir, "qc_trim", "")
    os.makedirs(trim_dir)
    trim_cmd = make_sickle_cmd(args, reads_ns, paired=reads_ns.paired, out_dir=trim_dir)
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
        if reads_ns.paired:
            args.fastq1 = os.path.join(trim_dir, "read1.fastq")
            args.fastq2 = os.path.join(trim_dir, "read2.fastq")
        else:
            args.fastq1 = os.path.join(trim_dir, "read1.fastq")
    kept = 0
    disc = 0
    report = []
    with open(os.path.join(trim_dir, "sickle.log"), "r") as inf:
        for line in inf:
            if line.strip() == "":
                continue
            splits = re.split(":\s", line.strip())
            report.append([splits[0], " ".join(splits[1:])])
            if "kept" in line:
                kept = kept + int(line.split(": ")[1].split("(")[0].strip())
            elif "discarded" in line:
                disc = disc + int(line.split(": ")[1].split("(")[0].strip())
            else:
                pass
    logger.info("Quality Trimming result:\n" + tabulate.tabulate(report))
    if float(disc / kept) * 100 > 10:
        logger.warning(">10\% of reads discarded during read trimming. "+
                       "Please examine the FastQC outputs...")


def make_seqtk_ds_cmd(args, reads_ns, new_coverage, outdir, config, logger):
    assert isinstance(new_coverage, int), "new_coverage must be an int"
    frac = float(new_coverage / reads_ns.coverage)
    cmd_list = []
    logger.info("downsampleing to %f X  to %f X (%f)",
                reads_ns.coverage, new_coverage , frac )
    for reads in [args.fastq1, args.fastq2]:
        if reads is not None and os.path.exists(reads):
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
    ds_dir = os.path.join(tmp_dir, "seqtk_dir", "")
    os.makedirs(ds_dir)
    ds_cmds = make_seqtk_ds_cmd(args=args, reads_ns=reads_ns, config=config,
                                new_coverage=new_cov, outdir=ds_dir, logger=logger)
    for cmd in ds_cmds:
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    args.fastq1 = os.path.join(ds_dir, os.path.basename(args.fastq1)) if args.fastq1 is not None else None
    args.fastq2 = os.path.join(ds_dir, os.path.basename(args.fastq2)) if args.fastq1 is not None else None
    return new_cov


def make_bwa_cmds(args, config, outdir, ref, reads_ns, fastq1, fastq2):
    """ map with bwa
    returns a list of cmds and the path to the mapping file, as a tuple
    """
    cmd_list = []
    index_cmd = "{0} index {1}  > {2}bwa_index.log 2>&1".format(
        config.bwa, ref, outdir)
    cmd_list.append(index_cmd)
    # Use bwa-bwt for 'short' reads less than 100 bp, and bwa-mem for longer reads
    if reads_ns.read_length_mean <= 100:
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


def make_samtools_cmds(config, sam, outdir, sorted_bam):
    # removed the -q 10
    convert = str("{0} view  -Shb {1} | {0} sort - >" +
                  "{3}").format(
                      config.samtools, sam, outdir, sorted_bam)
    index = str("{0} index {1} 2>{2}samtools_index.log").format(
                      config.samtools, sorted_bam, outdir)
    return [convert, index]


def align_reads(dirname, reads_ns,  downsample, args, config, logger):
    """
    Maps reads to reference using bwa

    required params: $ (tmp directory);
                 $ (reference);
                 $ (read length)
                 $ (flag to indicate downsampling required...)

               : $ (0)
    """
    align_dir = os.path.join(args.tmp_dir, "align_" + dirname, "")
    bwa_dir = os.path.join(align_dir, "bwa", "")
    samtools_dir = os.path.join(args.tmp_dir, "samtools" + dirname, "")
    seqtk_dir = os.path.join(args.tmp_dir, "seqtk" + dirname, "")
    for d in [align_dir, bwa_dir, samtools_dir, seqtk_dir]:
        os.makedirs(d)
    bwa_reference = os.path.join(bwa_dir, os.path.basename(args.reference))
    shutil.copyfile(args.reference, bwa_reference)
    if downsample:
        logger.info("Downsampling reads for insert-size estimation...")
        ds_cmds = make_seqtk_ds_cmd(args=args, reads_ns=reads_ns, config=config,
                                new_coverage=10, outdir=seqtk_dir, logger=logger)
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
    bwa_cmds, mapping_sam = make_bwa_cmds(args, config, outdir=bwa_dir, ref=bwa_reference,
                                      reads_ns=reads_ns, fastq1=fastq1,
                                      fastq2=fastq2)
    sorted_bam = os.path.splitext(mapping_sam)[0] + ".bam"
    samtools_cmds =  make_samtools_cmds(config, sam=mapping_sam,
                                        outdir=samtools_dir, sorted_bam=sorted_bam)
    for cmd in bwa_cmds + samtools_cmds:
        logger.debug(cmd)
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    return sorted_bam


def make_picard_stats_command(bam, config, picard_outdir):
    cmd = str("{0} CollectInsertSizeMetrics " +
              "INPUT={1} HISTOGRAM_FILE={2}insert_histogram.pdf " +
              "OUTPUT={2}insert_stats.txt QUIET=true VERBOSITY=ERROR " +
              "ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT > " +
              "{2}CollectInsertMetrics.log 2>&1 ").format(
                  config.picard, bam, picard_outdir)
    return(cmd, picard_outdir + "insert_stats.txt")


def get_insert_stats(bam, reads_ns, args, config, logger):
    """
    converts contigs.sam -> bam, sorts, indexes and generates
    insert stats with Picard

    required params: $ (tmp directory);
                      $ (reference);

    returns        : $ (insert size)
                   : $ (stddev)
    """
    logger.info("Collecting insert size statistics ")
    stat_dir = os.path.join(args.tmp_dir, "insert_stats")
    os.makedirs(stat_dir)
    pic_cmd, statfile = make_picard_stats_command(bam, config=config, picard_outdir=stat_dir)
    logger.debug(pic_cmd)
    subprocess.run(pic_cmd,
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    get_next = False
    with open(statfile, "r") as inf:
        for line in inf:
            if get_next:
                min_insert = line.split("\t")[2]
                max_insert = line.split("\t")[3]
                mean_insert = line.split("\t")[4]
                stddev_insert = line.split("\t")[5]
                break
            if "MEDIAN" in line:
                get_next=True
    return (mean_insert, stddev_insert)

def replace_placeholders(string, config=None, reads_ns=None, args=None, results=None):
    # key is the name of placholder, and value is a list [object, attr]
    # if object is None, just use attribute string
    #
    # none of the objects are mandatory
    replace_dict = {
        # assembler placeholders
        # "__BUGBUILDER_BIN__": [None, "/$FindBin::Bin/g"],
        # "__ASMDIR__": os.path.dirname(config.assemblers[assembler['name']]),
        "__MEMORY__": [args, "memory"],
        "__TMPDIR__": [args, "tmp_dir"],
        "__FASTQ1__": [args, "fastq1"],
        "__FASTQ2__": [args, "fastq2"],
        "__ORIG_FASTQ1__": [args, "untrimmed_fastq1"],
        "__ORIG_FASTQ2__": [args, "untrimmed_fastq2"],
        "__DE_FERE_CONTIGS__": [args, "de_fere_contigs"],
        "__LONGFASTQ__": [args, "long_fastq"],
        "__REFERENCE__": [args, "reference"],
        "__CATEGORY__": [reads_ns, "lib_type"],
        "__ENCODING__": [reads_ns, "encoding"],
        "__GENOME_SIZE__": [args, "genome_size"],
        "__PLATFORM__": [args, "platform"],
        "__READ_LENGTH__": [reads_ns, "read_length_mean"],
        "__INSSIZE__": [reads_ns, "insert_mean"],
        "__INSSD__": [reads_ns, "insert_stddev"],
        "__THREADS__": [args, "threads"],
        # scaffolder placeholders
        "__CONTIGS__": [results, "current_scaffolds"],
        "__SCAFFOLDS__": [results, "current_scaffolds"],
        "__SCAFFDIR__": [results, "current_scaffolds"]

    }
    for k, v in replace_dict.items():
        if v[0] is not None:
            try:
                string = string.replace(
                    k, str(getattr(v[0], v[1])) if v is not None else "")
            except: # so many ways this could go wrong: IndexError, KeyError
                string = string.replace(k, "")

    # for k,v in replace_dict.items():
    #     string = string.replace(k, str(v) if v is not None else "")
    if re.match(".*(__.*?__).*", string):
        raise ValueError("Unsubstituted placeholder found! \n%s" % string)
    return string


def get_assembler_cmds(assembler, assembler_args, args, config, reads_ns):
    # get the proper running command
    if reads_ns.lib_type == "hybrid":
        cmd = assembler['command_hybrid']
    elif reads_ns.lib_type == "de_fere":
        cmd = assembler['command_hybrid']
    elif args.fastq2 is not None:
        cmd = assembler['command_pe']
    else:
        cmd = assembler['command_se']
    # add args
    if assembler_args:
        cmd = cmd + " " + assembler_args
    elif assembler['default_args']:
        cmd = cmd + " " + assembler['default_args']
    # get output name
    contig_output = assembler['contig_output']
    scaffold_output = assembler['scaffold_output'] if assembler['scaffold_output'] else None
    CREATE = assembler['create_dir']
    if CREATE:
        os.makedirs(os.path.join(args.tmp_dir, assembler['name']))
    cmd = replace_placeholders(string=cmd, config=config, args=args,
                               reads_ns=reads_ns)
    # cmd = cmd + " > {0}/{1}.log 2>&1".format(args.tmp_dir, assembler['name'])
    main_function = cmd.split(" ")[0]
    cmd = " ".join(cmd.split(" ")[1:])
    return (main_function,
            cmd,
            replace_placeholders(contig_output, config, reads_ns, args),
            replace_placeholders(scaffold_output, config, reads_ns, args))


def standardize_fasta_output(infile, outfile, ctype):
    assert ctype in ['scaffolds', 'contigs'], "invalid contig type"
    count = 0
    with open(infile, "r") as inf:
        with open(outfile, "w") as outf:
            for rec in SeqIO.parse(inf, "fasta"):
                count = count + 1
                rec.id = "%s_%06d" % (ctype, count)
                rec.description = ""
                SeqIO.write(rec, outf, "fasta")


def get_L50_N50(lengths):
    lengths.sort()
    tot = sum(lengths)
    fifty = float(tot / 2)
    progress, N50 = 0, 0
    for L50 in lengths:
        N50 = N50 + 1
        progress = progress + L50
        if progress >= fifty:
            return (L50, N50)


def get_contig_stats(contigs, ctype):
    """
    Reports contig statistics on assembly. Reports on scaffolds or
    contigs depending upon 2nd argument passed - contigs gives values
    for all contigs and those >200bp

    required params: $ (path to contigs)
                 $ ('scaffolds'|'contigs')

    returns          $ (0)
    """
    assert ctype in ['scaffolds', 'contigs'], "invalid contig type"
    count, count_200 = 0, 0
    all_lengths, all_lengths_200 = [], []
    lengths, lengths_200 = [], []
    with open(contigs, "r") as inf:
        for rec in SeqIO.parse(contigs, "fasta"):
            length = len(rec.seq)
            count = count + 1
            all_lengths.append(length)
            if length > 200:
                lengths_200.append(length)
                count_200 = count_200 + 1
                all_lengths_200.append(length)

    sorted_lengths     = sorted(all_lengths)
    sorted_lengths_200 = sorted(all_lengths_200)
    L50, N50 = get_L50_N50(all_lengths)
    L50_200, N50_200 = get_L50_N50(all_lengths_200)
    return [
        ["", "All %s" % ctype, "%s >200bp" %ctype],
        [ctype + " count", count, count_200],
        ["Max Length", max(all_lengths), max(all_lengths_200)],
        ["Assembly size", sum(all_lengths), sum(all_lengths_200)],
        ["L50", L50, L50_200],
        ["N50", N50, N50_200 ]
    ]


def check_spades_kmers(assembler, cmd, readlen, min_diff=2, logger=None):
    """ warns and fixes bad kmers
    spades = "k=21,23,55"
    Just worry about SPAdes
    ABySS = "K=15 k=64", where K is span and k is kmer size, in de bruijn mode?
    canu = Picks kmers automatically
    """
    if assembler != "spades":
        return cmd
    prek_cmd = cmd.split("-k")[0].strip()
    k = cmd.split("-k")[1].strip().split(" ")[0].strip()
    postk_cmd = " ".join(cmd.split("-k")[1].split(" ")[2:])
    if k is 'auto':
        logger.debug("no need to check k, we let spades set k")
        return cmd
    try:
        klist = [int(x) for x in k.split(",")]
    except Exception as e:
        logger.error("error splitting kmers in %s by comma!", k)
        logger.error(e)
        raise ValueError
    logger.debug(klist)
    new_ks = []
    for i in klist:
        if i > readlen:
            logger.info("removing %d from list of kmers: exceeds read length",
                        i)
        elif readlen - i <= min_diff:
            logger.info("removing %d from list of kmers: too close " +
                        "to read length", i)
        elif i % 2 == 0:
            logger.info("removing %d from list of kmers: must be odd", i)
        else:
            new_ks.append(i)
    return "{0} -k {1} {2}".format(
        prek_cmd,
        ",".join([str(x) for x in new_ks]),
        postk_cmd)


def run_assembler(assembler, assembler_args, args, reads_ns, config, logger):
    """
    Runs specified assembler on fastq files

    required params: $ (tmpdir)
                 $ (assembler name)
                 $ (reference fasta sequece - probably not needed)
                 $ (arguments to pass to assembler)
                 $ (link - flag to indicate contigs should be symlinked into tmpdir)
                 $ (category - type of assembly to run, since an assembler may fall into more
                 than one category)
                 $ (encoding - some assemblers need explicitly telling)
                 $ (genome size)
                 $ (average read length)
	           $ (insert size)
                 $ (insert size stddev)
    returns          $ (0)
    """
    conf_assembler = [x for x in config.assemblers if x['name'] == assembler]
    assert len(conf_assembler) == 1, "multiple matches for assemblers"
    conf_assembler = conf_assembler[0]
    func, assembler_cmd, contigs_path, scaffolds_path = get_assembler_cmds(
        assembler=conf_assembler,
        assembler_args=assembler_args,
        args=args,
        config=config, reads_ns=reads_ns)
    logger.info(" Starting %s assembly ... ", assembler)
    # here we execute the "run" function from the appropriate runner script
    sub_mains[assembler](getattr(config, assembler), assembler_cmd, logger)


    logger.info("%s assembly statistics", assembler);
    report = get_contig_stats( contigs_path, 'contigs' )
    logger.info("CONTIG STATS:\n" + tabulate.tabulate(report))
    # Only generate scaffold stats for paired read alignments...
    if os.path.exists(scaffolds_path) and args.fastq2 is not None:
        report = get_contig_stats( scaffolds_path, 'scaffolds' )
        logger.info("SCAFFOLD STATS:\n" + tabulate.tabulate(report))
    else:
        logger.info("SCAFFOLD STATS: None")

    # rename contigs/scaffolds for consistent naming, since we need to retrieve by
    # id later, so it helps if we know what the ids look like...
    # if ( !$create ) {
    #     my $outdir = dirname($contig_output);
    #     chdir $outdir or die " Error changing to dir $outdir: $! ";
    # }
    renamed_contigs = os.path.join(os.path.dirname(contigs_path),
                                   "BugBuilder.contigs.fasta")
    standardize_fasta_output(
        infile=contigs_path,
        outfile=renamed_contigs,
        ctype="contigs")
    if os.path.exists(scaffolds_path):
        renamed_scaffolds = os.path.join(os.path.dirname(scaffolds_path),
                                           "BugBuilder.scaffolds.fasta")
        standardize_fasta_output(
            infile=scaffolds_path,
            outfile=renamed_scaffolds,
            ctype="contigs")
    else:
        renamed_scaffolds = None
    return(renamed_contigs, renamed_scaffolds)


def merge_assemblies(args, config, reads_ns, logger):
    """
    combines two assemblies using the selected merge-method

    required params: $ (tmpdir)
                 $ (arrayref of assemblers used)
                 $ (method)
                 $ (reference)

    returns        : $ (0)
    """
    logger.info("Merging assemblies (%s)...", args.merge_tool)
    try:
        tool = getattr(config, args.merge_tool)
    except AttributeError:
        raise ValueError("%s merge tool not available" % args.merge_tool)

    cmd = tool['command']
    contig_output = tool['contig_output']
    if tool['create_dir']:
        outdir = os.path.join(args.tmpdir, args.merge_tool)
        os.makedirs(outdir)
    merge_cmd = replace_placeholders(string=cmd, config=config,
                                     reads_ns=reads_ns, args=args)
    logger.debug(merge_cmd)
    subprocess.run(merge_cmd,
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)

    # symlink( "$method/$contig_output", "contigs.fasta" )
    #   or die " Error creating symlink : $! ";

    report = get_contig_stats(contig_output, 'contigs' )
    logger.info("MERGED ASSEMBLY STATISTICS:\n" + tabulate.tabulate(report))
    return contig_output


def check_id(args, contigs, config, logger):
    """
    Checks wether the assembled contigs have sufficient identity to the
    reference for reference-based scaffolding

    required params: $ (tmp directory);

    returns        : $ id_ok (boolean)
    """
    logger.info("Checking identity of assembly with reference...")
    ID_OK = False
    check_dir = os.path.join(args.tmp_dir, "id_check", "")
    os.makedirs(check_dir)
    cmd = str("{0} -query {1} -subject {2} -outfmt 6 -evalue 0.01 -out " +
              "{3}blastout.tab 2>&1 > {3}blastn.log").format(
        config.blastn, contigs, args.reference, check_dir)
    logger.debug(cmd)
    subprocess.run(cmd,
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    fake_aves = []
    #TODO correct this
    with open(check_dir + "blastout.tab", "r") as inf:
        for line in inf:
            fake_aves.append(float(line.strip().split("\t")[2]))
    try:
        percent_id = statistics.mean(fake_aves)
    except statistics.StatisticsError: # usually a "mean requires at least one data point"
        percent_id = 0.0
    logger.debug("percent id between reference and contigs: %d", percent_id)
    # system($cmd) == 0 or die "Error running $cmd: $!";
    # my $blio = Bio::SearchIO->new( -format => 'blastxml',
    #                                -file   => 'blastout.xml' );
    # my ( $aligned, $unaligned );
    # while ( my $result = $blio->next_result() ) {
    #     foreach my $hit ( $result->hits() ) {
    #         my $tiling = Bio::Search::Tiling::MapTiling->new($hit);
    #         $aligned   += $tiling->num_aligned();
    #         $unaligned += $tiling->num_unaligned();
    #     }
    # }
    # my $percent_id = sprintf( '%.2f', ( $aligned / ( $aligned + $unaligned ) * 100 ) );
    # print "ID=$percent_id %\n";
    if percent_id > 80:
        ID_OK = True
    else:
        logger.warning("Percentage identity of assembly with reference looks" +
                       " too low (%.2f)", percent_id)
        logger.warning("Reference will not be used for scaffolding or " +
                       "ordering contigs")
    return ID_OK


def check_already_assembled_dirs(args, config, logger):
    """Ensure that the directory(ies) providded have the expect output files
    Chances are, they do, but we need to be exceedingly sure that they do
    """
    res_dict = {"contig": [],
                "scaffold": []}
    for path in args.already_assembled_dirs:
        if not os.path.exists(path):
            raise FileNotFoundError("Could not find directory: %s" % path)
    if len(args.assemblers) != len(args.already_assembled_dirs):
        raise ValueError(
            "You must explicity specify the assemblers used to generate " +
            "the already-assembled results.  The lengths of those lists do " +
            "not match")
    outlist = []
    for i in range(0, len(args.assemblers)):
        assembler = args.assemblers[i]
        for f in ["contig", "scaffold"]:
            # get teh configuration for the selected assembler
            assembler_conf = [x for x in config.assemblers if x['name'] == assembler][0]
            res = assembler_conf[f + '_output']
            # res = res.replace("__TMPDIR__", args.already_assembled_dirs[i])
            res =  os.path.join(args.already_assembled_dirs[i],
                                os.path.basename(res))
            if not os.path.exists(res):
                raise FileNotFoundError(str(
                    "Expected %s output %s could not be found. Please " +
                    "check files and  re-assemble") %(args.assemblers[i], res))
            else:
                res_dict[f].append(res)
        outlist.append(
            {"name": assembler,
             "contigs": res_dict["contig"][i],
             "scaffolds": res_dict["scaffold"][i]})
    return outlist


def check_already_assembled_args(args, config, logger):
    """Ensure that the args (*_contigs or *_scaffolds) provided have the
    expect output files retuns a list of {assembler, contigs, scaffolds} dictionaries
    """
    res_dict = {"contig": [],
                "scaffold": []}
    for path in [args.already_assembled_contigs, args.already_assembled_scaffolds]:
        if not os.path.exists(path):
            raise FileNotFoundError("Could not find directory: %s" % path)
    if len(args.assemblers) != len(args.already_assembled_contigs) or\
        len(args.assemblers) != len(args.already_assembled_scaffolds):
        raise ValueError(
            "You must explicity specify the assemblers used to generate " +
            "the already-assembled results.  The lengths of those lists do " +
            "not match")
    outlist = []
    for i in range(0, len(args.assemblers)):
        assembler = args.assemblers[i]
        outlist.append(
            {"name": assembler,
             "contigs": args.already_assembled_contigs[i],
             "scaffolds": args.already_assembled_scaffolds[i]})
    return outlist

def check_args(args, config):
    # ensure reference-requiring things are accounted for
    if args.merge_method is None and len(args.assemblers) > 1:
        raise ValueError("Must provide merge method if using multiple assemblers")
    for thing in ["scaffolder", "finisher"]:
        if getattr(config, thing + "s") is None:
            continue
        for tool in getattr(config, thing + "s"):
            if tool['name'].lower() == getattr(args, thing):
                if args.reference is None and tool['ref_required']:
                    raise ValueError("%s %s requires a reference." % \
                                     (thing, getattr(args, thing)))
    # ensure exes are there for fastqc, seqtk
    if args.downsample != 0 and config.seqtk is None:
        raise ValueError("Seqtk is needed for downsampling")
    if not args.skip_fastqc != 0 and config.fastqc is None:
        raise ValueError("fastqc not found; must enable --skip-fastqc")
    if not args.skip_fastqc and config.fastqc is None:
        raise ValueError("fastqc not found; must enable --skip-fastqc")
    if not args.skip_trim  and config.sickle is None:
        raise ValueError("sickle-trim not found; must enable --skip-trim")


def make_empty_results_object():
    """ returns a namepace containing both run history and paths of res files
    This is to keep the args object from being mutilated over time, and to
    Ensure that we have a way of looking back at a run to see what happened.
    """
    results = Namespace(
        organism=None,
        assemblers_list=None, # to be filled by match_assembler args
        assemblers_results_dict_list=[], # {name:None, contigs:None, scaffold:None}
        ID_OK=False, # to be set by check_id
        current_contigs=None, current_scaffolds=None, current_reference=None,
        current_contigs_source=None, current_scaffolds_source=None,
        old_contigs=[[None, None]],  # [path, source]
        old_scaffolds=[[None, None]],  # [path, source]
        old_references=[[None, None]],  # [path, source]
    )
    return results


def log_read_and_run_data(reads_ns, args, results, logger):  # pragma nocover
    paired_str = "Paired" if reads_ns.paired else "Fragment"
    read_table=[
        ["Mean Read Length", reads_ns.read_length_mean],
        ["Read Length Standard Deviation", reads_ns.read_length_stddev ],
        ["Insert size", reads_ns.insert_mean],
        ["Insert size Standard Deviation", reads_ns.insert_stddev],
        ["Mean Long Read Length", reads_ns.long_read_length_mean],
        ["Mean Long Read Standard Deviation", reads_ns.long_read_length_stddev],
        ["Library type", paired_str],
        ["Platform", args.platform],
        ["Quality Encoding", reads_ns.encoding],
        ["Projected Coverage", str(reads_ns.coverage) + "x"],
        ["Projected Coverage (Downsampled)", str(reads_ns.downsampled_coverage) + "x"],
        ["Projected Long Read Coverage", str(reads_ns.long_read_coverage) + "x"]
        ]
    logger.info("LIBRARY DETAILS:\n" + tabulate.tabulate(read_table))

    run_table = [
        ["Assembly category",        reads_ns.lib_type],
        ["Selected assemblers",      " ".join(args.assemblers)],
        ["Selected assembly merger", args.merge_method],
        ["Selected scaffolder",      args.scaffolder],
        ["Selected finisher",        args.finisher],
        ["Selected variant caller",  args.varcaller],
        ["Trim QV",                  args.trim_qv],
        ["Trim length",              args.trim_length],
        ["Split Origin",             not args.skip_split_origin]
    ]
    logger.info("ASSEMBLER DETAILS:\n" + tabulate.tabulate(run_table))


def run_scaffolder(args, tools, config, reads_ns, results, run_id, logger=None):
    if tool.scaffolder['linkage_evidence'] == "paired-end":
        run_pe_scaffolder(args, tools, config, reads_ns, results, run_id, logger=logger)
    else:
        run_ref_scaffolder(args, tools, config, reads_ns, results, run_id, logger=logger)


def run_ref_scaffolder(args, tools, config, reads_ns, results, run_id, logger=None):
    """
    Runs specified scaffolder....

    Scaffolders which don't separate unscaffolded contigs can be wrapped
    in a script which created a "$scaffolder.contig_ids" file listing the
    IDs of contigs scaffolded. These will then be used following scaffolding
    to create our own file of unscaffolded contigs
    """
    logger.info(" Starting %s", tools.scaffolder['name'])
    exec_cmd = replace_placeholders(string=tools.scaffolder['command'],
                                    config=config, results=results,
                                    reads_ns=reads_ns, args=args)
    assembler_run = exec_cmd.split(" ")[0]
    exec_cmd = " ".join(exec_cmd.split(" ")[1:])
    scaffold_output = replace_placeholders(string=tools.scaffolder['scaffold_output'],
                                           config=config, results=results,
                                           reads_ns=reads_ns, args=args)
    unscaffolded_output = None
    if tools.scaffolder['unscaffolded_output']:
        unscaffolded_output = tools.scaffolder['unscaffolded_output']

    # no default args are currently implemeted
    default_args = tools.scaffolder['default_args']
    if args.scaffolder_args:
        exec_cmd = exec_cmd + args.scaffolder_args
    else:
        if default_args is not None:
            exec_cmd = exec_cmd + default_args

    run_dir = os.path.join(args.tmp_dir,
                           tools.scaffolder['name'] + "_" + str(run_id),
                           "")
    if True: # TODO make all scaffolder subscripts make a dir tools.scaffolder['create_dir']:
        os.makedirs(run_dir)

    # since we need to scaffold per-reference, which  doesn't work if your not using one..
    # list of paths of resulting scaffolds to be combined, if ref sequence has multiple records
    merge_these_scaffolds = []
    # this helps keep trak of resulting files
    ref_list = [] # [ref_id, ref_path]
    with open(args.reference, "r") as ref:
        for idx, rec in enumerate(SeqIO.parse(ref, "fasta")):
            ref_list.append([rec.id, os.path.join(run_dir, "reference_" + rec.id)])
            with open(ref_list[idx][1], "w") as outf:
                SeqIO.write( rec, outf, "fasta")
    # If the reference contains multiple contigs, we first neeed to align
    # our contigs to these
    # to identify which contigs to scaffold against which reference, since
    # some scaffolders targeted at bacteria don't handle multiple reference
    # sequences

    # moved ID parsing out, so we can keep naming the same
    if not len(ref_list) > 1:
        # if we don't have mulitple references, we just need to make the
        # reference and contigs available under consistent names
        reference = ref_list[0]
        # contigs = contigs_output
    else:
        # blast indexing doesn't produce a workable index from a symlink, so
        # need to copy the reference sequences
        reference_copy = os.path.join(run_dir, "reference_for_scaffolding.fasta")
        contigs_copy = os.path.join(run_dir, "contigs_to_be_scaffolded.fasta")
        shutil.copyfile(args.reference, reference_copy)
        shutil.copyfile(results.current_contigs, contigs_copy)
        blast_cmd = str("{0} -query {1} -subject {2}  -task blastn -out " +
                        "{3}clusters.blast 2>&1 > {3}blastn.log").format(
                            config.blastn, contigs_copy, args.reference, run_dir)
        blast_cmd2 = str("{0} -query {1} -subject {2} -outfmt 6 -out " +
                  "{3}clusters.tab 2>&1 > {3}blastn.log").format(
                      config.blastn, contigs_copy, args.reference, run_dir)
        blast_cmd3 = str("{0} -query {1} -subject {2} -outfmt 5 -evalue 0.01 -out " +
                  "{3}clusters.xml 2>&1 > {3}blastn.log").format(
                      config.blastn, contigs_copy, args.reference, run_dir)
        for cmd in [blast_cmd, blast_cmd2, blast_cmd3]:
            subprocess.run(cmd,
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
        # create dict of contigs per reference sequence, or unaligned
        # only need to worry about the top hit for each...
        ref_seqs = {"unaligned": []}
        with open(os.path.join(run_dir, "clusters.tab"), "r") as blast:
            # out fmt is queryID, subjectID, % match, length of match, and other stuff separated by tabs
            for line in blast:
                fields = line.strip().split("\t")
                query_id, subject_id = fields[0], fields[1]
                logger.debug("processing hit for %s", query_id)
                if subject_id not in ref_seqs.keys():
                    ref_seqs[subject_id] = [query_id]
                else:
                    ref_seqs[subject_id].append(query_id)
        logger.debug("generate a fasta file of contigs which align to each reference")
        for ref, cntgs in ref_seqs.items():
            out_name = os.path.join(run_dir, "reference_" + ref + "_contigs.fasta")
            with open(contigs_copy, "r") as contig_f:
                with open(out_name, "a") as out_f:
                    for rec in SeqIO.parse(contig_f, "fasta"):
                        if rec.id in cntgs:
                            SeqIO.write(rec, out_f,  "fasta")
    logger.debug("ref_list: %s", ref_list)
    for ref_id, ref_path in ref_list:
        scaff_contigs = os.path.join(run_dir, "reference_" + ref_id + "_contigs.fasta")
        if os.path.exists(scaff_contigs) and os.path.getsize(scaff_contigs) > 0:
            logger.info("Scaffolding vs. %s", ref_id)
            if unscaffolded_output:
                unscaffolded_output = replace_placeholders(
                    string=unscaffolded_output,config=config,
                    reads_ns=reads_ns, args=args)
            subscaff_dir = os.path.join(run_dir, "SIS_" + ref_id)
            os.makedirs(subscaff_dir)
            merge_these_scaffolds.append(os.path.join(subscaff_dir, "scaffolds.fasta"))
            # here we execute the "run" function from the appropriate runner script
            sub_mains[args.scaffolder](config=config, args=args, results=results,
                                       contigs=scaff_contigs,
                                       ref=ref_path,
                                       scaff_dir=subscaff_dir,
                                       logger=logger)
            resulting_scaffold = os.path.join(run_dir, "scaffolds.fasta")
    counter = 1
    with open(resulting_scaffold, "w") as outf:
        for scaf_path in merge_these_scaffolds:
            print(scaf_path)
            with open(scaf_path, "r") as scaf:
                for rec in SeqIO.parse(scaf, "fasta"):
                    rec.id = "scaffold_%05d" % counter
                    rec.desc = ""
                    SeqIO.write(rec, outf, "fasta")
                    counter = counter + 1

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


    print ("\nScaffolded assembly stats\n=========================\n\n")
    get_contig_stats(resulting_scaffold, 'scaffolds' );
    return resulting_scaffold

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
    # return linkage_evidence

# }

def run_pe_scaffolder(args, config, reads_ns, results, run_id, logger=None):
    """
    Runs specified scaffolder....

    Scaffolders which don't separate unscaffolded contigs can be wrapped
    in a script which created a "$scaffolder.contig_ids" file listing the
    IDs of contigs scaffolded. These will then be used following scaffolding
    to create our own file of unscaffolded contigs

    required params: $ (tmpdir)
                 $ (reference)
                 $ (scaffolder)
                 $ (scaffolder args)
                 $ (library insert size)
	           $ (library insert sd)
                 $ (run_id - appended to tmpdir to allow multiple runs)
                 $ (path to contigs to scaffold)
                 $ (mean read length)

    returns        : $ (linkage evidence type)
    """
    logger.info(" Starting $scaffolder");
    tool_name = args.scaffolder
    try:
        conf_scaffolder = [x for x in config.scaffolders if x['name'].lower() == args.scaffolder][0]
    except IndexError:
        raise ValueError("Scaffolder %s is not defined" % tool_name)
    exec_cmd = replace_placeholders(string=conf_scaffolder['command'],
                                    config=config, results=results,
                                    reads_ns=reads_ns, args=args)
    assembler_run = exec_cmd.split(" ")[0]
    exec_cmd = " ".join(exec_cmd.split(" ")[1:])
    scaffold_output = replace_placeholders(string=conf_scaffolder['scaffold_output'],
                                           config=config, results=results,
                                           reads_ns=reads_ns, args=args)

    unscaffolded_output = None
    if conf_scaffolder['unscaffolded_output']:
        unscaffolded_output = conf_scaffolder['unscaffolded_output']
    create = conf_scaffolder['create_dir']
    linkage_evidence = conf_scaffolder['linkage_evidence']
    # no default args are currently implemeted
    default_args = conf_scaffolder['default_args']
    if args.scaffolder_args:
        exec_cmd = exec_cmd + args.scaffolder_args
    else:
        if default_args is not None:
            exec_cmd = exec_cmd + default_args
    run_dir = os.path.join(args.tmp_dir, tool_name + "_" + str(run_id))
    os.makedirs(run_dir)
    # Treat reference-based scaffolder separately from paired-read scaffolders,
    # since we need to scaffold per-reference, which
    # doesn't work if your not using one...
    scaffolder_cmd_list = []  # commands to be subprocessed
    # list of paths of resulting scaffolds to be combined, if ref sequence has multiple records
    merge_these_scaffolds = []
    # this helps keep trak of resulting files
    ref_ids = []
    ref_paths = []
    with open(args.reference, "r") as ref:
        for idx, rec in enumerate(SeqIO.parse(ref, "fasta")):
            ref_ids.append(rec.id)
            ref_paths.append(os.path.join(run_dir, "reference_" + rec.id))
            with open(ref_paths[idx], "w") as outf:
                SeqIO.write( rec, outf, "fasta")
    # non reference-guided scaffolding, that is....
    # if no reference provided we won't have an estimate of insert size,
    # so need to get this by read alignment vs the assembly.
    if reads_ns.mean_insert is None:
        logger.info("Aligning reads to contigs to determine insert size")
        sorted_bam = align_reads(dirname="align_to_contigs", reads_ns=reads_ns,
                                 config=config, downsample=True, args=args,
                                 logger=logger)
        reads_ns.insert_mean, reads_ns.insert_stddev = get_insert_stats(
            bam=sorted_bam, config=config, args=args, reads_ns=reads_ns, logger=logger)

    run_scaffold_output = os.path.join(args.tmp_dir,
                                       args.scaffolder + "_no_reference",
                                       scaffold_output)

    # moved scaf_args out of if statement as well!
    exec_cmd = exec_cmd + " 2>&1 > " + os.path.join(
        args.tmpdir,
        args.scaffolder + "_" + str(run_id) + ".log")
    # scaffolder_cmd_list.append(exec_cmd)
    merge_these_scaffolds.append(os.path.join(run_scaffold_output, "scaffolds.fasta"))


    # when refactor, end get_scaffolder_cmds here
    # now we have are command(s), so run them

    # here we execute the "run" function from the appropriate runner script
    sub_mains[args.scaffolder](config=config, args=args, results=results,
                          scaff_dir=run_dir, logger=logger)
    # for cmd in scaffolder_cmd_list:
    #     pass
    # subprocess.run(exec_cmd,
    #                shell=sys.platform != "win32",
    #                stdout=subprocess.PIPE,
    #                stderr=subprocess.PIPE,
    #                check=True)
    resulting_scaffold = os.path.join(run_scaffold_output,
                                      "scaffolds.fasta")
    if len(merge_these_scaffolds) != 0:
        counter = 1
        with open(resulting_scaffold, "w") as outf:
            for scaf in merge_these_scaffolds:
                for rec in SeqIO.parse(scaf):
                    rec.id = "scaffold_%05d" % counter
                    SeqIO.write(rec, outf, fasta)
    return resulting_scaffold

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
    # return linkage_evidence

# }

def find_origin(args, config, results, logger):
    """
    Attempts to identify location of origin based upon contig overlapping
    base 1 of the reference sequence. This assumes  each reference sequence
    is a complete circular molecular i.e.a chromosome or a plasmid

    required params: $ (tmpdir)
                 $ (scaffolder)
                 $ (scaffolder_args)
                 $ (reference)
                 $ (insert_size)
                 $ (insert_stddev)
                 $ (read_length_mean)

    returns: $ (0)

    """
    ori_dir = os.path.join(args.tmp_dir, "origin")
    os.makedirs(ori_dir)
    logger.info("Attempting to identify origin...");
    cmds =  make_nucmer_delta_show_cmds(config=config, ref=args.reference,
                                       query=results.current_scaffolds,
                                       out_dir=ori_dir, prefix="ori", header=False)
    logger.debug("Running the following cmds:")
    for cmd in cmds:
        logger.debug(cmd)
        subprocess.run(cmd,
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    origin = None
    ori_coords = os.path.join(ori_dir, "ori.coords")
    with open(ori_coords, "r") as coords:
        for line in coords:
            line = line.strip()
            line = line.replace("|", "")
            fields = re.split("\s*", line)
            if not origin and fields[0] ==1:
                origin = "{0}:{1}".format(fields[8], fields[2])
                logger.info("Potential origin found at %s", origin )
    outIO = os.path.join(ori_dir, "split_ori_scaff.fasta")  # OutIO
    splitIO = os.path.join(ori_dir, "scaffolds_ori_split.fasta")  # SplitIO
    if origin is not None:
        ori_scaffold, pos = origin.split(":")
        with open(results.current_scaffolds, "r") as inscaff:
            for rec in SeqIO.parse(inscaff, "fasta"):
                if ori_scaffold == rec.id:
                    # write out both parts after splitting
                    with open(splitIO, "w") as splitscaff :
                        part_a = rec.seq[0: pos]
                        part_b = rec.seq[pos + 1: ]
                        ori_a = SeqRecord(id=rec.id + "_A",
                                          seq=Seq(part_a))
                        ori_b = SeqRecord(id=rec.id + "_B",
                                          seq=Seq(part_b))
                        SeqIO.write(ori_a, splitscaff, "fasta")
                        SeqIO.write(ori_b, splitscaff, "fasta")
                    #  now, rerun scaffolder on our split origins
                    results.old_scaffolds.append([results.current_scaffolds,
                                                  results.current_scaffolds_source])
                    results.current_scaffolds = split_ori_scaff
                    results.current_scaffolds_source = "split_ori_scaff"
                    new_scaffolds = run_scaffolder(
                        args=args, config=config, reads_ns=reads_ns, results=results,
                        run_id=2, logger=logger)
                    results.old_scaffolds.append([results.current_scaffolds,
                                                  results.current_scaffolds_source])
                    results.current_scaffolds = new_scaffolds
                    results.current_scaffolds_source = "find_origin_rescaffolded"

                    with open(new_scaffolds, "r") as infile, open(outIO, "a") as outfile:
                        for new_rec in SeqIO.parse(infile):
                            SeqIO.write(new_rec, outfile, "fasta")
                else:
                    # these dont need splitting
                    with open(outIO, "a") as outfile:
                        SeqIO.write(rec, outfile, "fasta")

    #renumber scaffolds to ensure they are unique...
    counter = 1
    renumbered_scaffs = os.path.join(
        os.path.dirname(splitIO), "scaffold.fasta")
    with open(splitIO, "r") as inf, open(renumbered_scaffs, "w") as outf:
        for rec in SeqIO.parse(inf):
            rec.id = "scaffold_%05d" % counter
            counter = counter + 1
            SeqIO.write(rec, outf, fasta)
    results.old_scaffolds.append([results.current_scaffolds, results.current_scaffolds_source])
    results.current_scaffolds = renumbered_scaffs
    results.current_scaffolds_source = "renumbering_after_origin_split"


def check_and_set_trim_length(reads_ns, args, logger):
    if reads_ns.read_length_mean is not None:
        if reads_ns.read_length_mean < 50 and reads_ns.read_length_mean < args.trim_length:
            logger.info("trim-length set to minimum of 25 due to mean read length: %d",
                        reads_ns.read_length_mean)
            return 25
    return args.trim_length

def use_already_assembled(args):
    for atr in ["already_assembled_dirs",
                "already_assembled_contigs",
                "already_assembled_scaffolds"]:
        if len(getattr(args, atr)) > 0:
            return True
    return False


def order_scaffolds():
    """
    Identifies origin based on homology with reference.
    Resulting scaffolds are then ordered and oriented relative to the reference...

    required params: $ (tmpdir)
                 $ (fasta reference)

    returns        : $ (0)
    """
    orient_dir = os.path.join(args.tmp_dir,  "orientating")
    os.makedirs(orient_dir)
    logger.info("Orienting scaffolds vs. reference...")
    cmds =  make_nucmer_delta_show_cmds(config=config, ref=args.reference,
                                       query=results.current_scaffolds,
                                       out_dir=orient_dir, prefix="ori2", header=False)
    logger.debug("Running the following cmds:")
    for cmd in cmds:
        logger.debug(cmd)
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)

    orient_coords = os.path.join(orient_dir, "ori2.coords")
    orient_count = { '+': 0,
                     '-': 0 }
    with open(orient_coords, "r") as inf:
        contig = '';
        for line in inf:
            line = line.strip()
            line = line.replace("|", "")
            fields = re.split("\s+", line)
            if contig != fields[8]:
                if contig != "":
                    store_orientation(orientations, contig, orient_count )
            contig = fields[8]
            start  = fields[2]
            end    = fields[3]
            length = ""
            if start < end:
                orient = '+'
                length = end - start;
            else:
                orient = '-'
                length = start - end
            if orient_count[orient]:
                orient_count[orient] = orient_count[orient] + length
            else:
                orient_count[orient] = length
    # record data for last contig...
    store_orientation(orientations, contig, orient_count )

    orig_length     = 0
    oriented_length = 0    #track how much sequence we align ok...
    unplaced_length = 0
    unplaced = []


    # Reorientate contigs and break origin, rewriting to a new file...
    new_scaffolds = os.path.join(orient_dir, "scaffolds.fasta")
    with open(results.current_scaffolds, "r") as inf, \
         open(new_scaffolds, "r") as outf:
        for i, rec in enumerate(SeqIO.parse(inf, "fasta")):
            orig_length = orig_length + len(rec.seq)
            rec.id = "scaffold_%5d"  % str(i + 1)
            if orientations[rec.id]:
                if orientations[rec.id] == '+':
                    SeqIO.write(rec, outf, "fasta")
                    oriented_length = oriented_length + len(rec.seq)
                else:
                    rc_rec = rec.reverse_complement()
                    SeqIO.write(rc_rec, outf, "fasta")
                    oriented_length = oriented_length + len(rec.seq)
            else:
                unplaced.append(rec)
                unplaced_length = unplaced_length + len(rec.seq)
        #  Now we have written the aligning contigs in order, write
        # the unaligned to the end of the scaffolds file
        # TODO check do these need names?
        for unpl_rec in unplaced:
            SeqIO.write(unpl_rec, outf, "fasta")
    # now we ned the length of the reference
    ref_length = 0
    with open(args.reference, "r") as ref:
        for rec in SeqIO.parse(ref, "fasta"):
            ref_length = ref_length + len(rec.seq)
    report = [
        ["", "Length (bp)"]
        ["Reference Sequence", ref_length],
        ["Assembly",           orig_length],
        ["Orientated contigs", oriented_length],
        ["Unaligned contigs",  unplaced_length]
    ]
    logger.info("ORIENTATING SCAFFOLDS RESULTS:\n%s", "\n".join(report))
    results.old_scaffolds.append([results.current_scaffolds, results.current_scaffolds_source])
    results.current_scaffolds = new_scaffolds
    results.current_scaffolds_source = "orient"


def store_orientation(orientations, contig, or_count):
    """
    save correct contig orientation in hash passed as first parameter

    since alignment fof scaffold can contain blocks in reverse orientation
    need to use the 'prevailing' orientation based on number of +/- blocks

    required parameters: $ (orientations hash)
                     $ (contig)
                     $ (hash of +/- base counts)

    returns           : none
    """
    if ((or_count['+'] and or_count['-']) and (or_count['+'] > or_count['-'])) \
       or (or_count['+'] and not or_count['-']):
        orientations[contig] = '+'
    else:
        orientations[contig] = '-'


def finish_assembly(args, config, reads_ns, results, logger):
    """
    Carried out assembly finishing using selected method
    required params: $ (tmpdir)
                 $ (finisher)
                 $ (insert size)
                 $ (insert stddev)
                 $ (base quality encoding)
                 $ (no. threads)

    returns        : $ (0)
    """
    logger.info("Finishing assembly with %s", args.finisher)
    tool = [x for x in config.finishers if x['name'].lower() == args.finisher.lower()]
    outdir = os.path.join(args.tmp_dir, args.finisher)
    if tool['create_dir']:
        os.makedirs(outdir)
    cmd = replace_placeholders(
        string=tool['command'],
        config=config,
        reads_ns=reads_ns, args=args)
    logger.debug("running the following command:\n %s". cmd)
    subprocess.run(cmd,
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    print ("Finished assembly statistics : \n============================\n")
    get_contig_stats( "$tmpdir/scaffolds.fasta", 'scaffolds' );


def main(args=None, logger=None):
    if args is None:
        args = get_args()
    check_files_present(args)
    results = make_empty_results_object()
    results.assemblers_list = match_assembler_args(args=args)
    results.organism = "{0}{1}{2}".format(args.genus, " " + args.species, " sp. " + args.strain)
    dt = arrow.now()
    if not args.outdir:
        if "unknown" in results.organism:
            args.outdir = "BugBuilder_" + dt.local.format('YYYY-MM-DD_HHmmss')
        else:
            args.outdir = "BugBuilder_" + results.organism.replace(" ", "_")
    args.outdir = os.path.abspath(os.path.expanduser(args.outdir))
    try:
        os.makedirs(args.outdir, exist_ok=False)
    except OSError:
        print("output directory exists! Exiting...")
        sys.exit(1)
    if logger is None:
        logger = set_up_logging(
            verbosity=args.verbosity,
            outfile=os.path.join(args.outdir, "BugBuilder.log"),
            name=__name__)
    for k, v in sorted(vars(args).items()):
        logger.debug("%s: %s", k, str(v))

    # are we dealing with a paired library?  set that within assess_reads
    logger.info("Welcome to BugBuilder")
    logger.info("Preparing to build your bug...")
    config_path = get_config_path()
    try:
        config = return_config(config_path, force=False, logger=logger)
    except:
        logger.error("Error reading config file! try running BugBuilder " +
                     "--configure to fix, or else delete the config file, " +
                     "reconfigure, and try again")
        sys.exit()
    logger.info("Using configuration from %s", config_path)
    logger.debug(config)
    # lets not be hasty
    check_args(args, config)
    setup_tmp_dir(args, output_root=args.outdir, logger=logger)
    # if args.reference is None:
    #     reference = None
    # else:
    #     reference = os.path.basename(args.reference)
    #  need to know a bit about the provided reads so we can act appropriately

    # copy paths to raw reads in case we need them for mascura
    args.untrimmed_fastq1 = args.fastq1
    args.untrimmed_fastq2 = args.fastq2

    logger.info("Assessing reads and library type")
    reads_ns  = assess_reads(args=args, config=config,
                             platform=args.platform, logger=logger)
    logger.debug(reads_ns)

    logger.info("Determining if a reference is needed")
    check_ref_needed(args=args, lib_type=reads_ns.lib_type)

    logger.info("Preparing config and tools")
    tools = select_tools(args, config=config, reads_ns=reads_ns, logger=logger)

    logger.debug("Determining if fastqc will be run")
    if not args.skip_fastqc and config.fastqc is not None:
        logger.info("Running fastqc on reads")
        run_fastqc(reads_ns, args, logger=logger)
    else:
        logger.info("Skipping fastqc")

    logger.debug("Ensuring appropriate trim length")
    args.trim_length = check_and_set_trim_length(reads_ns, args, logger)

    logger.debug("Determining whether to perfrom trimming")
    if reads_ns.lib_type != "long" and not args.skip_trim:
        logger.info("Trimming reads based on quality")
        quality_trim_reads(args, config, reads_ns, logger)

    logger.debug("Determining whether to downsample")
    if (reads_ns.coverage is not None and reads_ns.coverage > 100) and \
       (assembler_needs_downsampling(tools) or args.downsample != 0):
        logger.info("Downsampling reads to %dx coverage", args.downsample)
        reads_ns.downsampled_coverage = downsample_reads(
            args=args, reads_ns=reads_ns,
            config=config, new_cov=args.downsample)

    if args.fastq2 and args.reference:
        logger.info("Aligning reads to reference for determinging insert size")
        sorted_bam = align_reads(dirname="align", reads_ns=reads_ns,
                                 config=config, downsample=True, args=args, logger=logger)
        reads_ns.insert_mean, reads_ns.insert_stddev = get_insert_stats(
            bam=sorted_bam, config=config, args=args, reads_ns=reads_ns, logger=logger)

    # print out the run data we have so far
    log_read_and_run_data(reads_ns, args, results, logger)

    # Righto, now the fun starts. First up, we assemble with 1 or 2 assemblers
    # (or check and read in results from a prior assembly)

    if use_already_assembled(args):
        logger.info("Using provided assembly(s) results")
        # if provided contigs and scaffolds
        if len(args.already_assembled_dirs) == 0:
            results.assemblers_results_dict_list = check_already_assembled_args(
                args, config, logger)
        # if provided dir(s)
        else:
            results.assemblers_results_dict_list = check_already_assembled_dirs(
                args, config, logger)
    else:
        logger.info("Running Assembler(s)")
        for assembler, assembler_args in results.assemblers_list:
            logger.info("Assembling with %s", assembler)
            contigs_path, scaffolds_path  = run_assembler(
                assembler=assembler, assembler_args=assembler_args,
                args=args, reads_ns=reads_ns, config=config, logger=logger)
            results.assemblers_results_dict_list.append(
                {"name": assembler,
                 "contigs":contigs_path,
                 "scaffolds": scaffolds_path})

    # if we have more than one assembler, run merge and point results to merged
    if len(results.assemblers_results_dict_list) > 1:
        #run merge
        merged_contigs_path = merge_assemblies(args=args, config=config,
                                               reads_ns=reads_ns, logger=logger)
        results.current_contigs = merged_contigs_path
        results.current_scaffolds = None

    else:
        results.current_contigs = \
            results.assemblers_results_dict_list[0]['contigs']
        results.current_scaffolds = \
            results.assemblers_results_dict_list[0]['scaffolds']

    print(results.__dict__)
    print(tools.__dict__)
    print(args.__dict__)

    #TODO why do we check scaffolder here?
    if args.scaffolder and args.reference:
        results.ID_OK = check_id(args, contigs=results.current_contigs,
                                 config=config,
                                 logger=logger)
    if tools.scaffolder is not None:
        if \
           (results.ID_OK and tools.scaffolder['linkage_evidence'] == "align-genus") or \
             not (results.ID_OK and tools.scaffolder['linkage_evidence'] == "paired-ends"):
            results.old_scaffolds.append(["either nowhere or straight from assembler",
                                         results.current_scaffolds])
            results.current_scaffolds = run_scaffolder(
                args=args, config=config, reads_ns=reads_ns, results=results,
                run_id=1, logger=logger)

    if results.current_scaffolds is not None:
        # we may have been given a reference for a long read assembly but no
        # scaffolder is used for these by default
        args.scaffolder = "mauve" if args.scaffolder is None else args.scaffolder
        if not args.skip_split_origin and args.reference is not None:
            find_origin(args,config, results,  logger )
        if results.ID_OK:
            order_scaffolds(args, logger)
        if args.finisher and results.ID_OK:
            finish_assembly(args, config, reads_ns, results, logger=logger)

    if results.current_scaffolds is not None:
        pass
         # gaps = build_agp( $tmpdir, $organism, $mode, $scaffold_type );

#     # sequence stable from this point, only annotations altered
#     for my $i (qw(1 2)) {
# 	print "pm: $i \n";
#         $pm->start and next();
#         if ( $i == 1 ) {

#             ##amosvalidate fails if we don't have mate-pairs
#             if ( -e "$tmpdir/read2.fastq" && ( $mode eq 'draft' ) ) {
#                 my $seq_file;
#                 ( -e "$tmpdir/scaffolds.fasta" ) ? ( $seq_file = 'scaffolds.fasta' ) : ( $seq_file = 'contigs.fasta' );

#                 align_reads( $tmpdir, $seq_file, $read_length_mean );
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

#     run_varcaller( $tmpdir, $varcall, $threads, $read_length_mean ) if ($varcall);
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






















def build_agp():
    """
    Creates an AGP file from the scaffolds, while generating new
    contig/scaffold outputs meeting ENA requirements (no consecutive runs
    of >=10 N, minimum contig size of 200 bp) if running in 'submission'
    mode, otherwise leaves short contigs and gaps<100bp intact. The
scaffold_type argument is used to determine the linkage evidence type
    for scaffold gaps

    required parameters: $ (tmpdir)
		     : $ (organism description)
                   : $ (mode - submission or draft)
                   : $ (scaffold_type: align or mate_pair)

    returns            : $ (0)
    """


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
