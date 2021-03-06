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


import argparse
import os
import string
import hashlib
import re
import yaml
import sys
import time
import shutil
import logging
import statistics
import subprocess
import pkg_resources
import tabulate
import multiprocessing



from BugBuilder import __version__
from Bio import SeqIO # , SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from argparse import Namespace

from .program_info import all_programs_dict
from .shared_methods import make_nucmer_delta_show_cmds, run_nucmer_cmds
from .raw_config import __config_data__
# Assemblers
from .run_spades import run as run_spades
from .run_skesa import run as run_skesa
from .run_abyss import run as run_abyss
# Scaffolders
from .run_sis import run as run_sis
from .run_ragout import run as run_ragout
# finishers
from .run_pilon import run as run_pilon
from .run_fgap import run as run_fgap
from .run_sealer import run as run_sealer
# varcallers
from .run_pilon_var import run as run_pilon_var


sub_mains = {
    # "canu": run_canu,
    "spades": run_spades,
    "abyss": run_abyss,
    "sis": run_sis,
    "ragout": run_ragout,
    "pilon": run_pilon,
    "fgap": run_fgap,
    "abyss-sealer": run_sealer,
    "pilon_var": run_pilon_var,
}

##################################################################
# a few constants

# stages [stage_flag, name, code, prereq_stages
all_stages = [
    ["q", "QC reads", "QC", None],
    ["t", "trim reads", "TRIM", None],
    ["d", "downsample reads", "DOWNSAMPLE", None],
    ["a", "run assembler(s)", "ASSEMBLE", None],
    ["b", "break contig at the origin ", "BREAK", "sa"],
    ["s", "scaffold the contigs", "SCAFFOLD", "a"],
    ["f", "run genome polishing", "FINISH", "as"],
    ["v", "run variant caller", "VARCALL", "a"],
    ["g", "gene call with prokka", "GENECALL", "a"]
]

##################################################################

def parse_available(thing, path=None):
    if path is None: # path var is to allow easier testing
        config_path = get_config_path()
        try:
            config = return_config(config_path, force=False, hardfail=False, logger=None)
        except Exception as e:
            return []
    else:
        config = return_config(path, force=False, hardfail=False, logger=None)

    # get all the potential names from the config
    thing_list = [x['name'].lower() for x in getattr(config, thing)]
    # return all the names that have an executable available
    # sometimes, if configuring fails, this can result in an attribute arrow
    try:
        return [x for x in thing_list if getattr(config, x.replace("-", "_")) is not None]
    except AttributeError:
        print("It looks like there is a damaged config file; Fixing!")
        with open(get_config_path(), "w") as conf:
            conf.write("STATUS: INCOMPLETE")
        configure(get_config_path())
        return(parse_available(thing=thing, path=get_config_path()))


def check_version_from_cmd(
        exe,
        cmd, line,
        pattern=r"^__version__ = '(?P<version>[^']+)'$",
        where='stderr',
        min_version="0.0.0", logger=None,
        coerce_two_digit=False):
    """ returns version string from an system call
    Hacky, but better than nothing.
    line arg is 1-indexed
    .strip() is called on match to remove whitspaces
    20170920 changed to remove shutil.which call.
    That should be done outside of this funciton
    """
    assert logger is not None, "must use logging"
    from distutils.version import StrictVersion
    result = subprocess.run("{0} {1}".format(exe, cmd),
                             # is this a securiy risk?
                            shell=sys.platform != "win32",
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            check=False)
    logger.debug(result)
    try:
        if where == 'stderr':
            printout = result.stderr.decode("utf-8").split("\n")
        elif where == 'stdout':
            printout = result.stdout.decode("utf-8").split("\n")
        else:
            raise ValueError("where option can only be 'stderr' or 'stdout'")
    except Exception as e:
        raise e
    if logger:
        logger.debug(printout)
    this_version = None
    try:
        m = re.search(pattern, printout[line - 1])
    except IndexError as e:
        raise e
    if m:
        this_version = m.group('version').strip()
    if logger:
        logger.debug("this_version: %s", this_version)
    if coerce_two_digit:
        this_version = "0.{0}".format(this_version)
        if logger:
            logger.debug("coerced this_version: %s", this_version)
    if this_version is None:
        raise ValueError("No verison was captured with pattern" +
                         "{0}".format(pattern))
    try:
        if StrictVersion(this_version) < StrictVersion(min_version):
            raise ValueError("{0} version {1} must be greater than {2}".format(
                cmd, this_version, min_version))
    except Exception as e:
        raise e
    return(this_version)


def configure(config_path, hardfail=True):
    with open(config_path, 'w') as outfile:  # write header and config params
        for line in __config_data__:
            outfile.write(line)
    # imported from program_info.py
    output_dict = {}
    for category, programs in all_programs_dict.items():
        if len(programs) == 0: # pyyaml will write out curly brackets otherwise
            continue
        for prog, name in programs:
            # trim off extension, replace dashes with underscores
            clean_name = os.path.splitext(name)[0].lower().replace("-", "_")
            if shutil.which(prog):
                output_dict[clean_name] = shutil.which(prog)
            else:
                output_dict[clean_name] = None
                if "mandatory" in category:
                    if hardfail:
                        raise OSError(
                            "%s is a mandatory program; please install" % name)
    # get conda env too
    conda_env = subprocess.run("conda list",
                               shell=sys.platform != "win32",
                               stdout=subprocess.PIPE)
    with open(config_path, 'a') as outfile:  # write paths to exes in config
        yaml.dump(output_dict, outfile, default_flow_style=False)
        for line in conda_env.stdout.decode("utf-8").split("\n"):
            outfile.write("## " + line + "\n")


def return_config(config_path, force=False, hardfail=True, logger=None):
    """returns config namespace object, configuring if needed
    """
    this_print = logger.debug if logger else print

    try: # if the config file exists, parse it
        config = parse_config(config_path)
    except Exception as e:
        if e == yaml.YAMLError:
            this_print("Error parsing YAML file!")
        elif e == FileNotFoundError:
            this_print("Config file not found at %s! creating one..." % config_path)
        else:
            this_print(e)
        force = True
    # if config file is missing or broken, recreate the config file
    if force or config.STATUS != "COMPLETE":
        this_print("(Re-)Configuring BugBuilder config file: %s" % config_path)
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


class SmartFormatter(argparse.HelpFormatter):
    """ poached from https://stackoverflow.com/questions/3853722
    """
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


def get_args():  # pragma: no cover
    """
    """
    parser = argparse.ArgumentParser(
        prog="BugBuilder",
        description="Bugbuilder %s: " % __version__ +
        "Automated pipeline for assembly of draft quality " +
        "bacterial genomes with reference guided scaffolding and annotation." +
        "Please see accompanying userguide for full documentation" +
        "" +
        "Please report any bugs/issues via github:" +
        "https://github.com/jamesabbott/BugBuilder/issues/new." +
        "All bug reports should include the contents of your config file" +
        " (find when running `BugBuilder --configure`)" +
        "which reports on the installed software packages and versions",
        formatter_class=SmartFormatter,
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

    opt = parser.add_argument_group('optional arguments')
    opt.add_argument("--fastq1", dest='fastq1', action="store",
                     help="Path to first read of paired library, or fragment library",
                     type=str, required="--long-fastq" not in sys.argv)
    opt.add_argument("--fastq2", dest='fastq2', action="store",
                     help="Path to second read of paired library",
                     type=str)
    opt.add_argument("--untrimmed_fastq1", dest='untrimmed_fastq1', action="store",
                     # help="used to hold path of raw F reads for use " +
                     # "with mascura; dont set from command line",
                     type=str, help=argparse.SUPPRESS)
    opt.add_argument("--untrimmed_fastq2", dest='untrimmed_fastq2', action="store",
                     # help="used to hold path of raw R reads for use " +
                     # "with mascura; dont use from command line",
                     type=str, help=argparse.SUPPRESS)
    opt.add_argument("--de_fere_contigs", dest='de_fere_contigs',
                     action="store",
                     help="contigs to be used in de fere novo assembly!",
                     type=str)
    opt.add_argument("--long_fastq", dest='long_fastq', action="store",
                     help="Path to fastq file from long-read sequencer",
                     type=str, required="--fastq1" not in sys.argv)
    opt.add_argument("--references", dest='references',
                     action="store",
                     nargs="*", type=str,
                     default=[],
                     help="Path(s) to fasta formatted reference sequence(s)")
    opt.add_argument("--prefix", dest='prefix',
                     action="store",
                     help="Prefix to use for output file naming",
                     type=str)
    opt.add_argument("--assembler-args", dest='assembler_args',
                     action="store",
                     help="Any additional arguments to pass to the " +
                     "assembler. Default values are set in the " +
                     "'default_args' attribute of the configuration " +
                     "file. If running multiple assemblers, " +
                     "assembler_args should be specified twice, once " +
                     "for each assemler, in the same order than the " +
                     "assemblers are specified. enclose with quotes and " +
                     "start each with a space, like ' -sample'",
                     nargs="*", default=[],
                     type=str)
    opt.add_argument("--scaffolder", dest='scaffolder', action="store",
                     help="scaffolder to use",
                     choices=parse_available("scaffolders"),
                     type=str)
    opt.add_argument("--scaffolder-args", dest='scaffolder_args',
                     action="store",
                     help="Any additional arguments to pass to the " +
                     "scaffolder. Overides the setting of  'default_args' " +
                     "setting in the scaffolder configuration.  " +
                     "Use single quotes",
                     type=str)
    opt.add_argument("--merger", dest='merger', action="store",
                     help="Assembly merging tool to use",
                     choices=parse_available("merge_tools"),
                     type=str)
    opt.add_argument("--finisher", dest='finisher', action="store",
                     help="Method for assembly finishing/gap closure",
                     choices=parse_available("finishers"),
                     type=str)
    opt.add_argument("--varcaller", dest='varcaller', action="store",
                     help="Method for variant calling",
                     choices=parse_available("varcallers"),
                     type=str)
    opt.add_argument("--insert-size", dest='insert_size', action="store",
                     help="Size (bp) of insert in paired-read library.  " +
                     "This will be determined empircally if a reference  " +
                     "genome sequence is provided, so only needs specifying  " +
                     "when assembling paired-read sequences for which no  " +
                     "reference genome is available.",
                     type=int)
    opt.add_argument("--insert-stddev", dest='insert_stddev', action="store",
                     help="standard deviation of insert in paired-read  " +
                     "library. This will be determined empircally if a  " +
                     "reference genome sequence is provided, so only needs  " +
                     "specifying when assembling paired-read sequences for  " +
                     "which no reference genome is available ",
                     type=int)
    opt.add_argument("--genome-size", dest='genome_size', action="store",
                     help="Approximate genome size. Required for  " +
                     "PacBio/MinION assemblies when not using a reference ",
                     type=int, default=0) # 0 is better than None for addition :)
    opt.add_argument("--downsample", dest='downsample', action="store",
                     help="Downsample depth; default is 100x ",
                     type=int, default=100)
    opt.add_argument("--species", dest='species', action="store",
                     help="Specific name of species if known (ie pyogenes)." +
                     " Included in resutling EMBL file, and passed to  " +
                     "Prokka during annotation stage.",
                     default="unknown_species",
                     type=str)
    opt.add_argument("--genus", dest='genus', action="store",
                     help="Genus of organism sequenced (ie. Streptococcus). " +
                     " Included in resutling EMBL file, and passed to  " +
                     "Prokka during annotation stage.",
                     default="unknown_genus",
                     type=str)
    opt.add_argument("--strain", dest='strain', action="store",
                     help="Name of strain used for inclusion in annotation " +
                     " results", default="unknown_strain",
                     type=str)
    opt.add_argument("--locustag", dest='locustag', action="store",
                     help="Locustag argument to pass to Prokka. Use to  " +
                     "customise locus_tag in generated EMBL records.",
                     type=str)
    opt.add_argument("--centre", dest='centre', action="store",
                     help="Sequence centre argument to pass to Prokka. " +
                     "Used to customise locus_tag in generated EMBL records.",
                     type=str)
    opt.add_argument("--mode", dest='mode', action="store",
                     help="Mode to run in",
                     choices=["submission", "draft"], default="submission",
                     type=str)
    opt.add_argument("--trim-qv", dest='trim_qv', action="store",
                     help="Quality threshold for trimming reads.;  " +
                     "default: %(default)s  ", default=20,
                     type=int)
    opt.add_argument("--trim-length", dest='trim_length', action="store",
                     help="Min. length of read to retain following  " +
                     "trimming. Default: 50 (25 for reads <50bp)", default=50,
                     type=int)
    opt.add_argument("--keepall", dest='keepall',
                     action="store_true",
                     help="keep all intermediate files ", default=False)
    opt.add_argument("--threads", dest='threads', action="store",
                     help="Number of threads to use for multi-threaded  " +
                     "applications",
                     type=int, default=4)
    opt.add_argument("--out-dir", dest='out_dir', action="store",
                     help="dir for results",
                     type=str)
    opt.add_argument("--tmp-dir", dest='tmp_dir', action="store",
                     help="Path for the temporary working directory." +
                     "If not provided, it will be created under  " +
                     "the out-dir path.",
                     type=str)
    opt.add_argument("--already_assembled_dirs", dest='already_assembled_dirs',
                     action="store",
                     help="dir(s) with existing assembler results; must" +
                     " specify assembler(s), and lists much match " +
                     "length", default=[],
                     nargs="*", type=str)
    # These can be used instead of the --already-assembled_dirs arg
    opt.add_argument("--contigs",
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
    opt.add_argument("--scaffolds",
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
    opt.add_argument("--configure", action=JustConfigure, type=None,
                     help="Reconfigure and exit ", nargs=0)
    opt.add_argument("--memory", dest='memory', action="store",
                     help="how much memory to use",
                     type=int, default=4)
    opt.add_argument("--stages", dest='stages', action="store",
                     help="R| which stages of BugBuilder will run:\n" +
                     "  q = QC reads\n" +
                     "  t = trim reads\n" +
                     "  d = downsample reads\n" +
                     "  a = run assembler (s)\n" +
                     "  b = break contig at the origin\n" +
                     "  s = scaffold the contigs\n" +
                     "  f = run polishing finisher on the assembly\n"+
                     "  v = call variants\n"
                     "  g = gene-call with prokka\n" +
                     "default: %(default)s",
                     type=str, default="qtdabsfvg")
    opt.add_argument("-v", "--verbosity", dest='verbosity',
                     action="store",
                     default=2, type=int, choices=[1, 2, 3, 4, 5],
                     help="writes debug to file in output dir;\n" +
                     "this sets verbosity level sent to stderr. \n" +
                     "  1 = debug()\n  2 = info()\n  3 = warning()\n" +
                     "  4 = error() \n  5 = critical(); " +
                     "default: %(default)s")
    opt.add_argument("-h", "--help",
                     action="help", default=argparse.SUPPRESS,
                     help="Displays this help message")
    opt.add_argument('--version', action='version',
                     version='BugBuilder {}'.format(__version__))
    # args = parser.parse_args(sys.argv[2:])
    args = parser.parse_args()
    return args


def set_up_logging(verbosity, outfile, name): # pragma: no cover
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
    try:
        logfile_handler = logging.FileHandler(outfile, "w")
        logfile_handler.setLevel(logging.DEBUG)
        logfile_handler_formatter = \
            logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        logfile_handler.setFormatter(logfile_handler_formatter)
        logger.addHandler(logfile_handler)
    except:
        logger.error(e, exc_info=True)
        logger.error("Could not open {0} for logging".format(outfile))
        sys.exit(1)
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to given verbosity
    console_err = logging.StreamHandler(sys.stderr)
    console_err.setLevel(level=(verbosity * 10))
    console_err_format = logging.Formatter(
        str("%(asctime)s \u001b[3%(levelname)s\033[1;0m  %(message)s"), "%H:%M:%S")
    console_err.setFormatter(console_err_format)
    # set some pretty colors, shorten names of loggers to keep lines aligned
    #  the "4m " etc is a bit annoying, but I'd rather have a colored console
    #  than a spotless log file
    logging.addLevelName(logging.DEBUG,    "4m ..")
    logging.addLevelName(logging.INFO,     "2m --")
    logging.addLevelName(logging.WARNING,  "3m !!")
    logging.addLevelName(logging.ERROR,    "1m xx")
    logging.addLevelName(logging.CRITICAL, "1m XX")
    logger.addHandler(console_err)
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


def parse_stages(args, reads_ns, logger, all_stages):
    """ pare the "stages" arg to determine the steps to be performed
    We want to be able to run all the stages individualy if possible
    """
    allowed = [x[0] for x in all_stages]
    illegal_in_arg = [letter for letter in args.stages if letter not in allowed]
    if len(illegal_in_arg) > 0:
        raise ValueError("%s is not a valid stage; see --help" %
                         " and ".join(illegal_in_arg))
    logger.debug("Stages to be run:")
    for flag, message, attribute, prereqs in all_stages:
        if flag in args.stages:
            setattr(reads_ns, attribute, True)
            logger.debug(message)
        else:
            setattr(reads_ns, attribute, False)
    # ok, so our stages are decided, but lets see if they are valid


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
    file_list = args.references[:]  # beware this needs to be a deep copy
    file_list.extend([args.fastq1, args.fastq2])
    for f in file_list:
        if f is not None and not os.path.exists(f):
            raise ValueError("file %s does not exist!" % f)
    if args.fastq2 is not None and args.fastq1 == args.fastq2:
        raise ValueError("fastq1 and fastq2 are the same!")
    if args.mode == "draft" and len(args.references) == 0:
        raise ValueError("Draft mode requiresa reference!")


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
    logger.debug("Creating working directory: %s" % args.tmp_dir)
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
    if len(args.references) != 0:
        new_refs = []
        for reference in args.references:
            new_reference = os.path.join(args.tmp_dir, os.path.basename(reference))
            with open(reference, "r") as inf:
                for rec in SeqIO.parse(inf, "fasta"):
                    with open(new_reference, "a") as outf:
                        rec.desc = None
                        SeqIO.write(rec, outf, 'fasta')
            new_refs.append(new_reference)
        args.references = new_refs


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
        logger.debug("mean read length is {0} in {1}".format(
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
    logger.debug("Quality score Min: %d;  Max:  %d", min_val, max_val)
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
    logger.debug("Read encoding detection summary:\n" + tabulate.tabulate(
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
    logger.debug("Checking reads...");
    counts = []
    means = []
    stddevs = []
    long_mean = 0
    long_stddev = 0
    tot_length = 0
    long_tot_length = 0
    short_libs = 0 # to be used to calc coverage
    # get the lengths of all the short reads
    type_list = ["short", "short", "long"]
    for i, read in enumerate([args.fastq1, args.fastq2, args.long_fastq]):
        if  read is None or not os.path.exists(read):
            continue
        if i < 2:  # ie, not the long_fastq
            short_libs = short_libs + 1
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
        if len(args.references) != 0:
            logger.debug("Infering genome size from reference(s)")
            size_list = []  # [ref, size]
            for ref in args.references:
                this_size = 0
                with open(ref, "r") as inf:
                    for rec in SeqIO.parse(inf, "fasta"):
                        this_size = this_size + len(rec.seq)
                size_list.append([ref, this_size])
            logger.debug("Refence Genome sizes:")
            logger.debug(size_list)
            args.genome_size = int(statistics.mean([x[1] for x in size_list]))
    logger.info("Genome Size: %d", args.genome_size)
    logger.debug("Length of all reads: %d", tot_length)
    coverage, long_coverage = None, None
    try:
        # (total length of libraries/genome_size) / number of libraries
        # coverage = round(float((tot_length / args.genome_size) / short_libs), 3)
        # think this doesnt matter how many libries are involved
        coverage = round(float(tot_length / args.genome_size), 3)
        long_coverage = round(float(long_tot_length / args.genome_size), 3)
    except:
        logger.error("Cannot calculate coverage without a reference!")
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


def check_ref_required(args, lib_type):
    """ ensure either a reference or genome size is provided for long read assembly
    """
    if lib_type == "long" :
        if args.genome_size == 0 and len(args.references) == 0:
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
    logger.debug("Requested assemblers: %s", " ".join(args.assemblers))
    if len(args.assemblers) > 2:
        raise ValueError("A maximum of 2 assemblers can be requested")
    # this should have been checked on config, checking here makes testing hard
    # # sanity check assemblers where these have been manually requested...
    # for assembler in args.assemblers:
    #     if shutil.which(assembler) is None and \
    #        shutil.which(assembler + ".py") is None:
    #         raise ValueError("%s is nt in PATH! exiting" % assembler)
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
                                 (assembler, conf_assembler['min_length']))
            if conf_assembler['max_length'] and reads_ns.read_length_mean > conf_assembler['max_length']:
                raise ValueError("%s does not support reads greather than %d" % \
                                 (assembler, conf_assembler['max_length']))
            if conf_assembler['insert_size_required'] and not reads_ns.insert_mean and len(args.references) == 0:
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
    if args.scaffolder not in [k['name'].lower() for k in config.scaffolders]:
        raise ValueError("%s not an available scaffolder!" %args.scaffolder)
    conf_scaffolder = [x for x in config.scaffolders if \
                       x['name'].lower() == args.scaffolder][0]
    if len(conf_scaffolder) == 0:
        raise ValueError("%s no a configured scaffolder" % args.scaffolder)
    if conf_scaffolder['linkage_evidence'] == "paired-ends" and not paired:
        raise ValueError(str(
            "%s requires paired reads, but you only specified " +
            "one fastq file.") % args.scaffolder)
    elif "align" in conf_scaffolder['linkage_evidence'] and len(args.references) == 0:
        raise ValueError(str("%s requires a reference for alignment, " +
                             "but none is specified.") % args.scaffolder)
    else:
        return (conf_scaffolder, conf_scaffolder['linkage_evidence'])

def get_merger_tool(args, config, paired):
    """must be used before scaffolder getter
    """
    if args.merger is None:
        return None
    if args.merger not in [x['name'] for x in config.merge_tools]:
        raise ValueError("%s not an available merger!" % args.merger)
    linkage_evidence = None
    for merger_method in config.merge_tools:
        if merger_method['name'] == args.merger:
            if not merger_method['allow_scaffolding']:
                raise ValueError(
                    "%s prohibits scaffolding! Please remove the 's' from " +
                    "--stages, change scaffolder, or change merge tool!")
                args.scaffolder = None
            return merger_method

def get_finisher(args, config, paired):
    if args.finisher is None:
        return None
    # print(config.finishers)
    if args.finisher.lower() not in [x['name'].lower() for x in config.finishers]:
        raise ValueError("%s not an available finisher!" %args.finisher)
    for conf_finisher in config.finishers:
        if conf_finisher['name'] == args.finisher:
            if len(args.references) == 0 and conf_finisher['ref_required']:
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
    if args.varcaller.lower() not in [x['name'].lower() for x in config.varcallers]:
        raise ValueError("%s not an available varcaller!" %args.varcaller)
    for conf_varcallers in config.varcallers:
        if conf_varcallers['name'] == args.varcaller:
            if len(args.references) == 0 and conf_varcallers['ref_required']:
                raise ValueError("%s requires reference"  % args.varcaller)
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
    merger = get_merger_tool(args, config, paired=reads_ns.paired)
    scaffolder_tool, linkage_evidence = get_scaffolder_and_linkage(
        args=args, config=config, paired=reads_ns.paired, logger=logger)
    logger.debug("getting merger")
    finisher_tool = get_finisher(args, config, paired=reads_ns.paired)
    varcall_tool = get_varcaller(args, config, paired=reads_ns.paired)

    logger.debug("linkage: %s" % linkage_evidence)
    return Namespace(assemblers=assembler_tools,
                     scaffolder=scaffolder_tool,
                     scaffold_type=linkage_evidence,
                     merger=merger,
                     finisher=finisher_tool,
                     varcaller=varcall_tool)


def assembler_needs_downsampling(tools):
    for assembler in tools.assemblers:
        if assembler['downsample_reads']:
            return True
    return False


def make_fastqc_cmd(exe, args, outdir):
    cmd = \
        "{0} -t {5} --extract -o {1}{2}{3}{4} > {1}fastqc.log 2>&1".format(
            exe,
            outdir,
            " " + args.fastq1 if args.fastq1 is not None else "",
            " " + args.fastq2 if args.fastq2 is not None else "",
            " " + args.long_fastq if args.long_fastq is not None else "",
            args.threads)
    return cmd


def run_fastq_cmd(cmd, logger=None):
    logger.debug(cmd)
    fastqc_res = subprocess.run(cmd,
                                shell=sys.platform != "win32",
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                check=False)
    if fastqc_res.returncode != 0:
        logger.warning("Error running fastqc with command %s", fastqc_cmd)
        sys.exit(1)

def run_fastqc(reads_ns, args, config, logger=None):
    """
    Carries out QC analysis using fastqc...

    required params: $ (tmpdir)
                 $ (type)

    returns: $ (0)
    """
    logger.info("Running FastQC...")
    fastqc_dir = os.path.join(args.tmp_dir, "fastqc", "")
    os.makedirs(fastqc_dir)
    fastqc_cmd = make_fastqc_cmd(exe=config.fastqc, args=args, outdir=fastqc_dir)
    run_fastq_cmd(fastqc_cmd, logger=logger)
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
                        os.path.join(args.outdir, name + ".html"))
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


def report_trim_quality(trim_log, logger):
    """ look through the sickle log and warn if we have lost over 10% of reads
    """
    kept = 0
    disc = 0
    report = []
    with open(trim_log, "r") as inf:
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
    lost_percent = float(disc / kept) * 100
    if  lost_percent > 10:
        logger.warning(">10%% of reads (%d%%) discarded during read trimming. "+
                       "Please examine the FastQC outputs...", lost_percent)
    else:
        logger.debug("Discarded %d%% of the reads during trimming",
                     lost_percent)

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
    # reasign the reads with the trimmed treads
    if reads_ns.paired:
        args.fastq1 = os.path.join(trim_dir, "read1.fastq")
        args.fastq2 = os.path.join(trim_dir, "read2.fastq")
    else:
        args.fastq1 = os.path.join(trim_dir, "read1.fastq")

    report_trim_quality(trim_log=os.path.join(trim_dir, "sickle.log"), logger=logger)


def make_seqtk_ds_cmd(args, reads_ns, new_coverage, outdir, config, logger):
    """
    """
    assert isinstance(new_coverage, int), "new_coverage must be an int"
    frac = float(new_coverage / reads_ns.coverage)
    cmd_list = []
    logger.info("downsampleing from %f X  to %f X (%f)",
                reads_ns.coverage, new_coverage , frac )
    for reads in [args.fastq1, args.fastq2]:
        print(reads)
        if reads is not None and os.path.exists(reads):
            cmd_list.append("{0} sample -s 100 {1} {2} > {3}".format(
                config.seqtk,  reads, frac,
                os.path.join(outdir, os.path.basename(reads))))
    return cmd_list


def downsample_reads(args, reads_ns, config, new_cov=100, logger=None):
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
    ds_dir = os.path.join(args.tmp_dir, "seqtk_dir", "")
    os.makedirs(ds_dir)
    ds_cmds = make_seqtk_ds_cmd(args=args, reads_ns=reads_ns, config=config,
                                new_coverage=new_cov, outdir=ds_dir, logger=logger)
    for cmd in ds_cmds:
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    if args.fastq1 is not None:
        args.fastq1 = os.path.join(ds_dir, os.path.basename(args.fastq1))
    if args.fastq2 is not None:
        args.fastq2 = os.path.join(ds_dir, os.path.basename(args.fastq2))
    return new_cov


def make_bwa_cmds(args, config, outdir, ref, reads_ns, fastq1, fastq2):
    """ map with bwa
    we dont use the reference to args.fastq or args.fastq2 in case
    we are downsampling
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
                            "2> {2}samse.log > {2}mapping.sam").format(
                                config.bwa, ref, outdir, fastq1)
        cmd_list.append(sampe_cmd)

    else: # if long reads
        if fastq2 is None:
            mem_cmd = "{0} mem -t {1} -M {2} {3} > {4}mapping.sam 2> {4}bwa_mem.log".format(
                config.bwa, args.threads, ref, fastq1, outdir)
        else:
            mem_cmd = "{0} mem -t {1} -M {2} {3} {5} > {4}mapping.sam 2> {4}bwa_mem.log".format(
                config.bwa, args.threads, ref, fastq1, outdir, fastq2)
        cmd_list.append(mem_cmd)

    return (cmd_list, os.path.join(outdir, "mapping.sam"))


def make_samtools_cmds(exe, sam, outdir, out_bam):
    # removed the -q 10
    convert = str("{0} view  -Shb {1} | {0} sort - >" +
                  " {3}").format(
                      exe, sam, outdir, out_bam)
    index = str("{0} index {1} 2> {2}samtools_index.log").format(
                      exe, out_bam, outdir)
    return [convert, index]


def align_reads(ref, dirname, reads_ns,  downsample, args, config, logger):
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
    bwa_reference = os.path.join(bwa_dir, os.path.basename(ref))
    shutil.copyfile(ref, bwa_reference)
    if downsample:
        logger.debug("Downsampling reads for insert-size estimation...")
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
    logger.debug("BWA aligning reads vs reference...")
    bwa_cmds, mapping_sam = make_bwa_cmds(args, config, outdir=bwa_dir, ref=bwa_reference,
                                      reads_ns=reads_ns, fastq1=fastq1,
                                      fastq2=fastq2)
    sorted_bam = os.path.splitext(mapping_sam)[0] + ".bam"
    samtools_cmds =  make_samtools_cmds(exe=config.samtools,
                                        sam=mapping_sam,
                                        outdir=samtools_dir, out_bam=sorted_bam)
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


def run_picard_insert_stats(bam, config, picard_outdir, logger):  # pragma: no cover
    """ run picard cmd to get insert stats, and return path to the results
    """
    pic_cmd, statfile = make_picard_stats_command(bam, config=config,
                                                  picard_outdir=picard_outdir)
    logger.debug(pic_cmd)
    subprocess.run(pic_cmd,
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    return statfile


def parse_picard_insert_stats(statfile):
    """ look through a picard stats results, return mean and stddev of insert
    UPDATED 2018-10-05 syntax of picard has changed

    """
    # get_next, = False
    # with open(statfile, "r") as inf:
    #     for line in inf:
    #         if get_next:
    #             min_insert = line.split("\t")[3]
    #             max_insert = line.split("\t")[4]
    #             mean_insert = line.split("\t")[5]
    #             stddev_insert = line.split("\t")[6]
    #             break
    #         if "MEDIAN" in line:
    #             #  the next line will contain the results
    #             get_next=True
    data_dict = {}

    get_next = False
    with open(statfile, "r") as inf:
        for line in inf:
            if get_next :
                data = line.split("\t")
                break
            if "MEDIAN" in line:
                headers = line.split("\t")
                #  the next line will contain the data
                get_next = True
    for k, v in zip(headers, data):
        data_dict[k] = v
    return (data_dict['MEAN_INSERT_SIZE'], data_dict['STANDARD_DEVIATION'])



def get_insert_stats(bam, args, config, logger):  # pragma: no cover
    """
    converts contigs.sam -> bam, sorts, indexes and generates
    insert stats with Picard

    required params: $ (tmp directory);
                      $ (reference);

    returns        : $ (insert size)
                   : $ (stddev)
    """
    logger.debug("Collecting insert size statistics ")
    stat_dir = os.path.join(args.tmp_dir, "insert_stats", "")
    os.makedirs(stat_dir)
    statfile = run_picard_insert_stats(bam=bam, config=config, picard_outdir=stat_dir,
                                       logger=logger)
    mean_insert, stddev_insert = parse_picard_insert_stats(statfile)
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
        "__REFERENCE__": [results, "current_reference"],
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
    # CREATE = assembler['create_dir']
    # if CREATE:
    #     os.makedirs(os.path.join(args.tmp_dir, assembler['name']))
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


def get_contig_info(contigs, x):
    Ns, Ns_x = 0, 0
    count, count_x = 0, 0
    all_lengths, all_lengths_x = [], []
    if os.path.getsize(contigs) == 0:
        raise ValueError("contigs file is empty!")
    with open(contigs, "r") as inf:
        for rec in SeqIO.parse(contigs, "fasta"):
            length = len(rec.seq)
            Ns = Ns + rec.seq.lower().count("n")
            count = count + 1
            all_lengths.append(length)
            if length > x:
                all_lengths_x.append(length)
                Ns_x = Ns_x + rec.seq.lower().count("n")
                count_x = count_x + 1
    all_lengths.sort()
    all_lengths_x.sort()
    L50, N50 = get_L50_N50(all_lengths)
    L50_x, N50_x = get_L50_N50(all_lengths_x)
    return {"count": count, "count_x": count_x,
            "all_lengths": all_lengths,
            "all_lengths_x": all_lengths_x,
            "L50": L50, "N50": N50,
            "L50_x": L50_x, "N50_x": N50_x,
            "Ns": Ns, "Ns_x": Ns_x,
            #  I know this is bad, but if you have slashes in our names,
            # you had this coming
            "path": "/".join(contigs.split( "/")[-2:])}


def get_contig_stats(contigs, old_contigs=None):
    """
    Reports contig statistics on assembly. Reports on scaffolds or
    contigs depending upon 2nd argument passed - contigs gives values
    for all contigs and those >200xbp

    required params: $ (path to contigs)
                 $ ('scaffolds'|'contigs')

    returns          $ (0)
    """
    output = [
        [""],          #0
        ["count"],         #1
        ["Max Length"],    #2
        ["Assembly size"], #3
        ["L50"],           #4
        ["N50"],           #5
        ["N's"],           #6
        ["path"]           #7
    ]
    ctype= "seqs"
    if old_contigs is not None:
        results = get_contig_info(old_contigs, x=200)
        output[0].extend(["Old seqs", "Old seqs > %d" % 200])
        output[1].extend([results['count'], results['count_x']])
        output[2].extend([max(results['all_lengths']), max(results['all_lengths_x'])])
        output[3].extend([sum(results['all_lengths']), sum(results['all_lengths_x'])])
        output[4].extend([results['L50'], results['L50_x']])
        output[5].extend([results['N50'], results['N50_x']])
        output[6].extend([results["Ns"], results["Ns_x"]])
        output[7].extend([results['path'], ""])
    results = get_contig_info(contigs, x=200)
    output[0].extend(["New seqs", "New seqs > %d" % 200])
    output[1].extend([results['count'], results['count_x']])
    output[2].extend([max(results['all_lengths']), max(results['all_lengths_x'])])
    output[3].extend([sum(results['all_lengths']), sum(results['all_lengths_x'])])
    output[4].extend([results['L50'], results['L50_x']])
    output[5].extend([results['N50'], results['N50_x']])
    output[6].extend([results["Ns"], results["Ns_x"]])
    output[7].extend([results['path'], ""])
    return output


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
        logger.error(e, exc_info=True)
        logger.error("error splitting kmers in %s by comma!", k)
        logger.error(e)
        raise ValueError
    logger.debug(klist)
    new_ks = []
    for i in klist:
        if i > readlen:
            logger.warning("removing %d from list of kmers: exceeds read length",
                        i)
        elif readlen - i <= min_diff:
            logger.warning("removing %d from list of kmers: too close " +
                        "to read length", i)
        elif i % 2 == 0:
            logger.warning("removing %d from list of kmers: must be odd", i)
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
    logger.info("Starting %s assembly ... ", assembler)
    if conf_assembler['create_dir']:
        outdir = os.path.join(args.tmp_dir, conf_assembler['name'])
        os.makedirs(outdir)
    # here we execute the "run" function from the appropriate runner script
    sub_mains[assembler](getattr(config, assembler), assembler_cmd, config,
                         reads_ns, logger)

    logger.info("%s assembly statistics", assembler);
    report = get_contig_stats(contigs_path)
    logger.info("CONTIG STATS:\n" + tabulate.tabulate(report))
    # Only generate scaffold stats for paired read alignments...
    if os.path.exists(scaffolds_path) and args.fastq2 is not None:
        report = get_contig_stats( scaffolds_path)
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


def merge_assemblies(tool, args, config, reads_ns, logger):
    """
    combines two assemblies using the selected merge-method

    required params: $ (tmpdir)
                 $ (arrayref of assemblers used)
                 $ (method)
                 $ (reference)

    returns        : $ (0)
    """
    logger.info("Merging assemblies (%s)...", args.merger)
    cmd = tool['command']
    contig_output = tool['contig_output']
    if tool['create_dir']:
        outdir = os.path.join(args.tmp_dir, args.merger)
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

    report = get_contig_stats(contig_output )
    logger.info("\n\nMERGED ASSEMBLY STATISTICS:" +
                "\n=========================\n\n" + tabulate.tabulate(report))
    return contig_output


def check_id(ref, args, contigs, config, results, logger):
    """
    Checks wether the assembled contigs have sufficient identity to the
    reference for reference-based scaffolding

    required params: $ (tmp directory);

    returns        : $ id_ok (boolean)
    """
    logger.info("Checking identity of assembly with reference...")
    ID_OK = False
    check_dir = os.path.join(args.tmp_dir, "id_check_" +
                             os.path.splitext(os.path.basename(ref))[0], "")
    os.makedirs(check_dir)
    cmd = str("{0} -query {1} -subject {2} -outfmt 6 -evalue 0.01 -out " +
              "{3}blastout.tab 2>&1 > {3}blastn.log").format(
        config.blastn, contigs, ref, check_dir)
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
        percent_id = round(statistics.mean(fake_aves), 3)
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
    return ID_OK, percent_id


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
    if args.merger is None and len(args.assemblers) > 1:
        raise ValueError("Must provide merge method if using multiple assemblers")
    for thing in ["scaffolder", "finisher"]:
        if getattr(config, thing + "s") is None:
            continue
        for tool in getattr(config, thing + "s"):
            if tool['name'].lower() == getattr(args, thing):
                if len(args.references) == 0 and tool['ref_required']:
                    raise ValueError("%s %s requires a reference." % \
                                     (thing, getattr(args, thing)))
    # ensure exes are there for fastqc, seqtk
    if args.downsample != 0 and config.seqtk is None:
        raise ValueError("Seqtk is needed for downsampling")
    if "q" in args.stages.lower() and config.fastqc is None:
        raise ValueError("fastqc not found; remove 'q' from --stages")
    if "t" in args.stages.lower() and config.sickle_trim is None:
        raise ValueError("sickle-trim not found; remove 't' from --stages")
    if "s" in args.stages.lower() and args.scaffolder is None:
        raise ValueError("scaffolder arg empty; either provide " +
                         "scaffolder or remove 's' from --stages")
    if "f" in args.stages.lower() and args.finisher is None:
        raise ValueError("finisher arg empty; either provide " +
                         "finisher or remove 'f' from --stages")
    if "v" in args.stages.lower() and args.varcaller is None:
        raise ValueError("varcaller arg empty; either provide " +
                         "varcaller or remove 'v' from --stages")


def make_empty_results_object():
    """ returns a namepace containing both run history and paths of res files
    This is to keep the args object from being mutilated over time, and to
    Ensure that we have a way of looking back at a run to see what happened.
    """
    results = Namespace(
        organism=None,
        assemblers_list=None, # to be filled by match_assembler args
        assemblers_results_dict_list=[], # {name:None, contigs:None, scaffold:None}
        ID_OK=None, # to be set by check_id
        reference_percent_sim=None, # to be set by check_id
        current_agp=None,
        current_vcf=None,
        current_contigs=None,
        current_scaffolds=None,
        current_reference=None,
        current_contigs_source=None,
        current_scaffolds_source=None,
        current_embl=None, # run_prokka fills this in
        current_embl_source=None, # run_prokka fills this in
        old_embl=[(None, "init")],  # [path, source]
        old_contigs=[(None, "init")],  # [path, source]
        old_scaffolds=[(None, "init")],  # [path, source]
        old_references=[(None, "init")],  # [path, source]
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
        ["Initial Coverage", str(reads_ns.coverage) + "x"],
        ["Projected Coverage (Downsampled)", str(reads_ns.downsampled_coverage) + "x"],
        ["Projected Long Read Coverage", str(reads_ns.long_read_coverage) + "x"]
        ]
    logger.info("\n\nLIBRARY DETAILS:" +
                "\n=========================\n" +
                tabulate.tabulate(read_table) + "\n")

    run_table = [
        ["Assembly category",        reads_ns.lib_type],
        ["Selected assemblers",      " ".join(args.assemblers)],
        ["Selected assembly merger", args.merger],
        ["Selected scaffolder",      args.scaffolder],
        ["Selected finisher",        args.finisher],
        ["Selected variant caller",  args.varcaller],
        ["Trim QV",                  args.trim_qv],
        ["Trim length",              args.trim_length],
        ["Break Origin",             'b' in args.stages]
    ]
    logger.info("\n\nASSEMBLER DETAILS:" +
                "\n=========================\n" +
                tabulate.tabulate(run_table) + "\n")


def run_scaffolder(args, tools, config, reads_ns, results, run_id, logger=None):
    if tools.scaffolder['linkage_evidence'] == "paired-end":
        scaffs = run_pe_scaffolder(args, tools, config, reads_ns, results, run_id, logger=logger)
    else:
        scaffs = run_ref_scaffolder(args, tools, config, reads_ns, results, run_id, logger=logger)
    return scaffs


def run_ref_scaffolder(args, tools, config, reads_ns, results, run_id, logger=None):
    """
    Runs specified scaffolder....

    Scaffolders which don't separate unscaffolded contigs can be wrapped
    in a script which created a "$scaffolder.contig_ids" file listing the
    IDs of contigs scaffolded. These will then be used following scaffolding
    to create our own file of unscaffolded contigs
    """
    logger.info("Starting %s", tools.scaffolder['name'])
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
    with open(results.current_reference, "r") as ref:
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
    # if not len(ref_list) > 1:
    #     # if we don't have mulitple references, we just need to make the
    #     # reference and contigs available under consistent names
    #     reference = ref_list[0]
    #     # contigs = contigs_output
    # else:
    if True:
        reference_copy = os.path.join(run_dir, "reference_for_scaffolding.fasta")
        contigs_copy = os.path.join(run_dir, "contigs_to_be_scaffolded.fasta")
        shutil.copyfile(results.current_reference, reference_copy)
        shutil.copyfile(results.current_contigs, contigs_copy)
        # blast_cmd = str("{0} -query {1} -subject {2}  -task blastn -out " +
        #                 "{3}clusters.blast 2>&1 > {3}blastn.log").format(
        #                     config.blastn, contigs_copy, reference_copy, run_dir)
        blast_cmd = str("{0} -query {1} -subject {2} -outfmt 6 -out " +
                  "{3}clusters.tab 2>&1 > {3}blastn.log").format(
                      config.blastn, contigs_copy, results.current_reference, run_dir)
        # blast_cmd3 = str("{0} -query {1} -subject {2} -outfmt 5 -evalue 0.01 -out " +
        #           "{3}clusters.xml 2>&1 > {3}blastn.log").format(
        #               config.blastn, contigs_copy, args.references, run_dir)
        for cmd in [blast_cmd]:
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
            subscaff_dir = os.path.join(run_dir, "SIS_" + ref_id, "")
            os.makedirs(subscaff_dir)
            merge_these_scaffolds.append(os.path.join(subscaff_dir, "scaffolds.fasta"))
            # here we execute the "run" function from the appropriate runner script
            sub_mains[args.scaffolder](config=config, args=args, results=results,
                                       contigs=scaff_contigs,
                                       ref=ref_path,
                                       scaff_dir=subscaff_dir,
                                       logger=logger)
            resulting_scaffold = os.path.join(run_dir, "scaffolds.fasta")
        else:
            raise ValueError("%s not found, or empty!" % scaff_contigs)
    counter = 1
    with open(resulting_scaffold, "w") as outf:
        for scaf_path in merge_these_scaffolds:
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

    report = get_contig_stats(resulting_scaffold, old_contigs=results.current_contigs);
    logger.info ("\n\nScaffolded assembly stats" +
                 "\n=========================\n" +
                 tabulate.tabulate(report) + "\n")
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
    with open(results.current_reference, "r") as ref:
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
        sorted_bam = align_reads(ref=results.current_reference, dirname="align_to_contigs", reads_ns=reads_ns,
                                 config=config, downsample=True, args=args,
                                 logger=logger)
        reads_ns.insert_mean, reads_ns.insert_stddev = get_insert_stats(
            bam=sorted_bam, config=config, args=args, logger=logger)

    run_scaffold_output = os.path.join(args.tmp_dir,
                                       args.scaffolder + "_no_reference",
                                       scaffold_output)

    # moved scaf_args out of if statement as well!
    exec_cmd = exec_cmd + " 2>&1 > " + os.path.join(
        args.tmp_dir,
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
                    SeqIO.write(rec, outf, "fasta")
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
def parse_origin_from_coords(coords, flex, reference, logger):
    """    returns a dict {seqID: "contig_name:ori",...}
    """
    pattern = re.compile("\s+")
    origin_dict = {}
    with open(reference, "r") as inf:
        for rec in SeqIO.parse(inf, "fasta"):
            FOUND = False
            with open(coords, "r") as in_coords:
                for line in in_coords:
                    line = line.strip()
                    line = line.replace("|", "")
                    fields = pattern.split(line)
                    # if we haven't found one already
                    if not FOUND:
                        #  and if the start is near 1 and on the correct record
                        if int(fields[0]) <= flex and fields[7] == rec.id :
                            origin_dict[rec.id] = "{0}:{1}".format(fields[8], fields[2])
                            FOUND = True
                            logger.info("Potential origin found at %s",
                                        origin_dict[rec.id] )
    return origin_dict

def find_origin(config, args, results, ori_dir, flex, logger):  # pragma: no cover
    """ attempts to find the origin for each record in the reference
    returns a dict {seqID: "contig_name:ori",...}
    """
    logger.info("Attempting to identify origin...");
    cmds, ori_coords =  make_nucmer_delta_show_cmds(config=config, ref=results.current_reference,
                                       query=results.current_scaffolds,
                                       out_dir=ori_dir, prefix="ori", header=False)
    run_nucmer_cmds(cmds, logger)
    origin_dict = parse_origin_from_coords(coords=ori_coords, flex=flex,
                                           reference=results.current_reference,
                                           logger=logger)
    return origin_dict


def split_origin(origin_dict, ori_dir, scaffolds, min_len, logger):
    """ given a origin dict {seq_id, contig:pos}, split seqence there and write
    """
    splitIO = os.path.join(ori_dir, "scaffolds_ori_split.fasta")  # SplitIO
    counter = 0
    if len(origin_dict) == 0:
        logger.warning("no origin found! Continuing...")
        return scaffolds
    # for ref_id, origin in sorted(origin_dict.items())
    contigs_containing_scaffolds = \
        [x.split(":")[0] for x in origin_dict.values()]
    with open(scaffolds, "r") as inscaff, open(splitIO, "w") as splitscaff:
        for rec in SeqIO.parse(inscaff, "fasta"):
            counter = counter + 1
            if rec.id in contigs_containing_scaffolds:
                # get the proper entry from the dictionary, its a bit clunky
                origin = [x for x in origin_dict.values() if
                          x.split(":")[0] == rec.id][0]
                ori_scaffold, pos = origin.split(":")
                logger.debug("Splitting %s at origin %s", ori_scaffold, pos)
                # write out both parts after splitting
                part_a = rec.seq[0: int(pos)]
                # dont do +1 because slicing in not inclusive
                # part_b = rec.seq[int(pos) + 1: ]
                part_b = rec.seq[int(pos): ]
                ori_a = SeqRecord(id="Scaffold_%05d" % counter,
                                  description="{} part A split at {}".format(
                                      rec.id , int(pos)),
                                  seq=part_a)
                counter = counter + 1
                ori_b = SeqRecord(id="Scaffold_%05d" % counter,
                                  description="{} part B split at {}".format(
                                      rec.id , int(pos)),
                                  seq=part_b)
                for ori in [ori_a, ori_b]:
                    if not len(ori.seq) < min_len:
                        SeqIO.write(ori, splitscaff, "fasta")
                    else:
                        logger.warning("the following sequence was ignored as it was too short after spliting at origin:")
                        logger.warning(ori)
            else:
                rec.id = "Scaffold_%05d" % counter
                SeqIO.write(rec, splitscaff, "fasta")
    report = get_contig_stats(splitIO, old_contigs=scaffolds);
    logger.info ("\n\nScaffolded assembly stats post - breaking origin" +
                 "\n=========================\n" +
                 tabulate.tabulate(report) + "\n")
    return splitIO


def find_and_split_origin(args, config, results, tools, reads_ns, logger):
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

    origin_dict = find_origin(config, args, results, ori_dir, flex=10, logger=logger)
    outIO = os.path.join(ori_dir, "split_ori_scaff.fasta")  # OutIO
    splitIO = split_origin(origin_dict, ori_dir, min_len=10,
                           scaffolds=results.current_scaffolds, logger=logger)

    #                 #  now, rerun scaffolder on our split origins
    #                 results.old_scaffolds.append([results.current_scaffolds,
    #                                               results.current_scaffolds_source])
    #                 results.current_scaffolds = outIO
    #                 results.current_scaffolds_source = "split_ori_scaff"
    #                 if tools.scaffolder:
    #                     new_scaffolds = run_scaffolder(
    #                         args=args, config=config, tools=tools,
    #                         reads_ns=reads_ns, results=results,
    #                         run_id=2, logger=logger)
    #                     results.old_scaffolds.append([results.current_scaffolds,
    #                                                   results.current_scaffolds_source])
    #                     results.current_scaffolds = new_scaffolds
    #                     results.current_scaffolds_source = "find_origin_rescaffolded"
    #                 else:
    #                     new_scaffolds = outIO
    #                     results.old_scaffolds.append(
    #                         [results.current_scaffolds,
    #                          results.current_scaffolds_source])
    #                     results.current_scaffolds = new_scaffolds
    #                     results.current_scaffolds_source = "find_origin_not_scaffolded"
    #                 with open(new_scaffolds, "r") as infile, open(outIO, "a") as outfile:
    #                     for new_rec in SeqIO.parse(infile, "fasta"):
    #                         SeqIO.write(new_rec, outfile, "fasta")
    #             else:
    #                 # these dont need splitting
    #                 with open(outIO, "a") as outfile:
    #                     SeqIO.write(rec, outfile, "fasta")

    # #renumber scaffolds to ensure they are unique...
    # counter = 1
    # renumbered_scaffs = os.path.join(
    #     os.path.dirname(splitIO), "scaffold.fasta")
    # with open(splitIO, "r") as inf, open(renumbered_scaffs, "w") as outf:
    #     for rec in SeqIO.parse(inf, "fasta"):
    #         rec.id = "scaffold_%05d" % counter
    #         counter = counter + 1
    #         SeqIO.write(rec, outf, "fasta")
    # results.old_scaffolds.append([results.current_scaffolds, results.current_scaffolds_source])
    # results.current_scaffolds = renumbered_scaffs
    # results.current_scaffolds_source = "renumbering_after_origin_split"


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


def order_scaffolds(args, config, results, logger):
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
    cmds, orient_coords =  make_nucmer_delta_show_cmds(config=config, ref=results.current_reference,
                                       query=results.current_scaffolds,
                                       out_dir=orient_dir, prefix="ori2", header=False)
    run_nucmer_cmds(cmds, logger=logger)
    orient_count = { '+': 0,
                     '-': 0 }
    with open(orient_coords, "r") as inf:
        contig = '';
        for line in inf:
            line = line.strip().replace("|", "")
            fields = re.split("\s+", line)
            if contig != fields[8]:
                if contig != "":
                    orientations = store_orientation(contig, orient_count )
            contig = fields[8]
            start  = int(fields[2])
            end    = int(fields[3])
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
    orientations = store_orientation(contig, orient_count )

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
            rec.id = "scaffold_%05d"  % str(i + 1)
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
    with open(results.current_reference, "r") as ref:
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


def store_orientation(contig, or_count):
    """
    save correct contig orientation in hash passed as first parameter

    since alignment fof scaffold can contain blocks in reverse orientation
    need to use the 'prevailing' orientation based on number of +/- blocks

    required parameters: $ (orientations hash)
                     $ (contig)
                     $ (hash of +/- base counts)

    returns           : none
    """
    orientations = {}
    if ((or_count['+'] and or_count['-']) and (or_count['+'] > or_count['-'])) \
       or (or_count['+'] and not or_count['-']):
        orientations[contig] = '+'
    else:
        orientations[contig] = '-'


def run_finisher(args, config, reads_ns, tools, results, logger):
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
    logger.info("Finishing assembly with %s", tools.finisher['name'])
    finisher_dir = os.path.join(args.tmp_dir, tools.finisher['name'], "")
    os.makedirs(finisher_dir)
    finished_scaffolds = os.path.join(finisher_dir, tools.finisher['output_scaffolds'])
    # here we execute the "run" function from the appropriate runner script
    sub_mains[tools.finisher['name']](config=config, args=args, results=results,
                               reads_ns=reads_ns,
                               scaffolds=results.current_scaffolds,
                               finisher_dir=finisher_dir,
                               logger=logger)
    # resulting_scaffold = os.path.join(run_dir, "scaffolds.fasta")
    # if tool['create_dir']:
    #     os.makedirs(outdir)
    # cmd = replace_placeholders(
    #     string=tool['command'],
    #     config=config,
    #     reads_ns=reads_ns, args=args)
    # logger.debug("running the following command:\n %s". cmd)
    # subprocess.run(cmd,
    #                shell=sys.platform != "win32",
    #                stdout=subprocess.PIPE,
    #                stderr=subprocess.PIPE,
    #                check=True)
    report = get_contig_stats(finished_scaffolds, results.current_scaffolds )
    logger.info ("\n\nGenome finishing with " + tools.finisher['name'] +
                 " statistics:\n"
                 "=========================\n" +
                 tabulate.tabulate(report) + "\n")
    update_results(results=results, thing="scaffolds", path=finished_scaffolds,
                   source=tools.finisher['name'] )


def build_agp(args, results, reads_ns, evidence, logger):
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

    see https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
    """

    if evidence not in ['paired-ends', 'align_genus', 'align_xgenus']:
        raise ValueError("Unknown evidence type: %s", evidence)
    logger.info("Creating AGP file...");
    agp_dir = os.path.join(args.tmp_dir, "agp")
    agpfile_path = os.path.join(agp_dir, "scaffolds.agp")
    os.makedirs(agp_dir)
    lines = []
    lines.append( "##agp-version 2.0")
    lines.append("#{}".format(results.organism))
    contigs_path = os.path.join(agp_dir, "contigs.fasta")
    scaffolds_path = os.path.join(agp_dir, "scaffolds.fasta")
    contig_count = 0
    scaffold_count = 0
    gap_dict = {}    #overall per-scaffold gaps to return....
    with open(results.current_scaffolds, "r") as scaffold_inIO:
        for scaffold_count, scaffold in enumerate(SeqIO.parse(scaffold_inIO, "fasta")):
            scaffold_id = "scaffold_%06d" % (scaffold_count + 1)
            contig_start, contig_end, gap_start, gap_end = 0, None, None, None

            # location tracking for included contigs/gaps
            # changed from a dict to a list, cause we like order
            # [id, coords_string, legnth]
            contigs_list = []
            gaps =  []
            # this is kind of crude, but seems to work...
            for idx, base in enumerate(list(scaffold.seq)):
                # print("{}, {}, gs: {}; ge: {}".format(
                #     idx, base, gap_start, gap_end ))
                # if a normal base and it isnt immediately following a gap,
                # go to the next base, keeping track of where we were



                # if this base is an N
                if base == "N":
                    # and we are not in a gap yet
                    if gap_start is None:
                        # mark the end of the old contig (if needed)
                        # and the start of the gap
                        contig_end = idx
                        gap_start = idx
                    else:  # ie, a gap has already been opened
                        gap_end = idx
                else: # base is not an N
                    if gap_start is None:
                        # if this isnt first base after a gap
                        contig_end = idx
                    else:
                        # we just exited a gap!!!
                        # report the contig lengths, the gaps, etc
                        gap_end = idx
                        gap_length = gap_end - gap_start

                        if gap_length < 10 and args.mode == 'submission':
                            # skip gaps <10 bp which are acceptable by EMBL
                            gap_start = None
                        elif args.mode == 'draft' and gap_length <= 1:
                            # we can leave single ambiguous bases alone
                            gap_start = None
                        else:
                            gap_end = idx
                            if contig_end - contig_start < 200 and args.mode == "submission":
                                #need to extend previous gap to new position including short contig
                                if len(gaps) != 0:
                                    last_gap = gaps[len(gaps)]
                                    last_start, last_end = last_gaps.split("-")
                                    new_gap = "{0}-{1}".format(last_start, gap_end)
                                    gaps.append(new_gap)
                                    gap_start = None
                                    contig_start = idx
                                else:
                                    pass
                            else:
                                contig_count = contig_count + 1
                                contig_id = "contig_%06d" % contig_count
                                with open(contigs_path, "a") as  outf:
                                    SeqIO.write(
                                        SeqRecord(
                                            scaffold.seq[contig_start: contig_end],
                                            id=contig_id),
                                        outf, "fasta")
                                contigs_list.append(
                                    [contig_id,
                                     "{}-{}".format(contig_start, contig_end),
                                    contig_end - contig_start
                                     ])
                                gap = "{}-{}".format(gap_start, gap_end)
                                gaps.append(gap)
                                gap_start = None
                                contig_start = idx
            # Output last contig from last contig_start position to scaffold end if it is longer than
            # 200 bp and we are running in submission mode, otherwise remove the last gap to truncate the scaffold...
            if len(scaffold.seq) - contig_start  > 200 or args.mode == 'draft' or contig_count == 0:
                contig_count = contig_count + 1
                contig_id = "contig_%06d" % contig_count
                with open(contigs_path, "a") as  outf:

                    SeqIO.write(SeqRecord(scaffold.seq[contig_start: len(scaffold.seq)], id=contig_id), outf, "fasta")
                contigs_list.append([contig_id,
                                     "{}-{}".format(contig_start, contig_end),
                                     contig_end - contig_start])
            else:
                gaps.pop()

            gap_dict[scaffold_id] = gaps
            scaffold_part = 0

            #write AGP output and new scaffolds fasta file
            record_dict = SeqIO.index(contigs_path, "fasta")
            scaffold_seq = ""
            if len(contigs_list) > 0:
                # if ( $#contig_ids > -1 ) {
                for i, (contig_id, coords, length) in enumerate(contigs_list):
                    if contig_id != "" :  # dont know when this could possibly be true
                        contig = record_dict[contig_id]
                        contig_start, contig_end  = coords.split("-")
                        scaffold_part = scaffold_part + 1
                        lines.append(
                            "{ob}\t{ob_beg}\t{ob_end}\t{part_number}\t{comp_type}\t{comp_id}\t{comp_beg}\t{comp_end}\t{orientation}".format(
                                ob=scaffold_id,
                                ob_beg=int(contig_start) + 1,
                                ob_end=int(contig_end) + 1,
                                part_number=scaffold_part,
                                comp_type="W",
                                comp_id=contig_id,
                                comp_beg=1,
                                comp_end=length,
                                orientation="+"))
                        scaffold_seq = scaffold_seq + str(contig.seq)
                        if i < len(contigs_list) - 1:
                            gap_start, gap_end = gaps[i].split("-")
                            gap_size = int(gap_end) - int(gap_start)
                            scaffold_part = scaffold_part + 1
                            lines.append(
                                "{ob}\t{ob_beg}\t{ob_end}\t{part_number}\t{comp_type}\t{gap_length}\t{gap_type}\t{linkage}\t{linkage_ev}".format(
                                    ob=scaffold_id,
                                    ob_beg=int(gap_start) + 1,
                                    ob_end=int(gap_end) + 1,
                                    part_number=scaffold_part,
                                    comp_type="N",
                                    gap_length=gap_size,
                                    gap_type="scaffold",
                                    linkage="yes",
                                    linkage_ev=evidence))
                            scaffold_seq = scaffold_seq + ("N" * gap_size)
            if scaffold_seq is not None:
                with open(scaffolds_path, "a") as outf:
                    SeqIO.write(SeqRecord(Seq(scaffold_seq), id=scaffold_id),
                                outf, "fasta")
            record_dict.close()
        with open(agpfile_path, "w") as outf:
            for line in lines:
                outf.write(line + "\n")
    results.current_agp = agpfile_path
    return gaps


def get_prokka_cmd(exe, outdir, args, seqs):
    return str(
        "{0} --addgenes --outdir {1} --prefix prokka --genus {2} " +
        "--species {3} --strain {4} --locustag {5} --centre {6} " +
        "--cpus {7} {9} > {8}prokka.log 2>&1").format(
            exe,           #0
            outdir,        #1
            args.genus,    #2
            args.species,  #3
            args.strain,   #4
            args.locustag, #5
            args.centre,   #6
            args.threads,  #7
            args.tmp_dir,  #8
            seqs)


def make_embl_from_gbk(gbk, output_file):
    with open(gbk, "r") as inf, open(output_file, "w") as outf:
        for rec in SeqIO.parse(inf, "genbank"):
            SeqIO.write(rec, outf, "embl")


def run_prokka(config, args, results, logger):
    """
    generates annotation on the assembly using prokka

    required parameters: $ (tmpdir)
		       $ (genus)
                     $ (species)
                     $ (strain)
                     $ (locustag)
                     $ (centre)

    returns            : $ (0)
    """
    logger.info("Starting PROKKA...")
    #use scaffolds if we have them, otherwise contigs....
    if results.current_scaffolds is not None:
        seqs = results.current_scaffolds
    else:
        update_results(results, thing="scaffolds",
                       path=results.current_contigs,
                       source="current_contigs, pre prokka" )
        seqs = results.current_contigs

    prokka_dir = os.path.join(args.tmp_dir, "prokka", "")
    # os.makedirs(prokka_dir)
    cmd = get_prokka_cmd(exe=config.prokka, outdir=prokka_dir, args=args, seqs=seqs)
    logger.debug("running the following command:\n %s", cmd)
    subprocess.run(cmd,
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    update_results(results, thing="embl",
                   path=os.path.join(prokka_dir, "prokka.embl"), source="prokka" )
    update_results(results, thing="scaffolds",
                   path=os.path.join(prokka_dir, "prokka.fna"), source="prokka" )
    return prokka_dir


def update_results(results, thing, path, source ):
    """ this bumps a current thing to a list of old things,
    updating the attribute with a new thing

    """
    assert thing in ["embl", "contigs", "scaffolds", "reference"], \
        "cant update %s" % thing
    # if the current thing's path is  None, this is the first entry, so
    # don't need to bump the old value to the appropriate list
    if getattr(results, "current_" + thing) is not None:
        getattr(results, "old_" + thing).append(
            (getattr(results, "current_" + thing),
             getattr(results,"current_" + thing + "_source")))
    setattr(results, "current_" + thing, path)
    setattr(results, "current_" + thing + "_source", source)


def amosvalidate(args, results, config, reads_ns, logger):
    """
    Uses abyss's abyss-samtoafg script to convert spades sam and contigs
    to an amos bank

    required params: $ (tmpdir)
                 $ (insert size)
                 $ (insert size stddev)

    returns        : $ (0)
    """
    amos_dir = os.path.join(args.tmp_dir, "amos", "")
    os.makedirs(amos_dir)
    seq_file = results.current_contigs
    if results.current_scaffolds is not None:
        seq_file = results.current_scaffolds

    out_sequences = os.path.join(amos_dir, os.path.basename(seq_file))
    # copy sequences with a newline separating records
    with open(results.current_scaffolds, "r") as inf, open(out_sequences, "w") as outf:
        for rec in SeqIO.parse(inf, "fasta"):
            SeqIO.write(rec, outf, "fasta")
            outf.write("\n")
    assert reads_ns.insert_mean is not None and \
        reads_ns.insert_stddev is not None, "how did we get here with paired reads but no insert size calcualted?"
    logger.info("Converting to amos bank...");
    sam2afg_cmd = "{exe}  -m {ins_size} -s {ins_stddev} {infile} {sam} > {outdir}amos.afg".format(
        exe=config.sam2afg,
        ins_size=reads_ns.insert_mean,
        ins_stddev=reads_ns.insert_stddev,
        infile=out_sequences,
        sam=os.path.join(args.tmp_dir, "bwa",
                         os.path.splitext(os.path.basename(seq_file))[0]),
        outdir=amos_dir)

    bnk_cmd = "{exe} -cb {outdir}assembly.bnk -m {outdir}amos.afg  > {outdir}bank-transact.log 2>&1".format(
        exe=config['bank_transact'],
        outfile=amos_dir)
    logger.info("Running amosvalidate");

    # read amosvalidate script and comment out '4xx' lines, which run SNP checks and are extreeeeemly slow....
    new_exe = os.path.join(amos_dir, "amosvalidate")
    regex = re.compile(r"^4.*")
    with open(config.amosvalidate, "r") as source, open(new_exe, "w") as new_source:
        for line in source:
            line = regex.sub("#4", line)  # double check this line actully works
            new_source.write(line)
    os.chmod(new_exe, "0755")  # im twitching too...
    validate_cmd = "{exe} {outdir}assembly.bnk > {outdir}amosvalidate.log 2>&1".format(
        exe=new_exe,
        outdir=amos_dir)
    mapping_cmd = \
    "{exe} -i -p -b {outdir}/assembly.bnk CTG 2> /dev/null | cut -f 2,3 > {outdir}/ctg_to_iid.txt".format(
        exe=config.bank_report,
        outdir=amos_dir)

    for cmd in [sam2afg_cmd, bnk_cmd, validate_cmd, mapping_cmd]:
        logger.debug("running the following command:\n %s". cmd)
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)


def summarise_amosvalidate():
    """
    postprocesses amosvalidate outputs to make more readily digestible

    required params: $ (tmpdir)

    returns        : $ ($ - hashref to parsed results)
    """
    logger.info("processing amosvalidate results...");
    amosvalidate = os.path.join(args.tmp_dir, "amos", "")
    iid_mapping = os.path.join(args.tmp_dir, "amos", "ctg_to_iid.txt")
    results_dict = {}

    # Read Contig->iid mapping into a hash
    contig_to_iid_dict = {}
    pattern = re.compile(r"\s+")
    with open(iid_mapping, "r") as inf:
        for line in inf:
            contig, iid  = pattern.split(line)
            iid = iid.strip()
            contig_to_iid_dict[iid] = contig
            results_dict[contig]    = []

    outputs = glob.glob(os.path.join(amosvalidate, "*feat"))
    with open(os.path.join(amosvalidate, "assembly.all.feat")) as ALL:
        for line in ALL:
            fields    = patter.split(line)
            contig_id = contig_to_iid_dict[fields[0]]
            start     = fields[2]
            end       = fields[3]
            tpe      = fields[4]
            res_arr   = results[contig_id]
            #no idea how it ends up with a negative start...but it does...
            start = 0 if start < 0  else start
            #  making this a list?
            res_arr = [start, end, tpe]
            results[contig_id] = res_arr
    return results_dict


def finish_up(results, t0, all_stages, logger=None, args=None):
    """ print some summarizing information when we're done
    """
    logger.info("\n\nFinished BugBuilding!")
    return_results(args, results)
    stages_string = " - "
    for s in list(args.stages):
        stage_name = [x[1] for x in all_stages if x[0] == s][0]
        stages_string = stages_string + str(stage_name) + "\n - "
    logger.info("\n\nStages Executed:\n%s", stages_string)
    if "a" in args.stages:
        report_config = get_contig_stats(results.current_contigs)
        report_scaffold = get_contig_stats(results.current_scaffolds)
        logger.info ("\n\nFinal Contig Statistics:\n"
                     "=========================\n" +
                     tabulate.tabulate(report_config) + "\n")
        logger.info ("\n\nFinal Scaffold Statistics:\n"
                    "=========================\n" +
                     tabulate.tabulate(report_scaffold) + "\n")
    logger.info("Done building your Bug!")
    logger.info("Elapsed time: %.2fm" % ((time.time() - t0) / 60))
    sys.exit(0)


def finished_desired_stages(this_stage, args, all_stages, logger):
    """ Check if any remaining stages need to be/can be done
    """
    found_this = False
    # run through each of the stages, marking where we are now
    for s in all_stages:
        if this_stage == s[0]:
            found_this = True
        else:
            # for future stages, are they requested via args.stages?
            # and if so, have we completed the prereques
            if found_this:
                if s[0] in args.stages:
                    if s[3] is not None:
                        prereqs = list(s[3])
                        prereqs_requested = [x for x in prereqs if
                                             x in list(args.stages)]
                        if prereqs  == prereqs_requested:
                            return False
                        else:
                            logger.warning(str(
                                "Cannot run stage {0} or later stages " +
                                "becauseit requires stages {1}").format(
                                    s[1], s[3]))
                            return True
                    else: # if no prereqs
                        return False

    return True


def finish_if_done(this_stage, args, all_stages, logger, t0, results):
    if finished_desired_stages(
            this_stage=this_stage,
            args=args,
            all_stages=all_stages, logger=logger):
        finish_up(results=results, t0=t0, logger=logger,
                  args=args, all_stages=all_stages)


def make_gene_report(results, logger):
    if results.current_embl is None:
        return 1
    cds, tRNA, rRNA, gap = 0, 0, 0, 0
    with open(results.current_embl, "r") as inf:
        for rec in SeqIO.parse(inf, "embl"):
            for feature in rec.features:
                if feature.type=="CDS":
                    cds = cds + 1
                if feature.type=="assembly_gap":
                    gap = gap + 1
                if feature.type=="tRNA":
                    tRNA = tRNA + 1
                if feature.type=="rRNA":
                    rRNA = rRNA + 1
    genes = [
        ["Feature Type", "Number"],
        ["CDS", cds],
        ["tRNA", tRNA],
        ["rRNA", rRNA],
        ["gap", gap]
    ]
    logger.info("\n\nANNOTATION SUMMARY:" +
                "\n=========================\n" +
                tabulate.tabulate(genes) + "\n")
    return 1

def make_comparison_cmds(config, query, results, ref_id, comp_dir):
    cmds = []
    prefix = os.path.join(comp_dir, ref_id)
    cmds.append(str(
        "{config.blastn} -query {query} -subject {results.current_reference}" +
        " -out {prefix}_comparison.blastout -outfmt 6 > " +
        "{comp_dir}blast.log").format(**locals())
    )
    # also build a mummerplot in png format...
    cmds.append(str(
        "{config.nucmer} --prefix {prefix} {results.current_reference} " +
        "{query} > {comp_dir}nucmer.log "+"2>&1").format(**locals()))
    cmds.append(str(
        "{config.mummerplot} -large --filter --layout -p {prefix} -t png " +
        "-R {results.current_reference} -Q {query} {prefix}.delta  > " +
        "{comp_dir}mummerplot.log 2>&1").format(**locals()))
    cmds.append(str(
        "{config.mummerplot} -large --filter --layout -c -p {prefix}cover -t png " +
        "-R {results.current_reference} -Q {query} {prefix}.delta  > " +
        "{comp_dir}mummerplot_pip.log 2>&1").format(**locals()))
    return cmds


def build_comparisons(args, config, results, logger):
    """
    generates a comparison appropriate for viewing with ACT and a
    MUMmerplot to provide a quick overview

    required params: $ (tmpdir)
                     $ (reference)
                     $ (organism)

    returns        : $ (0)

    """
    logger.info("Building comparison plots with mummerplot")
    comp_dir = os.path.join(args.tmp_dir, "comparisons", "")
    os.makedirs(comp_dir)

    query = results.current_scaffolds
    comp_cmds = make_comparison_cmds(
        config=config, query=query,
        results=results,
        ref_id=os.path.splitext(
            os.path.basename(results.current_reference))[0],
        comp_dir=comp_dir)
    for cmd in comp_cmds:
        logger.debug(cmd)
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)


def return_results(args, results):
    """
    Copies results back from tmpdir, returning full working directory
    if dircopy argument specified

    required params: $ (tmpdir)
        $ (strain)
        $ (dircopy - flag to indicate entire directory
        should be returned)
        $ (mode - draft mode also needs amos bank copying)

    returns        : $ (0)
    """
    # get the ones in the results object
    targets = [
        "current_scaffolds", "current_embl", "current_contigs",
        "current_agp", "current_vcf"]
    target_paths = []
    for t in targets:
        if t is not None:
            tfile = getattr(results, t)
            target_paths.append(tfile)

    for t in target_paths:
            if os.path.exists(tfile):
                shutil.copyfile(
                    tfile,
                    os.path.join(args.outdir, os.path.basename(tfile)))
            else:
                raise FileNotFoundError
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


def run_cgview(results, args, config, logger):
    """
    run_cgview

    Runs cgview to generate a genome map from the annotations

    Required parameters: $ (tmpdir)

    Returns: $ (0)
    """
    cgview_dir = os.path.join(args.tmp_dir, "cgview", "")
    os.makedirs(cgview_dir)

    logger.info("Creating genome visualisaton...");
    embl = results.current_embl
    if embl is None:
        raise ValueError("No EMBL file currently defined")
    outfile = os.path.join(cgview_dir, "cgview.png")
    cmds = []
    cmds.append(str(
        "{config.perl} {config.cgview_xml_builder} -sequence {embl} " +
        "-output {cgview_dir}scaffolds_cgview.xml -gc_skew T > " +
        "{cgview_dir}xml_creator.log 2>&1").format(**locals()))
    cmds.append(str(
        "{config.cgview} -f png -i {cgview_dir}scaffolds_cgview.xml -o {outfile} > " +
        "{cgview_dir}cgview.log 2>&1").format(**locals()))
    for cmd in cmds:
        logger.debug(cmd)
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)


def run_varcaller(args, results, reads_ns, config, tools, logger):
    """
    Carries out variant calling using requested variant caller

    required params: $ (tmpdir)
                 $ (varcall)
                 $ (no. threads)
                 $ (read length)

    """

    # my $tmpdir      = shift;
    # my $varcall     = shift;
    # my $threads     = shift;
    # my $read_length = shift;

    logger.info("Running variant calling %s", tools.varcaller['name'])
    varcall_dir = os.path.join(args.tmp_dir, "varcaller_" + tools.varcaller['name'], "")
    os.makedirs(varcall_dir)

    map_cmds = []
    bwa_index = str("{config.bwa} index {results.current_reference} > " +
                    "{varcall_dir}bwa_index.log 2>&1").format(**locals())

    # Use bwa-bwt for 'short' reads less than 100 bp, and bwa-mem for longer reads
    if reads_ns.read_length_mean <= 100:
        bwa_map1 = str(
            "{config.bwa} aln -t {args.threads} " +
            "{results.current_reference} {args.fastq1} > " +
            "{varcall_dir}read1.sai 2> {varcall_dir}bwa_aln1.log").format(
                **locals())
        bwa_map2 =  str(
            "{config.bwa} aln -t {args.threads} " +
            "{results.current_reference} {args.fastq2} > " +
            "{varcall_dir}read2.sai 2> {varcall_dir}bwa_aln2.log").format(
                **locals())
        bwa_aln_single = str(
            "{config.bwa} sampe {results.current_reference} " +
            "{varcall_dir}read1.sai {args.fastq1} " +
            "2> {varcall_dir}sampe.log > {varcall_dir}reference.sam").format(
                **locals())
        bwa_aln_paired = str(
            "{config.bwa} sampe {results.current_reference} " +
            "{varcall_dir}read1.sai {varcall_dir}read2.sai  " +
            "{args.fastq1} {args.fastq2} " +
            "2> {varcall_dir}sampe.log > {varcall_dir}reference.sam").format(
                **locals())
        if reads_ns.paired:
            bwa_cmds = [bwa_index, bwa_map1, bwa_map2, bwa_aln_paired]
        else:
            bwa_cmds = [bwa_index, bwa_map1, bwa_aln_single]
    else:
        bwa_mem_single = str(
            "{config.bwa} mem -t {args.threads} {results.current_scaffolds} " +
            "{args.fastq1} {args.fastq2} > {varcall_dir} reference.sam "+
            "2> {varcall_dir}bwa_mem.log").format(**locals())
        bwa_mem_paired = str(
            "{config.bwa} mem -t {args.threads} {results.current_scaffolds} " +
            "{args.fastq1} {args.fastq2} > {varcaller_dir}reference.sam "+
            "2> {varcall_dir}bwa_mem.log").format(**locals())
        if reads_ns.paired:
            bwa_cmds = [bwa_index, bwa_mem_paired]
        else:
            bwa_cmds = [bwa_index, bwa_mem_single]

    samtools_view_cmd = str(
        "{config.samtools} view -q 10 -Sb " +
        "{varcall_dir}reference.sam 2> " +
        "{varcall_dir}samtoolsview.log >" +
        "{varcall_dir}reference.bam").format(
            **locals())
    samtools_sort_cmd = str(
        "{config.samtools} sort {varcall_dir}reference.bam -o " +
        "{varcall_dir}reference.sorted.bam > " +
        "{varcall_dir}samtools_sort.log 2>&1").format(**locals())
    samtools_index_cmd = str(
        "{config.samtools} index " +
        "{varcall_dir}reference.sorted.bam").format(
            **locals())
    map_cmds = bwa_cmds
    map_cmds.append(samtools_view_cmd)
    map_cmds.append(samtools_sort_cmd)
    map_cmds.append(samtools_index_cmd)
    for cmd in map_cmds:
        logger.debug(cmd)
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    sub_mains[tools.varcaller['name'] + "_var"](
        config=config, args=args, results=results,
        reads_ns=reads_ns,
        reference=results.current_reference,
        reference_bam=os.path.join(varcall_dir, "reference.sorted.bam"),
        varcall_dir=varcall_dir,
        logger=logger)

    results.current_vcf = os.path.join(varcall_dir, "var.filtered.vcf")
    varcount = 0
    with open(results.current_vcf, "r") as vcf:
        for line in vcf:
            if line.startswith("#"):
                varcount = varcount + 1

    logger.info("Identified {varcount} variants...".format(**locals()))
    return 0


def main(args=None, logger=None):
    if args is None:
        args = get_args()
    check_files_present(args)
    results = make_empty_results_object()
    results.assemblers_list = match_assembler_args(args=args)
    results.organism = "{0}{1}{2}".format(args.genus, " " + args.species, " sp. " + args.strain)
    t0 = time.time()
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
    logger.debug(
        "-----------------------  BEGIN CONFIG ------------------------------")
    for k, v in sorted(vars(config).items()):
        logger.debug("    %s: %s", k, str(v))
    logger.debug(
        "-----------------------  END CONFIG ------------------------------")
    # lets not be hasty
    check_args(args, config)

    setup_tmp_dir(args, output_root=args.outdir, logger=logger)

    # copy paths to raw reads in case we need them for mascura
    args.untrimmed_fastq1 = args.fastq1
    args.untrimmed_fastq2 = args.fastq2

    logger.info("Assessing reads and library type")
    reads_ns  = assess_reads(args=args, config=config,
                             platform=args.platform, logger=logger)
    logger.debug(reads_ns)
    # determine which steps we will be doing
    parse_stages(args, reads_ns, logger, all_stages=all_stages)
    logger.info("Determining if a reference is needed")
    check_ref_required(args=args, lib_type=reads_ns.lib_type)

    logger.info("Preparing config and tools")
    tools = select_tools(args, config=config, reads_ns=reads_ns, logger=logger)

    #-------------------------   QC
    logger.debug("Determining if fastqc will be run")
    if reads_ns.QC:
        if  config.fastqc is not None:
            run_fastqc(reads_ns, args, config, logger=logger)
        else:
            logger.info("Skipping fastqc, as no executable is in PATH")
        finish_if_done(this_stage="q", args=args, all_stages=all_stages,
                       logger=logger, t0=t0, results=results)

    #-------------------------- TRIM
    if reads_ns.TRIM:
        logger.debug("Ensuring appropriate trim length")
        args.trim_length = check_and_set_trim_length(reads_ns, args, logger)
        logger.debug("Determining whether to perform trimming")
        if reads_ns.lib_type != "long":  # and not args.skip_trim:
            logger.info("Trimming reads based on quality")
            quality_trim_reads(args, config, reads_ns, logger)
        finish_if_done(this_stage="t", args=args, all_stages=all_stages,
                       logger=logger, t0=t0, results=results)

    #-------------------------- DOWNSAMPLE
    if reads_ns.DOWNSAMPLE:
        # we need to check other potentially conflicitng args
        ACTUALLY_DOWNSAMPLE = True
        logger.debug("Determining whether to downsample")
        if reads_ns.coverage is None:
            ACTUALLY_DOWNSAMPLE = False
            raise ValueError(
                "Cannot downsample, as coverage could not be " +
                "calcuated! Try rerunning without 'd' in --stages")
        if (reads_ns.coverage < args.downsample):
            ACTUALLY_DOWNSAMPLE = False
            logger.info(str("Skipping downsampling; coverage %d is less " +
                            "than --downsample value of %d"),
                        reads_ns.coverage,
                        args.downsample)
        if not assembler_needs_downsampling(tools):
            ACTUALLY_DOWNSAMPLE = False
            logger.info("Skipping downsampling; one or more assemblers "+
                        "does not need it")
        if ACTUALLY_DOWNSAMPLE:
            logger.info("Downsampling reads to %dx coverage", args.downsample)
            reads_ns.downsampled_coverage = downsample_reads(
                args=args, reads_ns=reads_ns,
                config=config, new_cov=args.downsample, logger=logger)
        finish_if_done(this_stage="d", args=args, all_stages=all_stages,
                       logger=logger, t0=t0, results=results)

    #-------------------------- Report status, pre assembly
    if args.fastq2 and len(args.references) != 0:
        logger.info("Aligning reads to first reference " +
                    "%s for determinging insert size" % args.references[0])
        sorted_bam = align_reads(ref=args.references[0], dirname="align", reads_ns=reads_ns,
                                 config=config, downsample=True, args=args, logger=logger)
        reads_ns.insert_mean, reads_ns.insert_stddev = get_insert_stats(
            bam=sorted_bam, config=config, args=args, logger=logger)

    # print out the run data we have so far
    log_read_and_run_data(reads_ns, args, results, logger)

    #-------------------------- ASSEMBLE (and merge)
    # Righto, now the fun starts. First up, we assemble with 1 or 2 assemblers
    # (or check and read in results from a prior assembly)
    # if not reads_ns.ASSEMBLE:


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
        merged_contigs_path = merge_assemblies(tool=tools.merger, args=args, config=config,
                                               reads_ns=reads_ns, logger=logger)
        results.current_contigs = merged_contigs_path
        results.current_contigs_source = "merged"
        results.current_scaffolds = None

    else:
        results.current_contigs = \
            results.assemblers_results_dict_list[0]['contigs']
        results.current_contigs_source = \
            results.assemblers_results_dict_list[0]['name']
        results.current_scaffolds = \
            results.assemblers_results_dict_list[0]['scaffolds']
        results.current_scaffolds_source = \
            results.assemblers_results_dict_list[0]['name']
    logger.debug("RESULTS:")
    logger.debug(results.__dict__)
    logger.debug("TOOLS:")
    logger.debug(tools.__dict__)
    logger.debug("ARGS:")
    logger.debug(args.__dict__)
    finish_if_done(this_stage="a", args=args, all_stages=all_stages,
                   logger=logger, t0=t0, results=results)
    #-------------------------  Check our reference(s)
    if len(args.references) > 0:
        if len(args.references) == 1:
            results.current_reference = args.references[0]
            results.ID_OK, results.reference_percent_sim = check_id(
                results.current_reference, args, contigs=results.current_contigs,
                results=results, config=config, logger=logger)
        if len(args.references) > 1:
            # if multiple references provided, select the primary, most similar one
            percent_list = []  # [%id, ref]
            for reference in args.references:
                ID_OK, percent = check_id(reference, args, contigs=results.current_contigs,
                                     results=results, config=config, logger=logger)
                if ID_OK:  # ignore any dodgey ones from the list
                    percent_list.append([percent, reference])
                    # doesnt matter which one it is -- already passed
                    results.ID_OK = True
                else:
                    logger.warning(
                        "ignoring reference %s, which only bears ~ %d similarity" %
                        (reference, percent))
                    args.references = [x for x in args.references if x != reference]
            if len(args.references) == 0:  # if we filtered out all the refs
                results.ID_OK = False    #  than we should warn the user
                logger.warning("No appropriate reference given!")
            best_hit = sorted(percent_list)[-1]
            results.current_reference = best_hit[1]
            results.reference_percent_sim = best_hit[0]

    #-------------------------- BREAK at origin
    if reads_ns.BREAK:
    #TODO why do we check scaffolder here?
    #  ensure that our reference is close enough before continuing
        if args.scaffolder and results.current_reference is None:
            logger.error("Cannot break contigs at an origin without a reference")
        # if we have a reference, lets break the origin here, before scaffolding
        if results.ID_OK:
            find_and_split_origin(args=args, config=config, results=results,
                                  tools=None, reads_ns=None, logger=logger)
        finish_if_done(this_stage="b", args=args, all_stages=all_stages,
                       logger=logger, t0=t0, results=results)

    #-------------------------- SCAFFOLD
    # right.  Now our results should contain contigs and/or scaffolds.
    # Lets do some scaffolding
    if reads_ns.SCAFFOLD:
        if args.scaffolder is None:
            raise ValueError("No scaffolder provided! Please either " +
                             "select one using the --scaffolder arg, or " +
                             "remove 's' from the stages argument!")
        # if we didnt break scaffolds, we wont have already checked our reference
        if results.ID_OK is None:
            results.ID_OK, ref_percent = check_id(
                results.current_reference, args, contigs=results.current_contigs,
                results=results, config=config, logger=logger)

        if \
           (results.ID_OK and tools.scaffolder['linkage_evidence'] == "align-genus") or \
             not (results.ID_OK and tools.scaffolder['linkage_evidence'] == "paired-ends"):
            results.old_scaffolds.append(["either nowhere or straight from assembler",
                                         results.current_scaffolds])
            results.current_scaffolds = run_scaffolder(
                args=args, config=config, tools=tools, reads_ns=reads_ns, results=results,
                run_id=1, logger=logger)
        #TODO why is this here? I still dont understand so we will skip for now
        if results.current_scaffolds is not None and False:
            # we may have been given a reference for a long read assembly but no
            # scaffolder is used for these by default
            args.scaffolder = "mauve" if args.scaffolder is None else args.scaffolder
            if not args.skip_split_origin and results.current_reference is not None:
                find_origin(args,config, results,  tools, reads_ns, logger )
            if results.ID_OK:
                order_scaffolds(args, config, results, logger)
        finish_if_done(this_stage="s", args=args, all_stages=all_stages,
                       logger=logger, t0=t0, results=results)

    #-------------------------- FINISHING
    if reads_ns.FINISH:
        REFERENCE_BASED_FINISHER = True
        if results.ID_OK is None:
            results.ID_OK, ref_per = check_id(
                results.current_reference, args, contigs=results.current_contigs,
                results=results, config=config, logger=logger)
        if args.finisher and results.ID_OK:
            run_finisher(args, config, reads_ns, tools, results, logger=logger)
        finish_if_done(this_stage="f", args=args, all_stages=all_stages,
                       logger=logger, t0=t0, results=results)

    #-------------------------- ANNOTATE WITH PROKKA

    # sequence stable from this point, only annotations altered
    #  it looks this this is attempting to be multiprocessed.
    #  for now, lets just do this in serial, as this cant save much time;
    #  id rather prokka have all the cores it can
    # simul_cmds = []
    # pool = multiprocessing.Pool(processes=split_cores)
    # logger.debug("running the following commands:")
    # logger.debug("\n".join(simul_cmds))
    # spades_results = [
    #     pool.apply_async(subprocess.run,
    #                      (cmd,),
    #               { "shell": sys.platform != "win32",
    #                "stdout": subprocess.PIPE,
    #                "stderr": subprocess.PIPE,
    #                "check": True})
    #     for cmd in simul_cmds]
    # pool.close()
    # pool.join()
    # logger.debug(spades_results)
    # logger.info("Sum of return codes (should be 0):")
    # spades_results_sum = sum([r.get() for r in spades_results])

    ##amosvalidate fails if we don't have mate-pairs
    # print(results)
    if reads_ns.GENECALL:
        if reads_ns.paired and args.mode == 'draft':
            seq_file = results.current_contigs
            if results.current_scaffolds is not None:
                seq_file = results.current_scaffolds
            sorted_bam = align_reads(
                ref=results.current_reference,
                dirname="align_pregenecall", reads_ns=reads_ns,
                config=config, downsample=True, args=args, logger=logger)
            # amosvalidate(args, results, config, reads_ns, logger)
            # # get_contig_to_iid_mapping($tmpdir);
            # amosvalidate_results = summarise_amosvalidate(args.tmp_dir)
        try:
            evidence = tools.scaffolder['linkage_evidence']
        except:
            evidence = "paired-ends"
        gaps = build_agp(args, results, reads_ns,
                         evidence=evidence,
                         logger=logger)
        prokka_dir = run_prokka(config, args, results, logger)
        make_embl_from_gbk(gbk=os.path.join(prokka_dir, "prokka.gbk"),
                           output_file=os.path.join(prokka_dir, "prokka.embl"))
        amosvalidate_results = None
        if reads_ns.VARCALL:
            run_varcaller(args, results, reads_ns, config, tools, logger)
        # merge_annotations(tmpdir, amosvalidate_results, gaps, genus, species, strain )
        #     #  This kept throwing an error about Bio::SeqIO
        run_cgview(results=results, args=args, config=config, logger=logger)

        build_comparisons(args=args, config=config, results=results, logger=logger)
        make_gene_report(results=results, logger=logger)
        finish_if_done(this_stage="p", args=args, all_stages=all_stages,
                       logger=logger, t0=t0, results=results)
        # logger.info("\n\nFinal Assembly Statistics...");
        # report_config = get_contig_stats(results.current_contigs)
        # report_scaffold = get_contig_stats(results.current_scaffolds)
        # logger.info ("\n\nFinal Contig Statistics:\n"
        #              "=========================\n" +
        #              tabulate.tabulate(report_config) + "\n")
        # logger.info ("\n\nFinal Scaffold Statistics:\n"
        #              "=========================\n" +
        #              tabulate.tabulate(report_scaffold) + "\n")
        # logger.info("Done building your Bug!")
        # logger.info("Elapsed time: %.2fm" % ((time.time() - t0) / 60))




















def merge_annotations():
    """
    updates annotated generated embl file with amosvalidate results
    and gap locations

    required parameters: $ (tmpdir)
                     $ (hashref of amosvalidate results, keyed on contigid)
                     $ (hashref of scaffold gaps)
                     $ (genus)
                     $ (species)
                     $ (strain)

    returns            : $ (none)
    """

    logger.info("Merging annotations");
    merge_dir = os.path.join(args.tmp_dir, "annotation_merge")
    os.makedirs(merge_dir)
    filename = None
    new_embl = os.path.join(merge_dir, os.path.basename(results.current_embl))
    with open(results.current_embl, "r") as IO, open(new_embl, "w") as outIO:

        amos_colours = {
            'CE_STRETCH'     :  '0 128 128',
            'CE_COMPRESS'    :  '0 128 128',
            'HIGH_SNP'       :  '128 128 0',
            'HIGH_READ_CVG'  :  '255 0 0',
            'HIGH_KMER'      :  '255 0 0',
            'KMER_COV'       :  '255 0 0',
            'HIGH_OUTIE_CVG' :  '255 0 0',
            'HIGH_NORMAL_CVG':  '255 0 0',
            'LOW_GOOD_CVG'   :  '0 0 255',
        }

        amos_notes = {
            'CE_STRETCH'      : 'Stretched mate-pairs: Possible repeat copy number expansion or other insertion',
            'CE_COMPRESS'     :  'Compressed mate-pairs; Possible collapsed repeat',
            'HIGH_SNP'        :  'High SNP frequency',
            'HIGH_READ_CVG'   :'High read coverage; Possible collapsed repeat',
            'HIGH_KMER'       :  'High frequency of normalized kmers: Possible collapsed repeat',
            'KMER_COV'        :'High frequency of normalized kmers: Possible collapsed repeat',
            'LOW_GOOD_CVG'    : 'Low coverage',
            'HIGH_NORMAL_CVG' : 'High coverage',
            'HIGH_OUTIE_CVG'  : 'High outie coverage',
                     }
        for embl_record in SeqIO.parse(IO, "embl"):
            orig_id = embl_record.id
            embl_record.id = orig_id
            embl_record.accession_number = orig_id
            embl_record.division = 'PRO'
            embl_record.molecule = 'genomic DNA'
            embl_record.is_circular = 1
            embl_record.date = arrow.now().format("DD-MMM-YYYY").upper()
            embl_record.description = "{0} {1} {2} genome scaffold".format(
                args.genus, args.species, args.strain)

            embl_record.comment = "Assembled using BugBuilder from http://github.com/jamesabbott/BugBuilder"
            # remove source entries from feature table
            features = [x for x in  embl_record.features]
            # retrieve amosvalidate results for this contig and sort by start co-ordinate...
"""
            if ($amosvalidate_results) {
                    my @amos_features =
                    map  { $_->[0] }
                    sort { $a->[1] <=> $b->[1] }
                    map  { [ $_, $_->{'start'} ] } @{ $amosvalidate_results->{$orig_id} };

                foreach my $feature (@amos_features) {
                    if ( $feature->{'start'} != $feature->{'end'} ) {
                            my $colour = $amos_colours{ $feature->{'type'} };
                            my $note   = $amos_notes{ $feature->{'type'} };
                            my $feature =
                            new Bio::SeqFeature::Generic(
                                -start       => $feature->{'start'} + 1,
                                -end         => $feature->{'end'},
                                -primary_tag => 'misc_feature',
                                -tag         => {
                                    'note'   => $note,
                                    'colour' => $colour,
                                }
                            );
                            push @features, $feature;
####
        if ($gaps) {
            foreach my $scaffold ( keys(%$gaps) ) {
                if ( $orig_id eq $scaffold ) {
                    my $scaffold_gaps = $gaps->{$scaffold};
                    foreach my $gap (@$scaffold_gaps) {
                        my ( $gap_start, $gap_end ) = split( /-/, $gap );
                        my $est_length;
                        if ( $gap_end - $gap_start == 100 ) {
                            $est_length = 'unknown';
                        }
                        else {
                            $est_length = $gap_end - $gap_start;
                        }
                        my $feature =
                          new Bio::SeqFeature::Generic(
                                          -start       => $gap_start,
                                          -end         => $gap_end,
                                          -primary_tag => 'assembly_gap',
                                          -tag => { 'estimated_length' => $est_length, 'gap_type' => 'within_scaffold' }
                          );
                        push( @features, $feature );
                    }
                }
            }
        }

        my $source =
          new Bio::SeqFeature::Generic(
                                        -start       => 1,
                                        -end         => $embl_record->length(),
                                        -primary_tag => 'source',
                                        -tag         => {
                                                  'organism' => "$genus $species $strain",
                                                  'strain'   => $strain
                                                }
                                      );

        $embl_record->add_SeqFeature($source);
        @features = map { $_->[0] }
          sort { $a->[1] <=> $b->[1] }
          map { [ $_, $_->start() ] } @features;

        foreach my $feat (@features) {
            if ( $feat->primary_tag() ne 'source' ) {
                $embl_record->add_SeqFeature($feat);
            }
        }

        $outIO->write_seq($embl_record);
    }

    chdir $tmpdir      or die "Error chdiring to $tmpdir: $!";
    unlink "$filename" or die "Error unlinking $filename: $!";
    symlink( "annotation_merge/$filename", "$filename" )
      or die "Error creating $filename symlink: $!";

}
"""
