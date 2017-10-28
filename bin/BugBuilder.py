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
"""
=pod

=head1 NAME

BugBuilder

=head1 SYNOPSIS

  BugBuilder --fastq1 fastq_read1 [--fastq2 fastq_read2] --platform platform [--reference reference_genome.fa] [--help] [--manual]

=head1 DESCRIPTION

  Automated pipeline for assembly of draft quality bacterial genomes with reference guided scaffolding
  and annotation.

  Please see accompanying userguide for full documentation

=head1 REQUIRED ARGUMENTS

=over 4

=item B<platform>: Sequencing platform  used i.e. illumina, 454, iontorrent

=back

=head1 OPTIONAL ARGUMENTS

=over 4

=item B<fastq1>: Path to first read of paired library, or fragment library

=item B<fastq2>: Path to second read of paired library

=item B<longfastq>: Path to fastq file from long-read sequencer

=item B<prefix>: Prefix to use for output file naming

=item B<reference>: Path to fasta formatted reference sequence

=item B<assembler>: Assembler(s) to run - may be specified twice, in which case the two assemblers
will be run in parallel and the results merged using minimus. If no assembler is specified,
BugBuilder will try to select an appropriate assembler automatically

=item B<assembler-args>: Any additional arguments to pass to the assembler. Default values are set
in the 'default_args' attribute of the configuration file. If running multiple assemblers,
assembler_args should be specified twice, once for each assemler, in the same order than the
assemblers are specified.

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

use Archive::Zip qw(:ERROR_CODES :CONSTANTS);
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::Search::Tiling::MapTiling;
use Bio::DB::Fasta;
use Cwd;
use DateTime;
use File::Temp qw(tempdir);
use File::Path qw(make_path rmtree);
use File::Copy;
use File::Copy::Recursive qw(dircopy);
use File::Basename;
use File::Find::Rule;
use File::Tee qw(tee);
use File::Which qw(which);
use FindBin;
use Getopt::Long;
use Parallel::ForkManager;
use PerlIO::gzip;
use Pod::Usage;
use Statistics::Basic qw(:all);
use Term::ANSIColor qw(:constants);
use Text::ASCIITable;
use YAML::XS qw(LoadFile);
"""

import argparse
import os
import sys
import arrow
import logging

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
    
def parse_config():
    pass

def parse_available_platforms():
    plotform_list = ["short_illumina", "long_illumina", "de_fere"]
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


def get_args():  # pragma: no cover
    """
    """
    parser = argparse.ArgumentParser(prog="ribo scan",
        description="Build a bug!",
        add_help=False)  # to allow for custom help
    parser.add_argument("contigs", action="store",
                        help="either a (multi)fasta or a directory " +
                        "containing one or more chromosomal " +
                        "sequences in fasta format")

    # taking a hint from http://stackoverflow.com/questions/24180527
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("--platform", dest='platform', action="store",
                               help="sequencing platform",
                               choices=parse_available_platforms,
                               type=str, required=True)
    requiredNamed.add_argument("-o", "--output", dest='output', action="store",
                               help="output directory; " +
                               "default: %(default)s", default=os.getcwd(),
                               type=str, required=True)

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--fastq1", dest='fastq1', action="store",
                               help="forward reads. or single-end library; ",
                               type=str, required="--long-fastq" not in sys.argv)
    optional.add_argument("--fastq2", dest='fastq2', action="store",
                               help="reverse reads; ",
                               type=str)
    optional.add_argument("--de_fere_contigs", dest='de_fere_contigs',
                               action="store",
                               help="contigs to be used in de fere novo assembly!"
                               type=str)
    optional.add_argument("--long_fastq", dest='long_fastq', action="store",
                               help="long_fastq reads, such as from Pacbio",
                               type=str, required="--fastq1" not in sys.argv)
    optional.add_argument("--reference", dest='reference',
                          action="store",
                          help="reference fasta",
                          type=str)
    optional.add_argument("--prefix", dest='prefix',
                          action="store",
                          help="prefix for output fasta",
                          type=str)
    optional.add_argument("--assembler", dest='assembler', action="store",
                          help="assembler to use",
                          choices=parse_available_assemblers,
                          nargs="+",
                          type=str)
    optional.add_argument("--assembler-args", dest='assembler-args',
                          action="store",
                          help="args to pass to the assembler, in single quotes",
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
    optional.add_argument("--merge-method", dest='merge-method', action="store",
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
    optional.add_argument("--insert-size", dest='insert-size', action="store",
                          help="insert size (bp)",
                          type=int)
    optional.add_argument("--insert-stddev", dest='insert-stddev', action="store",
                          help="insert size standard deviation ",
                          type=int)
    optional.add_argument("--genome-size", dest='genome-size', action="store",
                          help="size of the genome ",
                          type=int)
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
    optional.add_argument("--skip-fastqc", dest='skip-fastqc', action="store_true",
                          help="size of the genome ", default=False)
    optional.add_argument("--skip-trim", dest='skip-time', action="store_true",
                          help="Quality threshold for trimming reads ", default=False)
    optional.add_argument("--trim-qv", dest='trim-qv', action="store",
                          help="quality threshold to trim  ", default=20,
                          type=int)
    optional.add_argument("--trim-length", dest='trim-length', action="store",
                          help="min read klength to retain", default=50,
                          type=int)
    optional.add_argument("--skip-split-origin", dest='skip-split-origin',
                          action="store_true",
                          help="split at origin ", default=False)
    optional.add_argument("--keepall", dest='keepall',
                          action="store_true",
                          help="keep all intermediate files ", default=False)
    optional.add_argument("--threads", dest='threads', action="store",
                          help="threads  ",
                          type=int)
    optional.add_argument("--out-dir", dest='out-dir', action="store",
                          help="dir for results",
                          type=str)
    optional.add_argument("--tmp-dir", dest='tmp-dir', action="store",
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
"""#  setup_tempdir creates a temporary directory and copies
#  over the relevent files
#
#  requried params: $ (path to fastq1)
#                   $ (path to fastq2)
#                   $ (path to longread fastq file)
#                   $ (path to fastafile)
#  optional params: $ (platform)
#                   $ (preset value for tempdir)
#
#
#  returns        : $ (path to tempdir)
#
"""
    if args.tmpdir is None:
        scratch_dir = os.path.join(output_root, "temp_dir")
    else:
        scratch_dir = args.tmpdir;

    print "Created working directory: ${tmpdir}\n";
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


def get_read_len_from_fastq(fastq1, logger=None):
    """from LP; return total read length and count in fastq1 file from all reads
    """
    lengths = []
    if os.path.splitext(fastq1)[-1] in ['.gz', '.gzip']:
        open_fun = gzip.open
    else:
        open_fun = open
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


def  assess_reads(args, tempdir, logger=None):
"""
Determine some characteristics of the provided reads. Need to determine
the mean read length, wheterh they are variable length or not, and the
quality encoding.

We'll also look at the predicted coverage if a reference genome is available

required params: $ (tmpdir)

optional params: $ (platform)

returns: $ (assembler type)
         $ (quality encoding)
         $ (mean length)
         $ (stddev of read length)
         $ (projected genome coverage)
         $ (number of reads)
"""
    logger.info("Checking reads...");

    my ( @mean, @stddev, $long_mean, $long_stddev, $tot_length, $long_tot_length );
    # get the lengths of all the short reads
    short_read_lengths, short_tot_lengths, short_counts = [], [], []
    for read in [args.fastq1, args.fastq2]:
        if not os.path.exists(read):
            next
        lengths  = get_read_len_from_fastq(read, logger=logger)
        short_lengths.append(lengths)
        short_tot_length.append(sum(lengths))
        short_count.append(len(lengths))
    # ensure paired reads have same numebr of reads
    if len(short_counts) > 1:
        if short_counts[0] != short_counts[1]:
            raise ValueError("Paired reads must have same number of reads!")
        
    if os.path.exists(args.long_fastq):
        long_lengths = get_read_lens_from_fastq(args.long_fastq, logger=logger)
    else:
        long_lengths, long_tot_length, long_count = [], 0, 0


    lengths = (tot_length + tot_long_length)
    mean = float((tot_length + tot_long_length) / (count + long_count)))
    
        my ( @lengths, $reads );
        my $target = readlink("$tmpdir/$read");
        if ( $target =~ /.gz$/ ) {
                open FASTQ, "<:gzip", "$tmpdir/$read" or die "Error opening $tmpdir/$read: $!";
        }
        else {
            open FASTQ, "$tmpdir/$read" or die "Error opening $tmpdir/$read: $!";
        }
        while ( my $line = <FASTQ> ) {
            $line = <FASTQ>;
            if ($line) {    #suppress undefined warnings for 0 length reads
                push @lengths, length($line);
		#<#
                $tot_length += length($line) if ( $read ne 'long.fastq' );
                    $long_tot_length += length($line) if ( $read eq 'long.fastq' );
                    $reads++;
                }

                #discard the spare lines....
                $line = <FASTQ>;
                $line = <FASTQ>;
            }
            close FASTQ;

            # remove any unusually long reads which mess up the statistics (yes, 454...were looking at you....);
            my $mean = mean(@lengths);
            my @tmp;
            foreach (@lengths) {
                push @tmp, $_ if ( $_ < 3 * $mean );
            }
            @lengths = @tmp;
	    #>#
            if ( $read eq 'long.fasta' ) {
                $long_mean   = sprintf( "%d", mean(@lengths) );
                $long_stddev = sprintf( "%d", stddev(@lengths) );
            }
            else {
                push @mean,   mean(@lengths);
                push @stddev, stddev(@lengths);

            }
        }
    }

    #determin mean read length and stddev....
    my ( $mean, $stddev );

    #  unless ( defined($long_mean) ) {
    if ( $#mean == 1 ) {
        $mean   = sprintf( "%d", ( ( $mean[0] + $mean[1] ) / 2 ) );
        $stddev = sprintf( "%d", ( ( $stddev[0] + $stddev[1] ) / 2 ) );
    }
    else {
        $mean   = sprintf( "%d", $mean[0] );
        $stddev = sprintf( "%d", $stddev[0] );

        #      }
    }

    my $type;

  TYPE: foreach my $category ( @{ $config->{'assembler_categories'} } ) {
        if ($platform) {
            next TYPE unless ( grep ( /^$platform$/i, @{ $category->{'platforms'} } ) );
            $type = $category->{'name'};
        }

        if ($mean) {
            if ( ( $mean > $category->{'min_length'} ) && ( $mean <= $category->{'max_length'} ) ) {
                $type = $category->{'name'};
                last TYPE;
            }
        }
        elsif ($long_mean) {

            if ( ( $long_mean > $category->{'min_length'} ) && ( $long_mean <= $category->{'max_length'} ) ) {
                $type = $category->{'name'};
                last TYPE;
            }
        }
    }

    my ( $coverage, $long_coverage );

    if ( -e "$tmpdir/reference.fasta" ) {
        my $inIO = Bio::SeqIO->new( -file => "$tmpdir/reference.fasta", -format => 'fasta' );
        while ( my $ref = $inIO->next_seq() ) {
            $genome_size = $genome_size + $ref->length();
        }
    }
    $coverage      = sprintf( "%d", $tot_length / $genome_size )      if ( $genome_size && $tot_length );
    $long_coverage = sprintf( "%d", $long_tot_length / $genome_size ) if ( $genome_size && $long_tot_length );

    my $encoding = id_fastq_encoding($tmpdir);

    return ( $type, $encoding, $mean, $stddev, $long_mean, $long_stddev, $coverage, $long_coverage, $tot_length );

}


def main(args):
    if args is None:
        args = get_args()
    assembler_list = match_assembler_args(
        assembler=args.assembler,
        assembler_args=args.assembler_args))
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
    PAIRED = False if not args.fastq2

    logger.debug( "Logging output to %s", log)

    logger.info("Welcome to BugBuilder")
    logger.info("\nPreparing to build your bug...\n\n")

    tmpdir = setup_tempdir(args, output_root)

    reference = os.path.basename(args.reference)

    #  need to know a bit about the provided reads so we can act appropriately
    my (
	$type,               $encoding,              $mean_read_length,
	$read_length_stddev, $mean_long_read_length, $long_read_length_stddev,
	$coverage,           $long_read_coverage,    $read_bases
	) = assess_reads( $tmpdir, $platform, $genome_size );

    die "\nPlease supply a genome size or reference sequence when running a long-read assembly"
	if ( $type eq 'long' && ( !$genome_size && !$reference ) );

    my ( $sel_assemblers, $sel_scaffolder, $scaffold_type, $sel_merge_method, $sel_finisher, $sel_varcall ) =
	select_tools( $config, \@assemblers, $scaffolder, $merge_method,     $finisher, $varcall,
		      $paired, $reference,   $type,       $mean_read_length, $insert_size );
    @assemblers   = @$sel_assemblers;
    $scaffolder   = $sel_scaffolder;
    $merge_method = $sel_merge_method;
    $finisher     = $sel_finisher;
    $varcall      = $sel_varcall;

    # pull downsampling attribute from assembler configuration
    my $assembler_conf_downsample = 0;
    my $assembler_conf            = $config->{'assemblers'};
    foreach my $conf_assembler (@$assembler_conf) {
	if ( $conf_assembler->{'name'} eq $assemblers[0] ) {
	    $assembler_conf_downsample = $conf_assembler->{'downsample_reads'};
	}
    }

    run_fastqc( $tmpdir, $type ) if ($run_fastqc);

    if ( $mean_read_length && ( $mean_read_length < 50 || $mean_read_length < $trim_length ) ) {
	$trim_length = 25;
	print "\ntrim-length reset to $trim_length due to mean read length ($mean_read_length)\n";
    }

    quality_trim_reads( $tmpdir, $encoding, $trim_qv, $trim_length )
	unless ( $type eq 'long' || !$trim_reads );

    my $downsampled_coverage;

    if ( ( $coverage && $coverage > 100 ) && ( $assembler_conf_downsample || $downsample ) ) {
	$downsampled_coverage =
	    downsample_reads( $tmpdir, $config, $downsample, $coverage, $read_bases, $mean_read_length );
    }

    if ( $fastq2 && $reference ) {
	## final arg to align_reads is flag to indicate downsamping required...
	align_reads( $tmpdir, $reference, $mean_read_length, 1 );
	( $insert_size, $stddev ) = get_insert_stats( "$tmpdir", $reference );
    }
    my $paired_str;
    ($paired) ? ( $paired_str = "Paired" ) : ( $paired_str = "Fragment" );
    my $tb = Text::ASCIITable->new( { 'headingText' => 'Read Information' } );
    $tb->setCols( 'Aaaaaaaaargh', 'Phwwwwwwweb' );
    $tb->setOptions( { 'hide_HeadRow' => 1 } );
    $tb->alignCol( 'Aaaaaaaaargh', 'left' );
    $tb->alignCol( 'Phwwwwwwweb',  'left' );
    $tb->addRow( "Mean Read Length",                  $mean_read_length )          if ($mean_read_length);
    $tb->addRow( "Read Length Standard Deviation",    $read_length_stddev )        if ($read_length_stddev);
    $tb->addRow( "Insert size",                       $insert_size )               if ($insert_size);
    $tb->addRow( "Insert size Standard Deviation",    $stddev )                    if ($stddev);
    $tb->addRow( "Mean Long Read Length",             $mean_long_read_length )     if ($mean_long_read_length);
    $tb->addRow( "Mean Long Read Standard Deviation", $long_read_length_stddev )   if ($long_read_length_stddev);
    $tb->addRow( "Library type",                      $paired_str );
    $tb->addRow( "Platform",                          $platform )                  if ($platform);
    $tb->addRow( "Quality Encoding",                  $encoding );
    $tb->addRow( "Projected Coverage",                "${coverage}x" )             if ($coverage);
    $tb->addRow( "Projected Coverage (Downsampled)",  "${downsampled_coverage}x" ) if ($downsampled_coverage);
    $tb->addRow( "Projected Long Read Coverage",      "${long_read_coverage}x" )   if ($long_read_coverage);
    print "$tb\n";

    $tb = Text::ASCIITable->new( { 'headingText' => 'Assembly Information' } );
    $tb->setCols( 'Aaaaaaaaargh', 'Phwwwwwwweb' );
    $tb->setOptions( { 'hide_HeadRow' => 1 } );
    $tb->alignCol( 'Aaaaaaaaargh', 'left' );
    $tb->alignCol( 'Phwwwwwwweb',  'left' );
    $tb->addRow( "Assembly category",        $type );
    $tb->addRow( "Selected assemblers",      join( ",", @assemblers ) );
    $tb->addRow( "Selected assembly merger", $merge_method ) if ($merge_method);
    $tb->addRow( "Selected scaffolder",      $scaffolder ) if ($scaffolder);
    $tb->addRow( "Selected finisher",        $finisher ) if ($finisher);
    $tb->addRow( "Selected variant caller",  $varcall ) if ($varcall);
    $tb->addRow( "Trim QV",                  $trim_qv ) if ($trim_reads);
    $tb->addRow( "Trim length",              $trim_length ) if ($trim_reads);
    $tb->addRow( "Split Origin",             $split_origin );
    print "$tb\n";

    my $pm = new Parallel::ForkManager(2);
    my $id_ok = 1;
    if ( !$already_assembled ){
	print "Running in $mode mode\n";
	my $ret;
	my $link = 0;

	#if we are running a single assembler we need to symlink its contig output
	# into $tmpdir - pass a $link flag to run_assembler
	$link = 1 if ( $#assemblers == 0 );
	for ( my $i = 0 ; $i <= $#assemblers ; $i++ ) {

	    my $assembler = $assemblers[$i];
	    my $args      = $assembler_args[$i];

	    run_assembler(
		$tmpdir,      $assembler, $reference,   $args,     $link,
		$type,        $encoding,  $genome_size, $platform, $mean_read_length,
		$insert_size, $stddev,    $threads
		);
	}
	merge_assemblies( $tmpdir, \@assemblers, $merge_method ) if ( $#assemblers == 1 );
	if ( $scaffolder && $reference ) {
	    $id_ok = check_id($tmpdir);
	}
    } else { # end if already assembled
    	print "using assemblies from another directory\n";
    	$id_ok = 1;  # it better be, that is
	$tmpdir = $scratchdir;
	my $oridir = $tmpdir . "/origin";
	rmtree "$tmpdir/origin" or print "Unable to remove origin dir: $!\n";
	rmtree "$tmpdir/agp" or print "Unable to remove agp dir: $!\n";
	unlink "$tmpdir/scaffolds.agp" or print "unable to delete scaffolds.agp\n";
	unlink "$tmpdir/contigs.embl" or print "unable to delete contigs.embl\n";
	unlink "$tmpdir/scaffolds.embl" or print "unable to delete scaffolds.embl\n";
	# unlink "$tmpdir/scaffolds.fasta" or print "unable to delete scaffolds.fasta\n";
	# unlink "$tmpdir/contigs.fasta" or print "unable to delete scaffolds.contigs\n";
	rmtree "$tmpdir/annotation_merge" or print "unable to delete merge\n";
	rmtree "$tmpdir/cgview" or print "unable to delete cgview\n";
	rmtree "$tmpdir/comparisons" or print "unable to delete comparisons\n";
	unlink "$tmpdir/contigs_cgview.png" or print "unable to delete contigs_cgview.png\n";
	unlink "$tmpdir/scaffolds_cgview.png" or print "unable to delete scaffolds_cgview.png\n";
	unlink glob "$tmpdir/comparison_vs*" or print "unable to delete comparisons";
    	print "idok: $id_ok\n";
    }
    run_scaffolder( $tmpdir, "$tmpdir/reference_parsed_ids.fasta",
                    $scaffolder, $scaffolder_args, $insert_size, $stddev, 1, "$tmpdir/contigs.fasta", $mean_read_length,
                    $threads )
      if ( $scaffolder && ( $id_ok && $scaffold_type eq 'align_genus' )
           || ( $id_ok && $scaffold_type eq 'paired-ends' )
           # why is that id_ok negated?
           || ( !$id_ok && $scaffold_type eq 'paired-ends' ) );

    if ( -e "$tmpdir/scaffolds.fasta" && $reference ) {

        #we may have been given a reference for a long read assembly but no scaffolder is used for these by default
        $scaffolder = "mauve" if ( !$scaffolder );

        find_origin( $tmpdir, $scaffolder, $scaffolder_args, "$tmpdir/reference_parsed_ids.fasta",
                     $insert_size, $stddev, $mean_read_length, $threads )
          if ($split_origin);
        order_scaffolds( $tmpdir, basename($reference) ) if ($id_ok);

        finish_assembly( $tmpdir, $finisher, $insert_size, $stddev, $encoding, $threads ) if ( $finisher && $id_ok );

    }

    my $gaps;

    if ( -e "$tmpdir/scaffolds.fasta" ) {
        $gaps = build_agp( $tmpdir, $organism, $mode, $scaffold_type );
    }

    # sequence stable from this point, only annotations altered
    for my $i (qw(1 2)) {
	print "pm: $i \n";
        $pm->start and next();
        if ( $i == 1 ) {

            ##amosvalidate fails if we don't have mate-pairs
            if ( -e "$tmpdir/read2.fastq" && ( $mode eq 'draft' ) ) {
                my $seq_file;
                ( -e "$tmpdir/scaffolds.fasta" ) ? ( $seq_file = 'scaffolds.fasta' ) : ( $seq_file = 'contigs.fasta' );

                align_reads( $tmpdir, $seq_file, $mean_read_length );
                amosvalidate( $tmpdir, $insert_size, $stddev );
            }
        }
        elsif ( $i == 2 ) {

            run_prokka( $tmpdir, $genus, $species, $strain, $locustag, $centre );
        }

        $pm->finish();
    }

    $pm->wait_all_children();

    my $amosvalidate_results;
    if ( -e "$tmpdir/read2.fastq" && ( $mode eq 'draft' ) ) {
        get_contig_to_iid_mapping($tmpdir);
        $amosvalidate_results = summarise_amosvalidate($tmpdir);
    }

    run_varcaller( $tmpdir, $varcall, $threads, $mean_read_length ) if ($varcall);
    merge_annotations( $tmpdir, $amosvalidate_results, $gaps, $genus, $species, $strain );
    #  This kept throwing an error about Bio::SeqIO 
    # run_cgview($tmpdir);

    build_comparisons( $tmpdir, basename($reference), $organism ) if ($reference);

    message("Final Assembly Statistics...");
    get_contig_stats( "$tmpdir/contigs.fasta", 'contigs' );
    get_contig_stats( "$tmpdir/scaffolds.fasta", 'scaffolds' ) if ( -e "$tmpdir/scaffolds.fasta" );
    my $emblfile;
    ( -e "$tmpdir/scaffolds.embl" ) ? ( $emblfile = "$tmpdir/scaffolds.embl" ) : ( $emblfile = "$tmpdir/contigs.embl" );

    my $io = Bio::SeqIO->new( -format => 'embl', -file => "$emblfile" );
    my ( $cds, $tRNA, $rRNA );
    while ( my $contig = $io->next_seq() ) {
        foreach my $feature ( $contig->get_SeqFeatures() ) {
            $cds++  if ( $feature->primary_tag eq 'CDS' );
            $tRNA++ if ( $feature->primary_tag eq 'tRNA' );
            $rRNA++ if ( $feature->primary_tag eq 'rRNA' );
        }
    }

    $tb = Text::ASCIITable->new();
    $tb->setCols( "Feature Type", "Number" );
    $tb->addRow( "CDS",  $cds );
    $tb->addRow( "tRNA", $tRNA );
    $tb->addRow( "rRNA", $rRNA );
    print "\nAnnotated features\n==================\n\n";
    print $tb, "\n";

    return_results( $tmpdir, $out_dir, $prefix, $keepall, $mode );
    chdir $orig_dir or warn "Failed to chdir:$ !";

    message("All done...");

}

######################################################################
#
# set_paths
#
# Setups the PATH, PYTHONPATH and PERL5LIB environmental variables
# according to the specifications in the config...
#
# required arguments: $ (config hashref)
#
# returns:            0
#
######################################################################

sub set_paths {

    my $config = shift;

    # Although we execute tools using fully qualified paths,
    # some of these expect certain things to be on path...
    $ENV{'PATH'} = $ENV{'PATH'} . ':' . $FindBin::Bin;
    $ENV{'PATH'} = $ENV{'PATH'} . ":" . "$FindBin::Bin/../packages/bin";
    $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'R_dir'} . '/bin'
      if ( $config->{'R_dir'} );
    $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'blast_dir'}
      if ( $config->{'blast_dir'} );
    $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'mummer_dir'}
      if ( $config->{'mummer_dir'} );
    $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'ncbi_utils_dir'}
      if ( $config->{'ncbi_utils_dir'} );
    $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'tbl2asn_dir'}
      if ( $config->{'tbl2asn_dir'} );
    $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'asn2gb_dir'}
      if ( $config->{'asn2gb_dir'} );
    $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'aragorn_dir'} . '/bin'
      if ( $config->{'aragorn_dir'} );
    $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'prodigal_dir'}
      if ( $config->{'prodigal_dir'} );
    $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'hmmer3_dir'} . '/bin'
      if ( $config->{'hmmer3_dir'} );

    #$ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'rnammer_dir'}
    #  if ( $config->{'rnammer_dir'} );
    $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'infernal_dir'} . '/bin'
      if ( $config->{'infernal_dir'} );
    $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'prokka_dir'}
      if ( $config->{'prokka_dir'} );
    $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'barrnap_dir'}
      if ( $config->{'barrnap_dir'} );
    $ENV{'PATH'} = $ENV{'PATH'} . ":" . $config->{'mauve_dir'} . '/linux-x64'
      if ( $config->{'mauve_dir'} );
    $ENV{'PATH'} = $ENV{'PATH'} . ':' . $config->{'abyss_sealer_dir'}
      if ( $config->{'abyss_sealer_dir'} );

    if ( $config->{'python_lib_path'} ) {
        if ( defined( $ENV{'PYTHONPATH'} ) ) {
            $ENV{'PYTHONPATH'} = "$ENV{'PYTHONPATH'}:" . $config->{'python_lib_path'};
        }
        else {
            $ENV{'PYTHONPATH'} = $config->{'python_lib_path'};
        }
    }

    if ( $config->{'perl_lib_path'} ) {
        if ( $ENV{'PERL5LIB'} ) {
            $ENV{'PERL5LIB'} = "$ENV{'PERL5LIB'}:" . $config->{'perl_lib_path'};
        }
        else {
            $ENV{'PERL5LIB'} = $config->{'perl_lib_path'};
        }
    }
    return;
}

######################################################################
#
# select_tools
# options, while checking for validity of selections...
#
# required params: $ ($config hash ref)
#                  $ (arrayref of requested assemblers)
#                  $ (name of requested scaffolder)
#                  $ (name of requested merge tool)
#                  $ (name of requested finishing tool)
#                  $ (name of requested variant caller)
#                  $ (paired - flag to indicate we are using paired reads
#                  $ (reference - indicates a reference has been provided)
#                  $ (type - category of sequence i.e. short_illumina)
#                  $ (mean read length)
#
# returns          $ (arrayref of assemblers to use)
#                  $ (name of scaffolder to use)
#
######################################################################

sub select_tools {

    my $config      = shift;
    my $assemblers  = shift;
    my $scaffolder  = shift;
    my $merger      = shift;
    my $finisher    = shift;
    my $varcall     = shift;
    my $paired      = shift;
    my $reference   = shift;
    my $type        = shift;
    my $read_length = shift;
    my $insert_size = shift;

    my @assemblers = @$assemblers;
    die "\nA maximum of 2 assemblers can be requested\n" if ( $#assemblers > 1 );

    # sanity check assemblers where these have been manually requested...
    if ( $#assemblers > -1 ) {
        my @available_assemblers = map { $_->{'name'} } @{ $config->{'assemblers'} };

        # first off, are the requested assemblers available
        foreach my $assembler (@assemblers) {
            unless ( grep ( /^$assembler$/i, @available_assemblers ) ) {
                die "\nError: $assembler is not a valid assembler\n\n" . "The following assemblers are configured: ",
                  join( ", ", @available_assemblers ), "\n\n";
            }
        }

        # then check the requested assemblers are appropriate for the sequence characteristics...
        foreach my $assembler (@assemblers) {
            my $assemblers = $config->{'assemblers'};
            foreach my $conf_assembler (@$assemblers) {
                my $name = lc( $conf_assembler->{'name'} );
                if ( $name eq $assembler ) {
                    if ( !$paired & !$conf_assembler->{'command_se'} ) {
                        die "\nError: $assembler requires paired reads, but you only specified one fastq file...";
                    }
                    if ( $conf_assembler->{'min_length'} && $read_length < $conf_assembler->{'min_length'} ) {
                        die "\nError: $assembler does not support reads less than "
                          . $conf_assembler->{'min_length'} . " bp";
                    }
                    if ( $conf_assembler->{'max_length'} && $read_length > $conf_assembler->{'max_length'} ) {
                        die "\nError: $assembler does not support reads longer than "
                          . $conf_assembler->{'max_length'} . " bp";
                    }
                    if ( $conf_assembler->{'insert_size_required'} && ( !$reference && !$insert_size ) ) {
                        die "\nError: $assembler requires the library insert size and stddev to be provided."
                          . "\nPlease add the --insert-size and --insert-stddev parameters, or provide a reference sequence";
                    }
                }
            }
        }
    }
    else {

        # If no assembler is requested, we need to select one based on the 'assembly_type' from the
        # configuration file....
        foreach my $category ( @{ $config->{'assembler_categories'} } ) {
            if ( $category->{'name'} eq $type ) {
                push @assemblers, $category->{'assemblers'}->[0];
            }
        }
        die "No valid assembler could be found" if ( $#assemblers == -1 );
    }

    # the scaffolder...
    # and the scaffolder...
    my @available_scaffolders = map { $_->[0] }
      sort { $a->[1] <=> $b->[1] }
      map { [ $_->{'name'}, $_->{'priority'} ] } @{ $config->{'scaffolders'} };
    my $linkage_evidence;

    if ($scaffolder) {
        unless ( grep ( /^$scaffolder$/i, @available_scaffolders ) ) {
            die "\nError: $scaffolder is not a valid scaffolder\n\n" . "The following scaffolders are configured: ",
              join( ", ", @available_scaffolders ), "\n\n";
        }
        my $conf_scaffolders = $config->{'scaffolders'};
        foreach my $conf_scaffolder (@$conf_scaffolders) {
            my $name = lc( $conf_scaffolder->{'name'} );
            if ( $name eq $scaffolder ) {
                if ( $conf_scaffolder->{'linkage_evidence'} eq 'paired-ends' & !$paired ) {
                    die "\nError: $scaffolder requires paired reads, but you only specified one fastq file...\n";
                }
                elsif ( $conf_scaffolder->{'linkage_evidence'} =~ /align/ & !$reference ) {
                    die "\nError: $scaffolder requires a reference for alignment, but none is specified...\n";
                }
                $linkage_evidence = $conf_scaffolder->{'linkage_evidence'};
            }

        }
    }

    my @available_mergers = map { $_->[0] }
      sort { $a->[1] <=> $b->[1] }
      map { [ $_->{'name'}, $_->{'priority'} ] } @{ $config->{'merge_tools'} };

    if ($merger) {
        unless ( grep ( /^$merger$/i, @available_mergers ) ) {
            die "\nError: $merger is not a valid merge-tool\n\n" . "The following tools are configured: ",
              join( ", ", @available_mergers ), "\n\n";
        }
        my $conf_mergers = $config->{'merge_tools'};
        foreach my $conf_merger (@$conf_mergers) {
            my $name = lc( $conf_merger->{'name'} );
            if ( $name eq $merger ) {
                undef($scaffolder) if $conf_merger->{'allow_scaffolding'} == 0;
            }
        }
    }

    my @available_finishers = map { $_->[0] }
      sort { $a->[1] <=> $b->[1] }
      map { [ $_->{'name'}, $_->{'priority'} ] } @{ $config->{'finishers'} };

    if ($finisher) {
        unless ( grep ( /^$finisher$/i, @available_finishers ) ) {
            die "\nError: $finisher is not a valid finishiing-tool\n\n" . "The following tools are configured: ",
              join( ", ", @available_finishers ), "\n\n";
        }
        my $conf_finishers = $config->{'finishers'};
        foreach my $conf_finisher (@$conf_finishers) {
            my $name = lc( $conf_finisher->{'name'} );
            if ( $name eq $finisher ) {
                die "Error: $finisher requires a reference sequence" unless ($reference);
                if ( $conf_finisher->{'paired_reads'} ) {
                    die "Error: $finisher requires paired reads" unless ($paired);
                }
            }
        }
    }

    my @available_varcallers = map { $_->[0] }
      sort { $a->[1] <=> $b->[1] }
      map { [ $_->{'name'}, $_->{'priority'} ] } @{ $config->{'varcallers'} };

    if ($varcall) {
        unless ( grep ( /^$varcall$/i, @available_varcallers ) ) {
            die "\nError: $varcall is not a valid variant caller\n\n" . "The following tools are configured: ",
              join( ", ", @available_varcallers ), "\n\n";
        }
        my $conf_varcallers = $config->{'varcallers'};
        foreach my $conf_varcaller (@$conf_varcallers) {
            my $name = lc( $conf_varcaller->{'name'} );
            if ( $name eq $varcall ) {
                die "Error: $varcall requires a reference sequence" unless ($reference);
            }
        }
    }
    print "linkage: $linkage_evidence";
    return ( \@assemblers, $scaffolder, $linkage_evidence, $merger, $finisher, $varcall );

}

#######################################################################
#
######################################################################
######################################################################
#
# return_results
#
# Copies results back from tmpdir, returning full working directory
# if dircopy argument specified
#
# required params: $ (tmpdir)
#                  $ (strain)
#                  $ (dircopy - flag to indicate entire directory
#                     should be returned)
#                  $ (mode - draft mode also needs amos bank copying)
#
# returns        : $ (0)
#
######################################################################

sub return_results {

    my $tmpdir  = shift;
    my $dir     = shift;
    my $prefix  = shift;
    my $dircopy = shift;
    my $mode    = shift;

    if ($dircopy) {
        dircopy( $tmpdir, "$dir" )
          or die "Error copying $tmpdir: $!";
    }
    else {
        my @files = qw(annotated.embl contigs.fasta scaffolds.fasta scaffolds.embl scaffolds.agp
          unplaced_contigs.fasta BugBuilder.log read1_fastqc.html read2_fastqc.html
          scaffolds_cgview.png contigs_cgview.png circleator.png circleator.svg reference.variants.vcf
          );

        opendir TMP, "$tmpdir" or die "Error opening $tmpdir: $!";
        my @all_files = readdir TMP;
        close TMP;

        foreach my $pattern (qw(blastout png)) {
            my @found = grep /$pattern/, @all_files;
            push @files, @found;
        }

        mkdir "$dir"
          or die "Error creating $dir: $!";

        foreach my $file (@files) {
            my $outfile;
            if ($prefix) {
                $outfile = "$dir/${prefix}_${file}";
            }
            else {
                $outfile = "$dir/$file";
            }
            my $target;
            if ( -l "$tmpdir/$file" ) {
                $target = readlink("$tmpdir/$file");
            }
            else {
                $target = "$tmpdir/$file";
            }
            copy( "$target", "$outfile" )
              or die "Error copying $file: $!"
              if ( -e "$tmpdir/$file" );
        }
        if ( $mode eq 'draft' ) {
            my $outfile;
            if ($prefix) {
                $outfile = "$dir/${prefix}_assembly.bnk";
            }
            else {
                $outfile = "$dir/assembly.bnk";
            }
            dircopy( "$tmpdir/amos/assembly.bnk", "$outfile" )
              or die "Error copying $tmpdir/amos/assembly.bnk: $!";
        }

    }

}

######################################################################
#
# run_fastqc
#
# Carries out QC analysis using fastqc...
#
# required params: $ (tmpdir)
#                  $ (type)
#
# returns: $ (0)
#
######################################################################

sub run_fastqc {

    my $tmpdir = shift;
    my $type   = shift;

    message("Running FastQC...");
    chdir $tmpdir          or die "Error chdiring to $tmpdir: $!";
    mkdir "$tmpdir/fastqc" or die "Error creating $tmpdir/fastqc: $!";

    my $fastq1     = readlink("$tmpdir/read1.fastq");
    my $fastq2     = readlink("$tmpdir/read2.fastq");
    my $long_fastq = readlink("$tmpdir/long.fastq");
    # my $de_fere_contigs = readlink("$tmpdir/de_fere_contigs.fasta");

    my $cmd = $config->{'fastqc_dir'} . "fastqc -t 4 -o $tmpdir/fastqc ";
    $cmd .= " $tmpdir/$fastq1"     if ($fastq1);
    $cmd .= " $tmpdir/$fastq2"     if ($fastq2);
    $cmd .= " $tmpdir/$long_fastq" if ($long_fastq);
    $cmd .= " >$tmpdir/fastqc/fastqc.log 2>&1";

    system($cmd) == 0      or die "Error running command $cmd: $!";
    chdir "$tmpdir/fastqc" or die "Error chdiring to $tmpdir/fastqc: $!";

    my $fail = 0;
    foreach my $file ( $fastq1, $fastq2, $long_fastq ) {
        if ($file) {
            $file =~ s/\.fastq(\.gz)?//;
            if ( -e "${file}_fastqc.zip" ) {
                my $archive = Archive::Zip->new("${file}_fastqc.zip");
                foreach my $member ( $archive->members ) {
                    next unless ( $member->fileName eq "${file}_fastqc/summary.txt" );
                    $member->extractToFileNamed("$tmpdir/${file}_fastqc.summary.txt");
                }
                open SUMMARY, "$tmpdir/${file}_fastqc.summary.txt" or die "Error opening $file.summary.txt: $!";
                print "\n\n$file FastQC report:\n";
                my $tb = Text::ASCIITable->new();
                $tb->setCols( "Result", "Metric" );
                while ( my $line = <SUMMARY> ) {
                    my @fields = split( /\t/, $line );
                    $tb->addRow( $fields[0], $fields[1] );
                    $fail++ if ( $line =~ /^FAIL/ );
                }
                close SUMMARY;
                print "$tb";
                copy( "$tmpdir/fastqc/${file}_fastqc.html", "$tmpdir/${file}_fastqc.html" )
                  or die "Error copying ${file}_fastqc.html: $!";
            }
        }
    }

    if ($fail) {
        if ( $type eq 'long' || $type eq 'hybrid' || $type eq 'de_fere' ) {
            print "NB: Reported quality issues with Per base sequence quality etc.\n";
            print "are normal when analysing PacBio sequence with FastQC...\n\n";
        }
        else {
            print "WARNING: some QC tests indicate quality issues with this data.\n";
            print "Please examine the fastqc outputs for these reads\n\n";
        }
    }

    chdir $tmpdir or die "Error chdiring to $tmpdir: $!";

    return (0);
}


######################################################################
#
# id_fastq_encoding
#
# Identifies fastq quality encoding type, based upon the information
# at http://en.wikipedia.org/wiki/FASTQ_format:
#
#  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
#  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
#  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
#  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
#  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
#  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
#  |                         |    |        |                              |                     |
# 33                        59   64       73                            104                   126
#
# S - Sanger        Phred+33,  raw reads typically (0, 40)
# X - Solexa        Solexa+64, raw reads typically (-5, 40)
# I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
# J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
#    with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
# L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
#
# required params: $ (tmpdir)
#
# returns:         $ (encoding)
#
######################################################################

sub id_fastq_encoding {

    my $tmpdir = shift;

    my $min = 104;
    my $max = 0;
    my $encoding;

    my $target;
    if ( -e "$tmpdir/long.fastq" ) {
        $target = readlink("$tmpdir/long.fastq");
    }
    else {
        $target = readlink("$tmpdir/read1.fastq");
    }

    if ( $target =~ /.gz$/ ) {
        open FASTQ, "<:gzip", "$tmpdir/$target" or die "Error opening $tmpdir/$target";
    }
    else {
        open FASTQ, "$tmpdir/$target" or die "Error opening $tmpdir/$target";
    }

    # looking at 1000 sequences should be plenty to work out the encoding
    for ( my $i = 0 ; $i <= 1000 ; $i++ ) {
        my $line;
        for ( my $j = 0 ; $j <= 3 ; $j++ ) {
            $line = <FASTQ>;
            if ( $j == 3 ) {
                if ($line) {
                    my @quals = split( //, $line );
                    foreach (@quals) {
                        my $value = ord($_);
                        $min = $value if ( $value < $min );
                        $max = $value if ( $value > $max );
                    }
                }
            }
        }
    }
    close FASTQ;

    if ( $min <= 59 && $max <= 74 ) {
        $encoding = 'sanger';
    }
    elsif ( $min > 59 && $min < 64 && $max < 104 ) {
        $encoding = 'solexa';
    }
    elsif ( $min > 64 ) {
        $encoding = 'illumina';
    }

    return ($encoding);

}

#######################################################################
#
# quality_trim_reads
#
# Quality trims reads using sickle
#
# required params: $ (tmpdir)
#                  $ (encoding)
#
# returns          $ (0)
#
#######################################################################

sub quality_trim_reads {

    my $tmpdir     = shift;
    my $encoding   = shift;
    my $qv         = shift;
    my $min_length = shift;

    if ( !$encoding ) {
        print "\n\nWARNING: Sequence quality encoding could not be determined\n";
        print "Sequence quality trimming will be skipped...\n\n";
    }
    else {
        mkdir "$tmpdir/qc_trim" or die "Error creating $tmpdir/qc_trim: $!";
        chdir "$tmpdir/qc_trim" or die "Error running chdir $tmpdir/qc_trim: $!";

        my $cmd = $config->{'sickle_dir'};
        if ( -e "$tmpdir/read2.fastq" ) {
            $cmd .=
                "sickle pe -f $tmpdir/read1.fastq "
              . "-r $tmpdir/read2.fastq -t $encoding -q $qv -l $min_length -o read1.fastq -p read2.fastq "
              . "-s singles.fastq > sickle.log";
        }
        else {
            $cmd .=
              "sickle se -f $tmpdir/read1.fastq " . "-t $encoding -q $qv -l $min_length -o read1.fastq > sickle.log";
        }

        system($cmd) == 0 or die "Error executing $cmd: $!";

        open LOG, "sickle.log" or die "Could not open sickle.log: $!";
        my $res;
        while ( my $line = <LOG> ) {
            $res .= $line;
        }
        close LOG;

        chdir $tmpdir or die "Error chdiring to $tmpdir: $!";

        unlink("$tmpdir/read1.fastq")
          or die "Error unlinking $tmpdir/read1.fastq: $!";
        symlink( "qc_trim/read1.fastq", "read1.fastq" )
          or die "Error creating read1.fastq symlink: $!";
        if ( -e "$tmpdir/read2.fastq" ) {
            unlink("$tmpdir/read2.fastq")
              or die "Error unlinking $tmpdir/read2.fastq: $!";
            symlink( "qc_trim/read2.fastq", "read2.fastq" )
              or die "Error creating read2.fastq symlink: $!";
        }

        print "\nQuality Trimming\n================\n";
        print $res;

        my ( $kept, $discarded );

        if ( -e "$tmpdir/read2.fastq" ) {
            $kept      = $1 if ( $res =~ /FastQ paired records kept: ([0-9]+)/ );
            $discarded = $1 if ( $res =~ /FastQ paired records discarded: ([0-9]+)/ );
        }
        else {
            $kept      = $1 if ( $res =~ /FastQ records kept: ([0-9]+)/ );
            $discarded = $1 if ( $res =~ /FastQ records discarded: ([0-9]+)/ );
        }

        print "\nWarning: >10% of reads discarded during read trimming.\nPlease examine the FastQC outputs... \n\n"
          if ( ( $discarded / $kept * 100 ) > 10 );
    }

    return ();
}

######################################################################
#
# downsample reads
#
# Downsamples reads by default to a maximum of 100x if higher
# coverages are present, or to specified coverage
#
# required params: $ (tmpdir)
#                  $ (config data)
#                  $ (coverage required)
#                  $ (original coverage)
#                  $ (number of bases in original reds)
#                  $ (mean read length)
#
# returns: 0
#
######################################################################

sub downsample_reads {

    my $tmpdir           = shift;
    my $config           = shift;
    my $coverage         = shift || 100;
    my $orig_coverage    = shift;
    my $read_bases       = shift;
    my $mean_read_length = shift;

    my $seqtk_dir = $config->{'seqtk_dir'};
    mkdir "$tmpdir/downsampling" or die "Error creating  $tmpdir/downsampling: $!";
    chdir "$tmpdir/"             or die "Error chdiring to $tmpdir/ downsampling : $! ";
    print "Projected Coverage: ${orig_coverage}x: Downsampling reads to ${coverage}x...\n";
    foreach my $fastq (qw(read1.fastq read2.fastq)) {
        if ( -e "$tmpdir/$fastq" ) {
            my $frac = ( $coverage / $orig_coverage );
            my $required_reads = ( ( $frac * $read_bases ) / $mean_read_length );
            # my $cmd = "$seqtk_dir/seqtk sample -s100 $tmpdir/$fastq $required_reads > $tmpdir/downsampling/$fastq";
            my $cmd = "$seqtk_dir/seqtk sample -s 100 $tmpdir/$fastq $required_reads > $tmpdir/downsampling/$fastq";
            system($cmd) == 0 or die " Error running $cmd: $! ";
        }
    }
    unlink("$tmpdir/read1.fastq")
      or die " Error unlinking $tmpdir/read1.fastq : $! ";
    symlink( "downsampling/read1.fastq", "read1.fastq" )
      or die " Error creating read1.fastq symlink : $! ";
    if ( -e "$tmpdir/read2.fastq" ) {
        unlink("$tmpdir/read2.fastq")
          or die "Error unlinking $tmpdir/read2.fastq: $! ";
        symlink( "downsampling/read2.fastq", "read2.fastq" )
          or die " Error creating read2.fastq symlink: $! ";
    }
    return ($coverage);
}

#######################################################################
#
# run_assembler
#
# Runs specified assembler on fastq files
#
# required params: $ (tmpdir)
#                  $ (assembler name)
#                  $ (reference fasta sequece - probably not needed)
#                  $ (arguments to pass to assembler)
#                  $ (link - flag to indicate contigs should be symlinked into tmpdir)
#                  $ (category - type of assembly to run, since an assembler may fall into more
#                  than one category)
#                  $ (encoding - some assemblers need explicitly telling)
#                  $ (genome size)
#                  $ (average read length)
#	           $ (insert size)
#                  $ (insert size stddev)
#
# returns          $ (0)
#
#######################################################################

sub run_assembler {

    my $tmpdir      = shift;
    my $assembler   = shift;
    my $reference   = shift;
    my $args        = shift;
    my $link        = shift;
    my $category    = shift;
    my $encoding    = shift;
    my $genome_size = shift;
    my $platform    = shift;
    my $read_length = shift;
    my $insert_size = shift;
    my $stddev      = shift;
    my $threads     = shift;

    my ( $cmd, $contig_output, $scaffold_output, $create );

    my $assemblers = $config->{'assemblers'};
    foreach my $conf_assembler (@$assemblers) {
        if ( $conf_assembler->{'name'} eq $assembler ) {
            if ( $category eq 'hybrid' ) {
                $cmd = $conf_assembler->{'command_hybrid'};
            }
            elsif ( $category eq 'de_fere' ) {
                $cmd = $conf_assembler->{'command_de_fere'};
            }
            elsif ( -e "$tmpdir/read2.fastq" ) {
                $cmd = $conf_assembler->{'command_pe'};
            }
            else {
                $cmd = $conf_assembler->{'command_se'};
            }
            if ($args) {
                $cmd .= " $args";
            }
            elsif ( $conf_assembler->{'default_args'} ) {
                $cmd .= " " . $conf_assembler->{'default_args'};
            }
            $contig_output   = $conf_assembler->{'contig_output'};
            $scaffold_output = $conf_assembler->{'scaffold_output'}
              if ( $conf_assembler->{'scaffold_output'} );
            $create = $conf_assembler->{'create_dir'};
        }
    }
    die "Assembler $assembler is not defined" unless ($cmd);

    if ($create) {
        mkdir "$tmpdir/$assembler"
          or die "could not create $tmpdir/$assembler: $! ";
        chdir "$tmpdir/$assembler"
          or die "could not chdir to $tmpdir/$assembler: $! ";
    }

    if ( $reference && !$genome_size ) {
        my $io = Bio::SeqIO->new( -file => "$tmpdir/$reference", -format => 'fasta' );
        while ( my $ref = $io->next_seq() ) {
            $genome_size += $ref->length();
        }
    }
    my $asm_dir = $config->{"${assembler}_dir"};
    foreach ( $cmd, $contig_output, $scaffold_output ) {
        s/__BUGBUILDER_BIN__/$FindBin::Bin/g;
        s/__ASMDIR__/$asm_dir/;
        s/__TMPDIR__/$tmpdir/g;
        s/__FASTQ1__/$tmpdir\/read1.fastq/;
        s/__FASTQ2__/$tmpdir\/read2.fastq/;
        s/__ORIG_FASTQ1__/$tmpdir\/orig_read1.fastq/;
        s/__ORIG_FASTQ2__/$tmpdir\/orig_read2.fastq/;
        s/__DE_FERE_CONTIGS__/$tmpdir\/de_fere_contigs.fasta/;
        s/__LONGFASTQ__/$tmpdir\/long.fastq/;
        s/__REFERENCE__/$reference/ if ($reference); #
        s/__CATEGORY__/$category/;
        s/__ENCODING__/$encoding/;
        s/__GENOME_SIZE__/$genome_size/ if defined($genome_size);
        s/__PLATFORM__/$platform/       if ($platform);
        s/__READ_LENGTH__/$read_length/ if ($read_length);
        s/__INSSIZE__/$insert_size/     if ($insert_size);
        s/__INSSD__/$stddev/            if ($stddev);
        s/__THREADS__/$threads/;
    }

    $cmd .= " > $tmpdir/$assembler.log 2>&1 ";

    message(" Starting $assembler assembly ... ");
    system($cmd) == 0 or die "Error running $cmd: $! ";

    my $header = "$assembler assembly statistics";
    print $header . "\n" . '=' x length($header) . "\n\n";
    get_contig_stats( "$contig_output", 'contigs' );

    # Only generate scaffold stats for paired read alignments...
    get_contig_stats( "$scaffold_output", 'scaffolds' )
      if ( $scaffold_output && -e $scaffold_output && -e "$tmpdir/read2.fastq" );

    # rename contigs/scaffolds for consistent naming, since we need to retrieve by
    # id later, so it helps if we know what the ids look like...
    if ( !$create ) {
        my $outdir = dirname($contig_output);
        chdir $outdir or die " Error changing to dir $outdir: $! ";
    }

    my $inIO = Bio::SeqIO->new( -format => 'fasta', -file => $contig_output );
    my $outIO =
      Bio::SeqIO->new( -format => 'fasta',
                       -file   => ">BugBuilder.contigs.fasta" );
    my $contig_count = 0;
    while ( my $seq = $inIO->next_seq() ) {
        my $contig_id = sprintf( "contig_%06s", ++$contig_count );
        $seq->id($contig_id);
        $outIO->write_seq($seq);
    }

    if ( $scaffold_output && ( -e $scaffold_output ) ) {    #&& ( -e "$tmpdir/read2.fastq" ) ) {
        my $inIO = Bio::SeqIO->new( -format => 'fasta', -file => $scaffold_output );
        my $outIO =
          Bio::SeqIO->new( -format => 'fasta',
                           -file   => ">BugBuilder.scaffolds.fasta" );
        my $scaffold_count = 0;
        while ( my $seq = $inIO->next_seq() ) {
            my $scaffold_id = sprintf( "scaffold_%06s", ++$scaffold_count );
            $seq->id($scaffold_id);
            $outIO->write_seq($seq);
        }
    }

    if ($link) {

        # need to use File::Find::Rule to get the right path since some outputs are more nested
        # than others...
        my $contigs = ( File::Find::Rule->file()->name("BugBuilder.contigs.fasta")->in("$tmpdir/$assembler") )[0];
        symlink( $contigs, "../contigs.fasta" )
          or die "Error creating symlink: $!";
        if ( $scaffold_output && ( -e $scaffold_output ) ) {
            my $scaffolds =
              ( File::Find::Rule->file()->name("BugBuilder.scaffolds.fasta")->in("$tmpdir/$assembler") )[0];
            symlink( $scaffolds, "../scaffolds.fasta" )
              or die "Error creating symlink: $!"
              if ($scaffolds);
        }
    }

    chdir $tmpdir or die "Could not chdir to $tmpdir: $!" if ($create);

    return ();
}

######################################################################
#
# merge_assemblies
#
# combines two assemblies using the selected merge-method
#
# required params: $ (tmpdir)
#                  $ (arrayref of assemblers used)
#                  $ (method)
#                  $ (reference)
#
# returns        : $ (0)
#
######################################################################

sub merge_assemblies {

    my $tmpdir     = shift;
    my $assemblers = shift;
    my $method     = shift;
    my $reference  = shift;

    message("Merging assemblies ($method)...");

    my ( $cmd, $create_dir, $contig_output );
    my $merge_tools = $config->{'merge_tools'};
    foreach my $tool (@$merge_tools) {
        my $name = $tool->{'name'};
        if ( $name eq $method ) {
            $cmd           = $tool->{'command'};
            $create_dir    = $tool->{'create_dir'};
            $contig_output = $tool->{'contig_output'};
        }
        if ($create_dir) {
            mkdir "$tmpdir/$method" or die "Error creating $tmpdir/$method: $! ";
            chdir "$tmpdir/$method" or die "Error chdiring to $tmpdir/$method: $!";
        }

        $cmd =~ s/__BUGBUILDER_BIN__/$FindBin::Bin/;
        $cmd =~ s/__TMPDIR__/$tmpdir/;
        $cmd =~ s/__ASSEMB1__/$assemblers->[0]/;
        $cmd =~ s/__ASSEMB2__/$assemblers->[1]/;
        $cmd =~ s/__REFERENCE__/${tmpdir}\/reference.fasta/;

        system($cmd) == 0 or die "Error running $cmd: $!";
        chdir $tmpdir     or die " Error chdiring to $tmpdir: $! ";
        symlink( "$method/$contig_output", "contigs.fasta" )
          or die " Error creating symlink : $! ";

        print "Merged assembly statistics : \n============================\n";
        get_contig_stats( "$tmpdir/$method/$contig_output", 'contigs' );

        return (0);

    }
}

######################################################################
#
# finish_assembly
#
# Carried out assembly finishing using selected method
#
# required params: $ (tmpdir)
#                  $ (finisher)
#                  $ (insert size)
#                  $ (insert stddev)
#                  $ (base quality encoding)
#                  $ (no. threads)
#
# returns        : $ (0)
#
######################################################################

sub finish_assembly {

    my $tmpdir        = shift;
    my $finisher      = shift;
    my $insert_size   = shift;
    my $insert_stddev = shift;
    my $encoding      = shift;
    my $threads       = shift;

    message("Finishing assembly ($finisher)...");

    my ( $cmd, $create_dir );
    my $finishers = $config->{'finishers'};
    foreach my $tool (@$finishers) {
        my $name = $tool->{'name'};
        if ( $name eq $finisher ) {
            $cmd        = $tool->{'command'};
            $create_dir = $tool->{'create_dir'};
        }
    }
    if ($create_dir) {
        mkdir "$tmpdir/$finisher" or die "Error creating $tmpdir/$finisher: $! ";
        chdir "$tmpdir/$finisher" or die "Error chdiring to $tmpdir/$finisher: $!";
    }

    $cmd =~ s/__BUGBUILDER_BIN__/$FindBin::Bin/;
    $cmd =~ s/__TMPDIR__/$tmpdir/;
    $cmd =~ s/__REFERENCE__/${tmpdir}\/reference.fasta/;
    $cmd =~ s/__INSSIZE__/$insert_size/;
    $cmd =~ s/__INSSD__/$insert_stddev/;
    $cmd =~ s/__ENCODING__/$encoding/;
    $cmd =~ s/__THREADS__/$threads/;

    system("perl $cmd") == 0 or die "Error running $cmd: $!";
    chdir $tmpdir     or die " Error chdiring to $tmpdir: $! ";

    print "Finished assembly statistics : \n============================\n";
    get_contig_stats( "$tmpdir/scaffolds.fasta", 'scaffolds' );

    return (0);
}

######################################################################
#
# align_reads
#
# Maps reads to reference using bwa
#
# required params: $ (tmp directory);
#                  $ (reference);
#                  $ (read length)
#                  $ (flag to indicate downsampling required...)
#
#                : $ (0)
#
######################################################################

sub align_reads {

    my $tmpdir      = shift;
    my $reference   = shift;
    my $read_length = shift;
    my $downsample  = shift;

    my $bwa_dir      = $config->{'bwa_dir'};
    my $samtools_dir = $config->{'samtools_dir'};
    my $seqtk_dir    = $config->{'seqtk_dir'};

    mkdir "$tmpdir/bwa"
      or die "Error creating $tmpdir/bwa: $!"
      if ( !-d "$tmpdir/bwa" );
    chdir "$tmpdir/bwa" or die "Error chdiring to $tmpdir/bwa: $!";
    copy( "$tmpdir/$reference", "$tmpdir/bwa/$reference" )
      or die "Error copying $tmpdir/$reference: $! ";

    if ($downsample) {

        message("Downsampling reads for insert-size estimation...");
        my $cmd = "$seqtk_dir/seqtk sample -s100 $tmpdir/read1.fastq 10000 > $tmpdir/bwa/read1.fastq";
        system($cmd) == 0 or die " Error running $cmd: $! ";
        if ( -e "$tmpdir/read2.fastq" ) {
            $cmd = "$seqtk_dir/seqtk sample -s100 $tmpdir/read2.fastq 10000 > $tmpdir/bwa/read2.fastq";
            system($cmd) == 0 or die " Error running $cmd: $!";

        }
    }
    else {
        symlink( "$tmpdir/read1.fastq", "$tmpdir/bwa/read1.fastq" ) or die "Error creating symlink: $! ";
        symlink( "$tmpdir/read2.fastq", "$tmpdir/bwa/read2.fastq" ) or die "Error creating symlink: $! ";
    }

    message("BWA aligning reads vs $reference...");

    my $cmd = "$bwa_dir/bwa index $tmpdir/bwa/$reference >$tmpdir/bwa/bwa_index.log 2>&1";
    system($cmd) == 0 or die " Error running $cmd";

    # Use bwa-bwt for 'short' reads less than 100 bp, and bwa-mem for longer reads
    if ( $read_length <= 100 ) {
        $cmd = "$bwa_dir/bwa aln -t 4 $tmpdir/bwa/$reference $tmpdir/bwa/read1.fastq > $tmpdir/bwa/read1.sai"
          . " 2> $tmpdir/bwa/bwa_sai1.log";
        system($cmd) == 0 or die "Error running $cmd";
        if ( -e "$tmpdir/read2.fastq" ) {
            $cmd = "$bwa_dir/bwa aln -t 4 $tmpdir/bwa/$reference $tmpdir/bwa/read2.fastq > $tmpdir/bwa/read2.sai"
              . " 2> $tmpdir/bwa/bwa_sai2.log";
            system($cmd) == 0 or die "Error running $cmd";
            $cmd = "$bwa_dir/bwa sampe $tmpdir/bwa/$reference $tmpdir/bwa/read1.sai $tmpdir/bwa/read2.sai "
              . "$tmpdir/bwa/read1.fastq $tmpdir/bwa/read2.fastq";
            $cmd .= " 2> $tmpdir/bwa/sampe.log > $tmpdir/bwa/$reference.sam";
            system($cmd) == 0 or die "Error running $cmd";
        }
        else {
            $cmd = "$bwa_dir/bwa samse $tmpdir/bwa/$reference $tmpdir/bwa/read1.sai $tmpdir/read1.fastq";
            $cmd .= "2> $tmpdir/bwa/samse.log > $tmpdir/bwa/$reference.sam";
            system($cmd) == 0 or die "Error running $cmd";
        }
    }
    else {
        if ( !-e "$tmpdir/read2.fastq" ) {

            # single-ended long reads
            $cmd = "$bwa_dir/bwa mem -t 4 -M $tmpdir/bwa/$reference $tmpdir/bwa/read1.fastq > $reference.sam "
              . "2>$tmpdir/bwa/bwa_mem.log";
            system($cmd) == 0 or die "Error running $cmd: $!";
        }
        else {

            # paired-end long reads
            $cmd =
"$bwa_dir/bwa mem -t 4 -M $tmpdir/bwa/$reference $tmpdir/bwa/read1.fastq $tmpdir/bwa/read2.fastq >$reference.sam "
              . "2>$tmpdir/bwa/bwa_mem.log";
            system($cmd) == 0 or die "Error running $cmd: $!";
        }
    }

    $cmd = "$samtools_dir/samtools view -q 10 -Sb $tmpdir/bwa/$reference.sam 2>$tmpdir/bwa/samtoolsview.log"
      # . "|$samtools_dir/samtools sort - $tmpdir/bwa/$reference";
      . "|$samtools_dir/samtools sort - > $tmpdir/bwa/$reference.bam";
    system($cmd) == 0 or die "Error running $cmd";

    $cmd = "$samtools_dir/samtools index $tmpdir/bwa/$reference.bam 2>$tmpdir/bwa/samtools_index.log";
    system($cmd) == 0                 or die "Error running $cmd";
    unlink("$tmpdir/bwa/read1.fastq") or die " Error unlinking $tmpdir/bwa/read1.fastq: $! ";
    unlink("$tmpdir/bwa/read2.fastq")
      or die " Error unlinking $tmpdir/bwa/read1.fastq : $! "
      if ( -e "$tmpdir/bwa/read2.fastq" );

    chdir($tmpdir) or die "Error changing to $tmpdir: $! ";

}

######################################################################
#
# get_insert_stats
#
# converts contigs.sam -> bam, sorts, indexes and generates
# insert stats with Picard
#
# required params: $ (tmp directory);
#                  $ (reference);
#
# returns        : $ (insert size)
#                : $ (stddev)
#
######################################################################

sub get_insert_stats {

    my $tmpdir    = shift;
    my $reference = shift;

    message(" Collecting insert size statistics ");

    mkdir "$tmpdir/insert_stats"
      or die "Error creating $tmpdir/insert_stats: $! "
      if ( !-d "$tmpdir/insert_stats" );
    chdir "$tmpdir/insert_stats"
      or die "Error chdiring to $tmpdir/insert_stats : $! ";

    # picard command for earlier versions....
    #my $cmd = "java -jar " . $config->{'picard_dir'} . "CollectInsertSizeMetrics.jar ";
    my $cmd = "java -jar " . $config->{'picard_dir'} . "picard.jar CollectInsertSizeMetrics ";
    $cmd .= "INPUT=$tmpdir/bwa/$reference.bam HISTOGRAM_FILE=insert_histogram.pdf OUTPUT=insert_stats.txt ";
    $cmd .=
      " QUIET=true VERBOSITY=ERROR ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT > CollectInsertMetrics.log 2>&1";
    system($cmd) == 0 or die " Error executing $cmd: $! ";

    open STATS, "insert_stats.txt"
      or die "Error opening insert_stats.txt: $! ";
    my ( $min_insert, $max_insert, $insert, $stddev );
  LINE: while ( my $line = <STATS> ) {
        if ( $line =~ /^MEDIAN/ ) {
            my @stats = split( /\t/, <STATS> );
            $min_insert = $stats[2];
            $max_insert = $stats[3];
            $insert     = sprintf( "%d", $stats[4] );
            $stddev     = sprintf( "%d", $stats[5] );
            last LINE;
        }
    }
    close STATS;

    chdir $tmpdir or die "Error chdiring to $tmpdir: $!";

    return ( $insert, $stddev );

}

######################################################################
#
# check_id
#
# Checks wether the assembled contigs have sufficient identity to the
# reference for reference-based scaffolding
#
# required params: $ (tmp directory);
#
# returns        : $ id_ok (boolean)
#
######################################################################

sub check_id {

    my $tmpdir = shift;

    message("Checking identity of assembly with reference...");

    my $id_ok = 0;

    mkdir "$tmpdir/id_check"
      or die "Error creating $tmpdir/id_check: $! "
      if ( !-d "$tmpdir/id_check" );
    chdir "$tmpdir/id_check"
      or die "Error chdiring to $tmpdir/id_check: $! ";

    my $cmd = $config->{'blast_dir'}
      . "/blastn -query $tmpdir/contigs.fasta -subject $tmpdir/reference_parsed_ids.fasta -outfmt 5 -evalue 0.01 -out blastout.xml 2>&1 > blastn.log";
    system($cmd) == 0 or die "Error running $cmd: $!";

    my $blio = Bio::SearchIO->new( -format => 'blastxml',
                                   -file   => 'blastout.xml' );

    my ( $aligned, $unaligned );
    while ( my $result = $blio->next_result() ) {
        foreach my $hit ( $result->hits() ) {
            my $tiling = Bio::Search::Tiling::MapTiling->new($hit);
            $aligned   += $tiling->num_aligned();
            $unaligned += $tiling->num_unaligned();
        }
    }
    my $percent_id = sprintf( '%.2f', ( $aligned / ( $aligned + $unaligned ) * 100 ) );
    print "ID=$percent_id %\n";

    if ( $percent_id > 80 ) {
        $id_ok = 1;
    }
    else {
        $id_ok = 0;
        print RED, "WARNING", RESET, ": Percentage identitiy of assembly with reference looks to low (${percent_id}%\n";
        print "Reference will not be used for scaffolding or ordering contigs\n";
    }

    chdir $tmpdir or die "Error chdiring to $tmpdir: $!";

    return ($id_ok);

}

######################################################################
#
# run_scaffolder
#
# Runs specified scaffolder....
#
# Scaffolders which don't separate unscaffolded contigs can be wrapped
# in a script which created a "$scaffolder.contig_ids" file listing the
# IDs of contigs scaffolded. These will then be used following scaffolding
# to create our own file of unscaffolded contigs
#
# required params: $ (tmpdir)
#                  $ (reference)
#                  $ (scaffolder)
#                  $ (scaffolder args)
#                  $ (library insert size)
#	           $ (library insert sd)
#                  $ (run_id - appended to tmpdir to allow multiple runs)
#                  $ (path to contigs to scaffold)
#                  $ (mean read length)
#
# returns        : $ (linkage evidence type)
#
######################################################################

sub run_scaffolder {

    my $tmpdir           = shift;
    my $reference        = shift;
    my $scaffolder       = shift;
    my $scaffolder_args  = shift;
    my $insert_size      = shift;
    my $insert_sd        = shift;
    my $run_id           = shift;
    my $contigs          = shift;
    my $mean_read_length = shift;
    my $threads          = shift;

    message(" Starting $scaffolder");

    my ( $cmd, $scaffold_output, $unscaffolded_output, $create, $linkage_evidence, $default_args );
    my $blast_dir = $config->{'blast_dir'};

    my $scaffolders = $config->{'scaffolders'};
    foreach my $conf_scaffolder (@$scaffolders) {
        if ( lc( $conf_scaffolder->{'name'} ) eq lc($scaffolder) ) {
            $cmd                 = $conf_scaffolder->{'command'};
            $scaffold_output     = $conf_scaffolder->{'scaffold_output'};
            $unscaffolded_output = $conf_scaffolder->{'unscaffolded_output'}
              if ( $conf_scaffolder->{'unscaffolded_output'} );
            $create           = $conf_scaffolder->{'create_dir'};
            $linkage_evidence = $conf_scaffolder->{'linkage_evidence'};
            $default_args     = $conf_scaffolder->{'default_args'};
        }
    }
    die "Scaffolder $scaffolder is not defined" unless ($cmd);
    my $run_dir .= "${tmpdir}/${scaffolder}_${run_id}";
    mkdir("$run_dir") or die "Error creating $run_dir: $! ";
    chdir("$run_dir") or die "Error in chdir $run_dir: $! ";

    # Treat reference-based scaffolder separately from paired-read scaffolders,
    # since we need to scaffold per-reference, which doesn't work if your not using one...
    if ( $linkage_evidence eq 'align_genus' ) {

        # If the reference contains multiple contigs, we first neeed to align out contigs to these
        # to identify which contigs to scaffold against which reference, since
        # some scaffolders targeted at bacteria don't handle multiple reference
        # sequences
        my $io = Bio::SeqIO->new( -format => 'fasta', -file => "$reference" );
        my @ref_ids;
        while ( my $seq = $io->next_seq() ) {
            my $ref_id = $seq->display_id();

            #$ref_id =~ s/lcl\|//;
            $ref_id = parse_ref_id($ref_id);

            #push @ref_ids, $seq->id();
            push @ref_ids, $ref_id;
            my $outIO = Bio::SeqIO->new( -format => 'fasta', -file => ">$run_dir/reference_${ref_id}" );
            $outIO->write_seq($seq);
        }

        if ( $#ref_ids > 0 ) {

            # blast indexing doesn't produce a workable index from a symlink, so need to copy the reference sequences
            copy( $reference, "$run_dir/reference.fasta" )
              or die "Error copying reference.fasta:$! ";
            symlink( "$contigs", "$run_dir/contigs.fasta" );
            my $fasta_db = Bio::DB::Fasta->new("$run_dir/contigs.fasta");

            my $blast_cmd =
              "$blast_dir/makeblastdb -in reference.fasta -dbtype nucl -parse_seqids 2>&1 > makeblastdb.log";
            system($blast_cmd) == 0 or die "Error building reference.fasta blast database: $!";
            $blast_cmd = "$blast_dir/blastn -query $contigs -task blastn -db reference.fasta "
              . "-out clusters.blast  2>&1 > blastn.log";
            system($blast_cmd) == 0 or die "Error executing $blast_cmd: $!";

            # create hash of contigs per reference sequence, or unaligned
            # only need to worry about the top hit for each...
            my %ref_seqs;
            $ref_seqs{'unaligned'} = [];
            my $io = Bio::SearchIO->new( -format => 'blast', -file => 'clusters.blast' );
            while ( my $result = $io->next_result() ) {
                my $top_hit = $result->next_hit();
                my $query   = $result->query_name();
                if ($top_hit) {
                    my $ref_seq = $top_hit->name();
                    if ( $ref_seqs{$ref_seq} ) {
                        push( @{ $ref_seqs{$ref_seq} }, $query );
                    }
                    else {
                        $ref_seqs{$ref_seq} = [$query];
                    }
                }
                else {
                    push( @{ $ref_seqs{'unaligned'} }, $query );
                }
            }

            # generate a fasta file of contigs which align to each reference
            foreach my $ref ( keys(%ref_seqs) ) {
                my $out_name = "$run_dir/reference_${ref}_contigs";
                $out_name =~ s/lcl\|//;
                my $outIO = Bio::SeqIO->new( -format => 'fasta', -file => ">$out_name" );
                foreach my $id ( @{ $ref_seqs{$ref} } ) {
                    my $contig = $fasta_db->get_Seq_by_id($id);
                    $outIO->write_seq($contig);
                }
            }

        }
        else {

            # if we don't have mulitple references, we just need to make the
            # reference and contigs available under consistent names
            symlink( $reference, "$run_dir/reference_" . $ref_ids[0] )
              or die "Error creating $run_dir/reference_${ref_ids[0]} symlink: $!"
              if ( !-e "$run_dir/reference_" . $ref_ids[0] );

            #my $src_contigs = "$run_dir/" . basename($contigs);
            symlink( $contigs, "$run_dir/reference_${ref_ids[0]}_contigs" )
              or die "Error creating $run_dir/reference_${ref_ids[0]}_contigs symlink";
        }

        my $mergedIO = Bio::SeqIO->new( -format => 'fasta', -file => ">$run_dir/scaffolds.fasta" );
        my $merged_scaffolds;    #for final merging of per-reference scaffolds
        my $merged_scaff_count = 0;

        # Now run the selecting scaffolder on each set of reference and contigs...
        foreach my $ref_id (@ref_ids) {

            my $scaff_contigs = "reference_${ref_id}_contigs";

            # Only if anything mapped to this reference
            next unless ( -e "${run_dir}/${scaff_contigs}" );

            print "\nScaffolding vs. $ref_id\n";

            my $exec_cmd  = $cmd;
            my $reference = "reference_${ref_id}";

            my @replace = ( $cmd, $scaffold_output );
            foreach ( $exec_cmd, $scaffold_output ) {
                s/__BUGBUILDER_BIN__/$FindBin::Bin/;
                s/__TMPDIR__/$tmpdir/g;
                s/__SCAFFDIR__/$run_dir\/${scaffolder}_${ref_id}/g;
                s/__RUN__/$run_id/;
                s/__FASTQ1__/$tmpdir\/read1.fastq/;
                s/__FASTQ2__/$tmpdir\/read2.fastq/;
                s/__REFERENCE__/${run_dir}\/${reference}/;
                s/__CONTIGS__/${run_dir}\/${scaff_contigs}/;
                s/__INSSIZE__/$insert_size/;
                s/__INSSD__/$insert_sd/;
                s/__THREADS__/$threads/;
            }
            if ( defined($unscaffolded_output) ) {
                $unscaffolded_output =~ s/__BUGBUILDER_BIN__/$FindBin::Bin/;
                $unscaffolded_output =~ s/__TMPDIR__/$tmpdir/g;
                $unscaffolded_output =~ s/__SCAFFDIR__/$run_dir\/$ref_id/g;
                $unscaffolded_output =~ s/__RUN__/$run_id/g;
                $unscaffolded_output =~ s/__FASTQ1__/$tmpdir\/read1.fastq/;
                $unscaffolded_output =~ s/__FASTQ2__/$tmpdir\/read2.fastq/;
                $unscaffolded_output =~ s/__REFERENCE__/${run_dir}_${reference}/;
                $unscaffolded_output =~ s/__CONTIGS__/$run_dir\/$scaff_contigs/;
                $unscaffolded_output =~ s/__INSSIZE__/$insert_size/;
                $unscaffolded_output =~ s/__INSSD__/$insert_sd/;
            }
            my $run_scaffold_output = "$run_dir/${scaffolder}_${ref_id}/$scaffold_output";

            if ($scaffolder_args) {
                $exec_cmd .= "$scaffolder_args";
            }
            elsif ($default_args) {
                $exec_cmd .= "$default_args" if ($default_args);
            }

            $exec_cmd .= " 2>&1 >$tmpdir/${scaffolder}_${run_id}_${ref_id}.log";

            mkdir("$run_dir/${scaffolder}_${ref_id}") or die "Error creating $run_dir/${scaffolder}_${ref_id}: $! ";
            chdir("$run_dir/${scaffolder}_${ref_id}") or die "Error in chdir $run_dir/${scaffolder}_${ref_id}: $! ";
            symlink( "$run_dir/${reference}", "$run_dir/${scaffolder}_${ref_id}/$reference" )
              or die "Error linking $reference:$! ";

            #symlink( "$run_dir/reference_${ref_ids[0]}_contigs", "$run_dir/${scaffolder}_${ref_id}/$contigs" )
            #or die "Error linking $contigs:$! ";

            system("perl $exec_cmd") == 0 or die "Error executing $exec_cmd: ";
            my $scaffIO =
              Bio::SeqIO->new( -format => 'fasta', -file => "$run_dir/${scaffolder}_${ref_id}/scaffolds.fasta" )
              or die "Error opening $run_dir/${scaffolder}_${ref_id}/scaffolds.fasta: $!";
            while ( my $scaff = $scaffIO->next_seq() ) {
                $scaff->display_id( 'scaffold_' . ++$merged_scaff_count );
                $mergedIO->write_seq($scaff);
            }
        }
    }
    else {    #non reference-guided scaffolding....

        my $exec_cmd = $cmd;

        # A kludge to work when tmpdir is not the top run dir - needed
        # because align_reads concatenates path from tmpdir and contigs.fasta
        if ( !-e "$tmpdir/contigs.fasta" ) {
            symlink( $contigs, "$tmpdir/contigs.fasta" );
        }

        # if no reference provided we won't have an estimate of insert size,
        # so need to get this by read alignment vs the assembly.
        if ( !$insert_size ) {
            align_reads( $tmpdir, "contigs.fasta", $mean_read_length, 1 );
            ( $insert_size, $insert_sd ) = get_insert_stats( "$tmpdir", "contigs.fasta" );
        }

        my @replace = ( $cmd, $scaffold_output );
        foreach ( $exec_cmd, $scaffold_output ) {
            s/__BUGBUILDER_BIN__/$FindBin::Bin/;
            s/__TMPDIR__/$tmpdir/g;
            s/__SCAFFDIR__/$run_dir/g;
            s/__RUN__/$run_id/;
            s/__FASTQ1__/$tmpdir\/read1.fastq/;
            s/__FASTQ2__/$tmpdir\/read2.fastq/;
            s/__CONTIGS__/$contigs/;
            s/__INSSIZE__/$insert_size/;
            s/__INSSD__/$insert_sd/;
            s/__THREADS__/$threads/;
        }
        if ( defined($unscaffolded_output) ) {
            $unscaffolded_output =~ s/__BUGBUILDER_BIN__/$FindBin::Bin/;
            $unscaffolded_output =~ s/__TMPDIR__/$tmpdir/g;
            $unscaffolded_output =~ s/__SCAFFDIR__/$run_dir/g;
            $unscaffolded_output =~ s/__RUN__/$run_id/g;
            $unscaffolded_output =~ s/__FASTQ1__/$tmpdir\/read1.fastq/;
            $unscaffolded_output =~ s/__FASTQ2__/$tmpdir\/read2.fastq/;
            $unscaffolded_output =~ s/__CONTIGS__/$contigs/;
            $unscaffolded_output =~ s/__INSSIZE__/$insert_size/;
            $unscaffolded_output =~ s/__INSSD__/$insert_sd/;
        }
        my $run_scaffold_output = "$run_dir/${scaffolder}/$scaffold_output";

        if ($scaffolder_args) {
            $exec_cmd .= "$scaffolder_args";
        }
        elsif ($default_args) {
            $exec_cmd .= "$default_args" if ($default_args);
        }

        $exec_cmd .= " 2>&1 >$tmpdir/${scaffolder}_${run_id}.log";

        system("perl $exec_cmd") == 0 or die "Error executing $exec_cmd:$!";

        #mkdir("$run_dir/${scaffolder}") or die "Error creating $run_dir/${scaffolder}: $! ";
        #chdir("$run_dir/${scaffolder}") or die "Error in chdir $run_dir/{scaffolder}: $! ";
        #symlink( "$run_dir/${reference}", "$run_dir/${scaffolder}/$reference" )
        #  or die "Error linking $reference:$! ";

        #symlink( "$run_dir/reference_${ref_ids[0]}_contigs", "$run_dir/${scaffolder}_${ref_id}/$contigs" )

    }

    # Create a fasta file of unplaced contigs if files of contig_ids are generated by the scaffolder wrapper
    my @id_files = File::Find::Rule->file()->name("${scaffolder}.contig_ids")->in($run_dir);
    if ( scalar(@id_files) ) {
        my %used_contigs;
        foreach my $id_file (@id_files) {
            open CONTIG_IDS, $id_file or die "Error opening $id_file: $!";
            while (<CONTIG_IDS>) {
                chomp;
                $used_contigs{$_}++;
            }
            close CONTIG_IDS;
        }
        my $inIO  = Bio::SeqIO->new( -format => 'fasta', -file => "$contigs" );
        my $outIO = Bio::SeqIO->new( -format => 'fasta', -file => ">$run_dir/unplaced_contigs.fasta" );
        while ( my $contig = $inIO->next_seq() ) {
            $outIO->write_seq($contig) unless ( $used_contigs{ $contig->display_id() } && $contig->length() > 200 );
        }
    }

    #renumber scaffolds to ensure they are unique...
    my $count = 0;
    my $inIO  = Bio::SeqIO->new( -file => "$run_dir/scaffolds.fasta", -format => "fasta" );
    my $outIO = Bio::SeqIO->new( -file => ">$run_dir/scaffolds_renumbered.fasta", -format => "fasta" );

    while ( my $seq = $inIO->next_seq() ) {
        $seq->display_id( "scaffold_" . ++$count );
        $outIO->write_seq($seq);
    }

    print "\nScaffolded assembly stats\n=========================\n\n";
    get_contig_stats( "$run_dir/scaffolds.fasta", 'scaffolds' );

    if ( $run_id == 1 ) {
        chdir $tmpdir or die "Error chdiring to $tmpdir: $! ";
        unlink "scaffolds.fasta"
          or die "Error removing scaffolds.fasta: $! "
          if ( -l "scaffolds.fasta" );
        symlink( "$run_dir/scaffolds_renumbered.fasta", "scaffolds.fasta" )
          or die "Error creating symlink: $! ";
        if ( defined($unscaffolded_output) ) {
            symlink( "$run_dir/$unscaffolded_output", "unplaced_contigs.fasta" )
              or die "Error creating symlink: $!";
        }
        elsif ( -e "$tmpdir/${scaffolder}/unplaced_contigs.fasta" ) {
            symlink( "${scaffolder}/unplaced_contigs.fasta", "unplaced_contigs.fasta" )
              or die "Error creating symlink: $!";
        }
    }
    return ($linkage_evidence);

}

######################################################################
#
# build_agp
#
# Creates an AGP file from the scaffolds, while generating new
# contig/scaffold outputs meeting ENA requirements (no consecutive runs
# of >=10 N, minimum contig size of 200 bp) if running in 'submission'
# mode, otherwise leaves short contigs and gaps<100bp intact. The
# scaffold_type argument is used to determine the linkage evidence type
# for scaffold gaps
#
# required parameters: $ (tmpdir)
#		     : $ (organism description)
#                    : $ (mode - submission or draft)
#                    : $ (scaffold_type: align or mate_pair)
#
# returns            : $ (0)
#
######################################################################

sub build_agp {

    my $tmpdir   = shift;
    my $organism = shift;
    my $mode     = shift;
    my $evidence = shift;

    die "Unknown evidence type: $evidence"
      unless (    $evidence eq 'paired-ends'
               || $evidence eq 'align_genus'
               || $evidence eq 'align_xgenus' );

    message("Creating AGP file...");

    mkdir "$tmpdir/agp" or die "Error creating agp dir: $!";
    chdir "$tmpdir/agp" or die "Error chdiring to $tmpdir/agp: $!";

    open AGP, ">scaffolds.agp"
      or die "Error opening scaffolds.agp for writing: $! ";
    print AGP "##agp-version 2.0\n";
    print AGP "#$organism\n";

    my $scaffold_inIO  = Bio::SeqIO->new( -file => "$tmpdir/scaffolds.fasta", -format => 'fasta' );
    my $scaffold_outIO = Bio::SeqIO->new( -file => ">scaffolds.fasta",        -format => "fasta" );
    my $contig_outIO   = Bio::SeqIO->new( -file => ">contigs.fasta",          -format => 'fasta' );

    my $contig_count = 0;
    my %gaps;    #overall per-scaffold gaps to return....

    while ( my $scaffold = $scaffold_inIO->next_seq() ) {

        my $contig_start = 0;
        my $scaffold_loc = 1;
        my $scaffold_id  = $scaffold->display_id();
        my $scaff_count  = $1 if ( $scaffold_id =~ /([0-9]+)$/ );
        $scaffold_id = sprintf( "scaffold_%06s", $scaff_count );
        my $scaffold_end = $scaffold->length();
        my ( $contig_end, $gap_start, $gap_end );

        my @bases = split( //, $scaffold->seq() );
        my ( %contigs, @gaps );    #location tracking for included contigs/gaps

        # this is kind of crude, but seems to work...
      BASE: for ( my $b_count = 0 ; $b_count <= $#bases ; $b_count++ ) {
            if ( ( $bases[$b_count] ne 'N' ) && ( !$gap_start ) ) {
                $contig_end = $b_count;
                next BASE;
            }
            elsif ( ( $bases[$b_count] eq 'N' ) & !($gap_start) ) {
                $gap_start = $b_count if ( !$gap_start );
                next BASE;
            }
            elsif ( ( $bases[$b_count] ne 'N' ) && ($gap_start) ) {

                #we have left the gap
                $gap_end = $b_count;
                my $gap_length = $gap_end - $gap_start;
                if ( ( $gap_length < 10 ) && ( $mode eq 'submission' ) ) {

                    # skip gaps <10 bp which are acceptable by EMBL
                    $gap_start = undef;
                    next BASE;
                }
                elsif ( ( $mode eq 'draft' ) && ( $gap_length <= 1 ) ) {

                    # we can leave single ambiguous bases alone
                    $gap_start = undef;
                    next BASE;
                }
                else {
                    $gap_end = $b_count;

                    # Output contigs only > 200 bp
                    if (    ( ( $contig_end - $contig_start ) < 200 )
                         && ( $mode eq 'submission' ) )
                    {

                        #need to extend previous gap to new position including short contig
                        if ( scalar(@gaps) ) {
                            my $last_gap = pop @gaps;
                            my ( $last_start, $last_end ) = split( /-/, $last_gap );
                            my $new_gap = $last_start . '-' . $gap_end;
                            push @gaps, $new_gap;
                            $gap_start    = undef;
                            $contig_start = $b_count;
                            next BASE;
                        }
                    }
                    else {
                        my $contig_id = sprintf( "contig_%06s", ++$contig_count );
                        my $contig_seq =
                          Bio::Seq->new( -display_id => $contig_id,
                                         -seq        => join( '', @bases[ $contig_start .. $contig_end ] ) );
                        $contig_outIO->write_seq($contig_seq);
                        $contigs{$contig_id} = {
                                                 'coords' => $contig_start . '-' . $contig_end,
                                                 length   => $contig_seq->length()
                                               };
                        my $gap = $gap_start . '-' . $gap_end;
                        push @gaps, $gap;

                        $gap_start    = undef;
                        $contig_start = $b_count;
                    }
                }
            }
        }

        # Output last contig from last contig_start position to scaffold end if it is longer than
        # 200 bp and we are running in submission mode, otherwise remove the last gap to truncate the scaffold...
        if (    ( ( $scaffold_end - $contig_start ) > 200 )
             || ( $mode eq 'draft' )
             || $contig_count == 0 )
        {
            my $contig_id = sprintf( "contig_%06s", ++$contig_count );
            my $contig_seq =
              Bio::Seq->new( -display_id => $contig_id,
                             -seq        => join( '', @bases[ $contig_start .. $#bases ] ) );
            $contig_outIO->write_seq($contig_seq);
            $contigs{$contig_id} = {
                                     'coords' => $contig_start . '-' . $scaffold_end,
                                     length   => $contig_seq->length()
                                   };
        }
        else {
            pop @gaps;
        }

        $gaps{$scaffold_id} = \@gaps;

        my $scaffold_part = 0;

        #write AGP output and new scaffolds fasta file
        unlink("contigs.fasta.index") if ( -e "contigs.fasta.index" );
        my $contig_db = Bio::DB::Fasta->new("contigs.fasta");
        my $scaffold_seq;

        my @contig_ids = map { $_->[0] }
          sort { $a->[1] <=> $b->[1] }
          map { [ $_, /(\d+)$/ ] } keys(%contigs);

        if ( $#contig_ids > -1 ) {
            for ( my $i = 0 ; $i <= $#contig_ids ; $i++ ) {
                my $contig_id   = $contig_ids[$i];
                my $contig_data = $contigs{$contig_id};
                my $coords      = $contig_data->{'coords'};
                my $length      = $contig_data->{'length'};

                if ( $contig_id ne "" ) {
                    my $contig = $contig_db->get_Seq_by_id($contig_id);

                    my ( $contig_start, $contig_end ) = split( /-/, $coords );

                    print AGP "$scaffold_id\t"
                      . ( $contig_start + 1 ) . "\t"
                      . ( $contig_end + 1 ) . "\t"
                      . ++$scaffold_part
                      . "\tW\t$contig_id\t1\t$length\t+\n";
                    $scaffold_seq .= $contig->seq();
                    if ( $i < $#contig_ids ) {

                        my $gap = $gaps[$i];
                        my ( $gap_start, $gap_end ) = split( /-/, $gap );
                        my $gap_size = $gap_end - $gap_start;

                        print AGP "$scaffold_id\t"
                          . ( $gap_start + 1 )
                          . "\t$gap_end\t"
                          . ++$scaffold_part
                          . "\tN\t$gap_size\tscaffold\tyes\t$evidence\n";
                        $scaffold_seq .= 'N' x $gap_size;
                    }
                }
            }
        }
        if ($scaffold_seq) {
            my $scaffold_seqobj = Bio::Seq->new( -display_id => $scaffold_id, -seq => $scaffold_seq );
            $scaffold_outIO->write_seq($scaffold_seqobj);
        }
    }

    close AGP     or warn "Error closring $tmpdir/scaffolds.agp: $!";
    chdir $tmpdir or die "Error chdiring to $tmpdir: $!";

    unlink "scaffolds.fasta" or die "Error removing scaffolds.fasta: $!";
    unlink "contigs.fasta"   or die "Error removing contigs.fasta; $!";

    symlink( "agp/scaffolds.agp", "scaffolds.agp" )
      or die "Error creating scaffolds.agp symlink: $!";
    symlink( "agp/scaffolds.fasta", "scaffolds.fasta" )
      or die "Error creating scaffolds.fasta symlink: $!";
    symlink( "agp/contigs.fasta", "contigs.fasta" )
      or die "Error creating contigs.fasta symlink: $!";

    return ( \%gaps );

}

######################################################################
#
# run_prokka
#
# generates annotation on the assembly using prokka
#
# required parameters: $ (tmpdir)
#		       $ (genus)
#                      $ (species)
#                      $ (strain)
#                      $ (locustag)
#                      $ (centre)
#
# returns            : $ (0)
#
######################################################################

sub run_prokka {

    my $tmpdir   = shift;
    my $genus    = shift;
    my $species  = shift;
    my $strain   = shift;
    my $locustag = shift;
    my $centre   = shift;

    message("Starting PROKKA...");

    #use scaffolds if we have them, otherwise contigs....
    my $seqs;
    if ( -e "$tmpdir/scaffolds.fasta" ) {
        $seqs = "$tmpdir/scaffolds.fasta";
    }
    else {
        $seqs = "$tmpdir/contigs.fasta";
    }

    my $type = $1 if ( $seqs =~ /\/(contigs|scaffolds).fasta/ );

    my $cmd = "prokka --addgenes --outdir $tmpdir/prokka --prefix prokka ";
    $cmd .= "--genus $genus " if ( $genus && ( $genus ne "unknown_genus" ) );
    $cmd .= "--species $species "
      if ( $species && ( $species ne "unknown_species" ) );
    $cmd .= "--strain $strain "
      if ( $strain && ( $strain ne "unknown_strain" ) );
    $cmd .= "--locustag $locustag " if ($locustag);
    $cmd .= "--centre $centre "     if ($centre);
    $cmd .= " $seqs ";
    $cmd .= " >prokka.log 2>&1";

    system($cmd);
    if ( $? != 0 ) { print "prokka exited with $?...\n" }

    # my $inIO        = Bio::SeqIO->new( -file => "$tmpdir/prokka/prokka.gbf",     -format => 'genbank' );
    my $inIO        = Bio::SeqIO->new( -file => "$tmpdir/prokka/prokka.gbk",     -format => 'genbank' );
    my $embl_outIO  = Bio::SeqIO->new( -file => ">$tmpdir/prokka/prokka.embl",   -format => 'embl' );
    my $fasta_outIO = Bio::SeqIO->new( -file => ">$tmpdir/prokka/${type}.fasta", -format => 'fasta' );

    # something is losing the scaffold/contig naming so needs to be regenerated...
    # at least ordering should be conserved
    my $count = 0;
    while ( my $seq = $inIO->next_seq() ) {
        my $label;
        ( $type eq 'scaffolds' ) ? ( $label = 'scaffold' ) : ( $label = 'contig' );
        my $id = sprintf( "${label}_%06s", ++$count );
        $seq->display_id($id);
        $seq->accession($id);
        $embl_outIO->write_seq($seq);
        $fasta_outIO->write_seq($seq);
    }

    chdir $tmpdir or die "Error chdiring to $tmpdir: $!";
    unlink "scaffolds.fasta" or die "Error removing scaffolds.fasta: $!" if ( $type eq 'scaffolds' );

    symlink( "prokka/prokka.embl", "$type.embl" )
      or die "Error creating $type.embl symlink: $!";
    symlink( "prokka/scaffolds.fasta", "scaffolds.fasta" )
      or die "Error creating scaffolds.fasta symlink: $!"
      if ( $type eq 'scaffolds' );

    return (0);
}

######################################################################
#
# amosvalidate
#
# Uses abyss's abyss-samtoafg script to convert spades sam and contigs
# to an amos bank
#
# required params: $ (tmpdir)
#                  $ (insert size)
#                  $ (insert size stddev)
#
# returns        : $ (0)
#
######################################################################

sub amosvalidate {

    my $tmpdir        = shift;
    my $insert_size   = shift;
    my $insert_stddev = shift;

    mkdir "$tmpdir/amos" or die "Error creating $tmpdir/amos: $!";
    chdir "$tmpdir/amos" or die "Error running chdir $tmpdir/amos: $!";

    my $seq_file;
    ( -e "$tmpdir/scaffolds.fasta" ) ? ( $seq_file = "scaffolds.fasta" ) : ( $seq_file = "contigs.fasta" );

    open SEQS, "$tmpdir/$seq_file"
      or die "Error opening $tmpdir/$seq_file: $!";
    open OUT, ">$tmpdir/amos/$seq_file"
      or die "Error opening $tmpdir/amos/contigs.fasta: $!";

    my $first_line = 0;    #kludgetastic...
    while ( my $line = <SEQS> ) {
        if ( $line =~ /^>/ ) {
            if ( $first_line > 0 ) {
                print OUT "\n$line";
            }
            else {
                print OUT "$line";
                $first_line++;
            }
        }
        else {
            chomp $line;
            print OUT $line;
        }
    }
    close SEQS;
    close OUT;

    message("Converting to amos bank...");
    my $cmd = $config->{'sam2afg'};
    $cmd .= " -m $insert_size -s $insert_stddev " if ( $insert_size && $insert_stddev );
    $cmd .= " $tmpdir/amos/${seq_file} $tmpdir/bwa/${seq_file}.sam";
    $cmd .= " > amos.afg";

    system($cmd) == 0 or die "Error executing $cmd: $!";

    $cmd = $config->{'amos_dir'}
      . "bank-transact -cb $tmpdir/amos/assembly.bnk -m $tmpdir/amos/amos.afg  > $tmpdir/amos/bank-transact.log 2>&1";
    system($cmd) == 0 or die "Error executing $cmd: $!";

    message("Running amosvalidate");

    # read amosvalidate script and comment out '4xx' lines, which run SNP checks and are extreeeeemly slow....
    open AMOSVALIDATE, $config->{'amos_dir'} . "/amosvalidate"
      or die "Error opening " . $config->{'amos_dir'} . "/amosvalidate: $!";
    open SCRIPT, ">$tmpdir/amos/amosvalidate" or die "Error opening $tmpdir/amos/amosvalidate: $!";
    while ( my $line = <AMOSVALIDATE> ) {
        $line =~ s/^4/#4/;
        print SCRIPT $line;
    }
    close AMOSVALIDATE;
    close SCRIPT;
    chmod 0755, "$tmpdir/amos/amosvalidate" or die "Error running chmod $tmpdir/amos/amosvalidate: $!";

    $cmd = "$tmpdir/amos/amosvalidate $tmpdir/amos/assembly.bnk > $tmpdir/amos/amosvalidate.log 2>&1";
    system($cmd) == 0 or die "Error executing $cmd: $!";

    chdir $tmpdir or die "Error chainging to $tmpdir: $!";
    return (0);

}

######################################################################
#
# find_origin
#
# Attempts to identify location of origin based upon contig overlapping
# base 1 of the reference sequence. This assumes  each reference sequence
# is a complete circular molecular i.e.a chromosome or a plasmid
#
# required params: $ (tmpdir)
#                  $ (scaffolder)
#                  $ (scaffolder_args)
#                  $ (reference)
#                  $ (insert_size)
#                  $ (insert_stddev)
#                  $ (mean_read_length)
#
# returns: $ (0)
#
######################################################################

sub find_origin {

    my $tmpdir           = shift;
    my $scaffolder       = shift;
    my $scaffolder_args  = shift;
    my $reference        = shift;
    my $insert_size      = shift;
    my $stddev           = shift;
    my $mean_read_length = shift;
    my $threads          = shift;

    my $ori_dir = $tmpdir . "/origin";
    mkdir "$ori_dir" or die "Error creating $ori_dir: $!";
    chdir "$ori_dir" or die "Error changing to $ori_dir: $!";
    symlink( "$tmpdir/read1.fastq", "read1.fastq" ) or die "Error creating read1.fastq symlink: $!";
    symlink( "$tmpdir/read2.fastq", "read2.fastq" ) or die "Error creating read1.fastq symlink: $!";

    message("Attempting to identify origin...");

    my $cmd = $config->{'mummer_dir'}
      . "/nucmer $tmpdir/reference_parsed_ids.fasta $tmpdir/scaffolds.fasta -p $ori_dir/ori > $ori_dir/nucmer.log 2>&1";
    system($cmd) == 0 or die "Error executing $cmd: $!";
    $cmd =
      $config->{'mummer_dir'} . "/delta-filter -1 $ori_dir/ori.delta 2>$ori_dir/delta-filter.log > $ori_dir/ori.filter";
    system($cmd) == 0 or die "Error executing $cmd: $!";
    $cmd =
      $config->{'mummer_dir'} . "/show-coords -H $ori_dir/ori.filter 2>$ori_dir/show-coords.log > $ori_dir/ori.coords";
    system($cmd) == 0 or die "Error executing $cmd: $!";

    open COORDS, "$ori_dir/ori.coords"
      or die "Error opening ori.coords: $!";

    my $origin;
  LINE: while ( my $line = <COORDS> ) {
        chomp $line;
        $line =~ s/\|//g;
        $line =~ s/^ *//;
        my @fields = split( /\s+/, $line );
        if ( !$origin && $fields[0] == 1 ) {
            $origin = "$fields[8]:$fields[2]";
            print "Potential origin found at $origin...\n\n";
        }
    }
    close COORDS;

    if ($origin) {
        my $io      = Bio::SeqIO->new( -file => "$tmpdir/scaffolds.fasta",    -format => 'fasta' );
        my $outIO   = Bio::SeqIO->new( -file => '>scaffolds_ori_split.fasta', -format => 'fasta' );
        my $splitIO = Bio::SeqIO->new( -file => ">split_ori_scaff.fasta",     -format => 'fasta' );
      SCAFFOLD: while ( my $scaffold = $io->next_seq() ) {
            my ( $ori_scaffold, $pos ) = split( /:/, $origin );
            if ( $ori_scaffold eq $scaffold->display_id() ) {
                my $part_a = $scaffold->subseq( 1, $pos );
                my $part_b = $scaffold->subseq( $pos + 1, $scaffold->length() );
                my $ori_a = Bio::Seq->new( -display_id => $scaffold->display_id . "_A", -seq => $part_a );
                my $ori_b = Bio::Seq->new( -display_id => $scaffold->display_id . "_B", -seq => $part_b );

                $splitIO->write_seq($ori_a);
                $splitIO->write_seq($ori_b);
                undef($splitIO);

                # Now rerun scaffolder on split sequence containing origin.
                run_scaffolder(
                                $ori_dir,          "$tmpdir/reference_parsed_ids.fasta",
                                $scaffolder,       $scaffolder_args,
                                $insert_size,      $stddev,
                                2,                 "$tmpdir/origin/split_ori_scaff.fasta",
                                $mean_read_length, $threads
                              );
                my $scaffIO =
                  Bio::SeqIO->new( -format => 'fasta', -file => "$ori_dir/${scaffolder}_2/scaffolds.fasta" );
                while ( my $scaffold = $scaffIO->next_seq() ) {
                    $outIO->write_seq($scaffold);
                }
                next SCAFFOLD;
            }
            else {
                $outIO->write_seq($scaffold);
            }
        }

        #renumber scaffolds to ensure they are unique...
        my $count = 0;
        my $inIO = Bio::SeqIO->new( -file => "$tmpdir/origin/scaffolds_ori_split.fasta", -format => "fasta" );
        $outIO = Bio::SeqIO->new( -file => ">$tmpdir/origin/scaffolds_renumbered.fasta", -format => "fasta" );

        while ( my $seq = $inIO->next_seq() ) {
            $seq->display_id( "scaffold_" . ++$count );
            $outIO->write_seq($seq);
        }

        chdir $tmpdir                     or die "Error changing to $tmpdir: $!";
        unlink("$tmpdir/scaffolds.fasta") or die "Error unlinking $tmpdir/scaffolds.fasta:$!";
        symlink( "$tmpdir/origin/scaffolds_renumbered.fasta", "$tmpdir/scaffolds.fasta" )
          or die "Error symlinking scaffolds_ori_split.fasta:$!";
    }

    return (0);
}

######################################################################
#
# order_scaffolds
#
# Identifies origin based on homology with reference.
# Resulting scaffolds are then ordered and oriented relative to the reference...
#
# required params: $ (tmpdir)
#                  $ (fasta reference)
#
# returns        : $ (0)
#
######################################################################

sub order_scaffolds {

    my $tmpdir    = shift;
    my $reference = shift;

    mkdir "$tmpdir/orientating" or die "Error creating $tmpdir/orientating: $!";
    chdir "$tmpdir/orientating"
      or die "Error changing to $tmpdir/orientating: $!";

    message("Orienting scaffolds vs. reference...");

    my $cmd = $config->{'mummer_dir'}
      . "/nucmer $tmpdir/$reference $tmpdir/scaffolds.fasta -p $tmpdir/orientating/ori2 > $tmpdir/orientating/nucmer2.log 2>&1";
    system($cmd) == 0 or die "Error executing $cmd: $!";
    $cmd = $config->{'mummer_dir'}
      . "/delta-filter -1 $tmpdir/orientating/ori2.delta 2>$tmpdir/orientating/delta-filter2.log > $tmpdir/orientating/ori2.filter";
    system($cmd) == 0 or die "Error executing $cmd: $!";
    $cmd = $config->{'mummer_dir'}
      . "/show-coords -H $tmpdir/orientating/ori2.filter 2>$tmpdir/orientating/show-coords2.log > $tmpdir/orientating/ori2.coords";
    system($cmd) == 0 or die "Error executing $cmd: $!";

    open COORDS, "$tmpdir/orientating/ori2.coords"
      or die "Error opening ori.coords: $!";
    my ( %orientations, $start, $end, $orient );
    my %orient_count = ( '+' => 0, '-' => 0 );
    my $contig = '';
  LINE: while ( my $line = <COORDS> ) {
        chomp $line;
        $line =~ s/\|//g;
        $line =~ s/^ *//;
        my @fields = split( /\s+/, $line );

        if ( ( $contig ne $fields[8] ) ) {
            if ( $contig ne '' ) {    #end of previous contig
                store_orientation( \%orientations, $contig, \%orient_count );
            }
        }
        $contig = $fields[8];
        $start  = $fields[2];
        $end    = $fields[3];
        my $length;
        if ( $start < $end ) {
            $orient = '+';
            $length = $end - $start;
        }
        else {
            $orient = '-';
            $length = $start - $end;
        }
        if ( $orient_count{$orient} ) {
            $orient_count{$orient} = $orient_count{$orient} + $length;
        }
        else {
            $orient_count{$orient} = $length;
        }
    }

    close COORDS;

    # record data for last contig...
    store_orientation( \%orientations, $contig, \%orient_count );

    my $orig_length     = 0;
    my $oriented_length = 0;    #track how much sequence we align ok...
    my $unplaced_length = 0;
    my @unplaced;

    my $io    = Bio::SeqIO->new( -file => '../scaffolds.fasta', -format => 'fasta' );
    my $outIO = Bio::SeqIO->new( -file => '>scaffolds.fasta',   -format => 'fasta' );

    # Reorientate contigs and break origin, rewriting to a new file...
    my $i = 0;                  #for renumbering scaffolds...
    while ( my $scaffold = $io->next_seq() ) {
        $orig_length += $scaffold->length();

        my $id = "scaffold_" . ++$i;
        if ( $orientations{ $scaffold->display_id() } ) {
            if ( $orientations{ $scaffold->display_id() } eq '+' ) {
                $scaffold->display_id($id);
                $outIO->write_seq($scaffold);
                $oriented_length += $scaffold->length();
            }
            else {
                my $rev = $scaffold->revcom();
                $rev->display_id($id);
                $outIO->write_seq($rev);
                $oriented_length += $scaffold->length();
            }
        }
        else {
            push @unplaced, $scaffold;
        }

    }

    foreach my $scaffold (@unplaced) {
        $outIO->write_seq($scaffold);
        $unplaced_length += $scaffold->length();
    }

    my $ref_length;
    my $refIO = Bio::SeqIO->new( -format => 'fasta', -file => "$tmpdir/$reference" );
    while ( my $seq = $refIO->next_seq() ) {
        $ref_length += $seq->length();
    }

    my $tb = Text::ASCIITable->new();
    $tb->setCols( "", "Length (bp)" );
    $tb->addRow( "Reference Sequence", $ref_length );
    $tb->addRow( "Assembly",           $orig_length );
    $tb->addRow( "Orientated contigs", $oriented_length );
    $tb->addRow( "Unaligned contigs",  $unplaced_length );

    print $tb, "\n";

    chdir $tmpdir or die "Error changing to $tmpdir: $!";
    unlink("scaffolds.fasta")
      or die "Error removing scaffolds.fasta symlink";
    symlink( "orientating/scaffolds.fasta", "scaffolds.fasta" )
      or die "Error linking orientating/scaffolds.fasta: $!";

    return (0);

}

######################################################################
#
# store_orientation
#
# save correct contig orientation in hash passed as first parameter
#
# since alignment fof scaffold can contain blocks in reverse orientation
# need to use the 'prevailing' orientation based on number of +/- blocks
#
# required parameters: $ (orientations hash)
#                      $ (contig)
#                      $ (hash of +/- base counts)
#
# returns           : none
#
######################################################################

sub store_orientation {

    my $orientations = shift;
    my $contig       = shift;
    my $orient_count = shift;

    if (
         (
              ( defined( $orient_count->{'+'} ) && defined( $orient_count->{'-'} ) )
           && ( $orient_count->{'+'} > $orient_count->{'-'} )
         )
         || ( defined( $orient_count->{'+'} ) && !defined( $orient_count->{'-'} ) )
       )
    {
        $orientations->{$contig} = '+';
    }
    else {
        $orientations->{$contig} = '-';
    }

}

######################################################################
#
# run_varcaller
#
# Carries out variant calling using requested variant caller
#
# required params: $ (tmpdir)
#                  $ (varcall)
#                  $ (no. threads)
#                  $ (read length)
#
# returns        : $ (0)
#
######################################################################

sub run_varcaller {

    my $tmpdir      = shift;
    my $varcall     = shift;
    my $threads     = shift;
    my $read_length = shift;

    message("Running variant calling ($varcall)...");

    my ( $cmd, $caller_cmd, $create_dir );
    my $varcallers = $config->{'varcallers'};
    foreach my $caller (@$varcallers) {
        my $name = $caller->{'name'};
        if ( $name eq $varcall ) {
            $caller_cmd = $caller->{'command'};
            $create_dir = $caller->{'create_dir'};
        }
    }
    my $vardir = "$tmpdir/var_${varcall}/";
    if ($create_dir) {
        mkdir "$tmpdir/var_${varcall}" or die "Error creating $tmpdir/var_${varcall}: $! ";
        chdir "$tmpdir/var_${varcall}" or die "Error chdiring to $tmpdir/var_${varcall}: $!";
    }

    symlink( "$tmpdir/reference_parsed_ids.fasta", "$vardir/reference.fasta" )
      or die "Error creating $vardir/reference.fasta symlink: $!";
    $cmd = $config->{'bwa_dir'} . symlink( "$tmpdir/read1.fastq", "$vardir/read1.fastq" )
      or die "Error creating symlink: $! ";
    symlink( "$tmpdir/read2.fastq", "$vardir/read2.fastq" )
      if ( -e "$tmpdir/read2.fastq" )
      or die "Error creating symlink: $! ";

    print "BWA aligning reads to assembly...\n";
    my $samtools_dir = $config->{'samtools_dir'};

    $cmd = $config->{'bwa_dir'} . "/bwa index $vardir/reference.fasta >$vardir/bwa_index.log 2>&1";
    system($cmd) == 0 or die " Error running $cmd";

    # Use bwa-bwt for 'short' reads less than 100 bp, and bwa-mem for longer reads
    if ( $read_length <= 100 ) {
        $cmd =
            $config->{'bwa_dir'}
          . "/bwa aln -t $threads $vardir/reference.fasta $vardir/read1.fastq > $vardir/read1.sai"
          . " 2> $vardir/bwa_sai1.log";
        system($cmd) == 0 or die "Error running $cmd";
        if ( -e "$vardir/read2.fastq" ) {
            $cmd =
                $config->{'bwa_dir'}
              . "/bwa aln -t $threads $vardir/reference.fasta $vardir/read2.fastq > $vardir/read2.sai"
              . " 2> $vardir/bwa_sai2.log";
            system($cmd) == 0 or die "Error running $cmd";
            $cmd =
                $config->{'bwa_dir'}
              . "/bwa sampe $vardir/reference.fasta $vardir/read1.sai $vardir/read2.sai "
              . "$vardir/read1.fastq $vardir/read2.fastq";
            $cmd .= " 2> $vardir/sampe.log > $vardir/scaffolds.sam";
            system($cmd) == 0 or die "Error running $cmd";
        }
        else {
            $cmd = $config->{'bwa_dir'} . "/bwa samse $vardir/reference.fasta $vardir/read1.sai $vardir/read1.fastq";
            $cmd .= "2> $vardir/samse.log > $vardir/scaffolds.sam";
            system($cmd) == 0 or die "Error running $cmd";
        }
    }
    else {
        if ( !-e "$vardir/read2.fastq" ) {

            # single-ended long reads
            $cmd =
                $config->{'bwa_dir'}
              . "/bwa mem -t $threads -M $vardir/reference.fasta $vardir/read1.fastq > reference.sam "
              . "2>$vardir/bwa_mem.log";
            system($cmd) == 0 or die "Error running $cmd: $!";
        }
        else {

            # paired-end long reads
            $cmd =
                $config->{'bwa_dir'}
              . "/bwa mem -t $threads -M $vardir/reference.fasta $vardir/read1.fastq $vardir/read2.fastq >reference.sam "
              . "2>$vardir/bwa_mem.log";
            system($cmd) == 0 or die "Error running $cmd: $!";
        }
    }

    $cmd = "$samtools_dir/samtools view -q 10 -Sb $vardir/reference.sam 2>$vardir/samtoolsview.log"
      . "|$samtools_dir/samtools sort - $vardir/reference";
    system($cmd) == 0 or die "Error running $cmd";

    $cmd = "$samtools_dir/samtools index $vardir/reference.bam 2>$vardir/samtools_index.log";
    system($cmd) == 0 or die "Error running $cmd";

    $caller_cmd =~ s/__BUGBUILDER_BIN__/$FindBin::Bin/;
    $caller_cmd =~ s/__TMPDIR__/$tmpdir/;
    $caller_cmd =~ s/__THREADS__/$threads/;

    system($caller_cmd) == 0 or die "Error running $cmd: $!";
    chdir $tmpdir            or die " Error chdiring to $tmpdir: $! ";

    symlink( "var_${varcall}/var.filtered.vcf", "reference.variants.vcf" )
      or die "Error creating $tmpdir/reference.variants.vcf: $!";
    my $varcount = `grep -vc ^# $tmpdir/reference.variants.vcf`;
    chomp $varcount;

    print "\nIdentified $varcount variants...\n";

    return (0);
}

######################################################################
#
# build_comparisons
#
# generates a comparison appropriate for viewing with ACT and a
# MUMmerplot to provide a quick overview
#
# required params: $ (tmpdir)
#                  $ (reference)
#                  $ (organism)
#
# returns        : $ (0)
#
######################################################################

sub build_comparisons {

    my $tmpdir    = shift;
    my $reference = shift;

    mkdir "$tmpdir/comparisons" or die "Error creating $tmpdir/comparisons: $!";
    chdir "$tmpdir/comparisons"
      or die "Error changing to $tmpdir/comparisons: $!";
    my $io     = Bio::SeqIO->new( -format => 'fasta', -file => "$tmpdir/$reference" );
    my $ref    = $io->next_seq();
    my $ref_id = $ref->display_id();
    $ref_id = ( split( /\|/, $ref_id ) )[0] if ( $ref_id =~ /\|/ );

    my $query;
    if ( -l "$tmpdir/scaffolds.fasta" ) {
        $query = "$tmpdir/scaffolds.fasta";
    }
    else {
        $query = "$tmpdir/contigs.fasta";
    }
    my $cmd =
        $config->{'blast_dir'}
      . "/blastn  -query $query -subject $tmpdir/$reference "
      . "-out $tmpdir/comparisons/comparison_vs_$ref_id.blastout -outfmt 6 > $tmpdir/comparisons/blast.log";
    system($cmd) == 0 or die "Error running $cmd";

    # also build a mummerplot in png format...
    $cmd = $config->{'mummer_dir'} . "nucmer --prefix $ref_id $tmpdir/$reference $query >nucmer.log 2>&1";
    system($cmd) == 0 or die "Error running $cmd";
    $cmd = $config->{'mummer_dir'}
      . "mummerplot -large --filter --layout -p $ref_id -t png -R $tmpdir/$reference -Q $query $ref_id.delta  >mummerplot.log 2>&1";
    system($cmd) == 0 or die "Error running $cmd";

    $cmd = $config->{'mummer_dir'}
      . "mummerplot -large --filter --layout -c -p ${ref_id}_pip -t png -R $tmpdir/$reference -Q $query $ref_id.delta  >mummerplot_pip.log 2>&1";
    system($cmd) == 0 or die "Error running $cmd";

    chdir $tmpdir or die "Error chdiring to $tmpdir: $!";
    symlink( "comparisons/comparison_vs_$ref_id.blastout", "comparison_vs_$ref_id.blastout" )
      or die "Error creating symlink: $!";
    symlink( "comparisons/$ref_id.png", "comparison_vs_$ref_id.png" )
      or die "Error creating symlink: $!";
    symlink( "comparisons/${ref_id}_pip.png", "comparison_vs_${ref_id}_pip.png" )
      or die "Error creating symlink: $!";

    return (0);

}

######################################################################
#
# get_contig_to_iid_mapping
#
# generates a mapping of contig ids to amos IID
#
# required params: $ (tmpdir)
#
# returns        : $ (0)
#
######################################################################

sub get_contig_to_iid_mapping {

    my $tmpdir = shift;
    message("Extracting contig -> AMOS iid mapping...");

    my $amos_dir = $config->{'amos_dir'};
    my $mapping =
`$amos_dir/bank-report -i -p -b $tmpdir/amos/assembly.bnk CTG 2> /dev/null|cut -f2,3 > $tmpdir/amos/ctg_to_iid.txt`;

    return (0);
}

######################################################################
#
# summarise_amosvalidate
#
# postprocesses amosvalidate outputs to make more readily digestible
#
# required params: $ (tmpdir)
#
# returns        : $ ($ - hashref to parsed results)
#
######################################################################

sub summarise_amosvalidate {

    my $tmpdir = shift;

    message("processing amosvalidate results...");

    my $amosvalidate = $tmpdir . "/amos";
    my $iid_mapping  = $tmpdir . "/amos/ctg_to_iid.txt";
    my %results;

    # Read Contig->iid mapping into a hash
    my %contig_to_iid;
    open IID, $iid_mapping or die "Error opening $iid_mapping: $!";
    while (<IID>) {
        my ( $contig, $iid ) = split(/\s+/);
        chomp $iid;
        $contig_to_iid{$iid} = $contig;
        $results{$contig}    = [];
    }
    close IID;

    opendir AMOSVALIDATE, $amosvalidate
      or die "Error opening $amosvalidate: $!";
    my @outputs = grep /feat$/, readdir AMOSVALIDATE;
    close AMOSVALIDATE;

    open ALL, "$amosvalidate/assembly.all.feat" or die "Could not open $amosvalidate/assembly.all.feat: $!";
    while (<ALL>) {
        my @fields    = split(/\s+/);
        my $contig_id = $contig_to_iid{ $fields[0] };
        my $start     = $fields[2];
        my $end       = $fields[3];
        my $type      = $fields[4];
        my $res_arr   = $results{$contig_id};
        $start = 0 if ( $start < 0 );    #no idea how it ends up with a negative start...but it does...
        push @$res_arr, { 'start' => $start, 'end' => $end, type => $type, };
        $results{$contig_id} = $res_arr;
    }
    close ALL;
    return ( \%results );
}
######################################################################
#
# merge_annotations
#
# updates annotated generated embl file with amosvalidate results
# and gap locations
#
# required parameters: $ (tmpdir)
#                      $ (hashref of amosvalidate results, keyed on contigid)
#                      $ (hashref of scaffold gaps)
#                      $ (genus)
#                      $ (species)
#                      $ (strain)
#
# returns            : $ (none)
#
######################################################################

sub merge_annotations {

    my $tmpdir               = shift;
    my $amosvalidate_results = shift;
    my $gaps                 = shift;
    my $genus                = shift;
    my $species              = shift;
    my $strain               = shift;

    message("Merging annotations");

    mkdir "$tmpdir/annotation_merge"
      or die "Error creating $tmpdir/annotation_merge: $!";
    chdir "$tmpdir/annotation_merge"
      or die "Error chdiring to $tmpdir/annotation_merge: $!";
    my $filename;

    ( -e "$tmpdir/scaffolds.embl" ) ? ( $filename = "scaffolds.embl" ) : ( $filename = "contigs.embl" );

    my $IO = Bio::SeqIO->new( -format => 'embl',
                              -file   => "$tmpdir/$filename" );

    my $outIO =
      Bio::SeqIO->new( -format => 'embl',
                       -file   => ">$filename" );

    my %amos_colours = (
                         'CE_STRETCH'      => '0 128 128',
                         'CE_COMPRESS'     => '0 128 128',
                         'HIGH_SNP'        => '128 128 0',
                         'HIGH_READ_CVG'   => '255 0 0',
                         'HIGH_KMER'       => '255 0 0',
                         'KMER_COV'        => '255 0 0',
                         'HIGH_OUTIE_CVG'  => '255 0 0',
                         'HIGH_NORMAL_CVG' => '255 0 0',
                         'LOW_GOOD_CVG'    => '0 0 255',
                       );

    my %amos_notes = (
                       'CE_STRETCH' => 'Stretched mate-pairs: Possible repeat copy number expansion or other insertion',
                       'CE_COMPRESS'     => 'Compressed mate-pairs; Possible collapsed repeat',
                       'HIGH_SNP'        => 'High SNP frequency',
                       'HIGH_READ_CVG'   => 'High read coverage; Possible collapsed repeat',
                       'HIGH_KMER'       => 'High frequency of normalized kmers: Possible collapsed repeat',
                       'KMER_COV'        => 'High frequency of normalized kmers: Possible collapsed repeat',
                       'LOW_GOOD_CVG'    => 'Low coverage',
                       'HIGH_NORMAL_CVG' => 'High coverage',
                       'HIGH_OUTIE_CVG'  => 'High outie coverage',
                     );

    while ( my $embl_record = $IO->next_seq() ) {
        my $orig_id = $embl_record->display_id();
        $embl_record->display_id("$orig_id");
        $embl_record->accession_number("$orig_id");
        $embl_record->division('PRO');
        $embl_record->molecule('genomic DNA');
        $embl_record->is_circular(1);

        $embl_record->add_date( get_embl_date() );
        $embl_record->description("$genus $species $strain genome scaffold");

        my @comments = $embl_record->annotation->get_Annotations('comment');
        my $annot    = new Bio::Annotation::Collection;
        my $comment  = Bio::Annotation::Comment->new;
        $comment->text("Assembled using BugBuilder from http://github.com/jamesabbott/BugBuilder");
        push @comments, $comment;
        foreach my $c (@comments) {
            $annot->add_Annotation( 'comment', $c );
        }
        $embl_record->annotation($annot);

        # remove source entries from feature table
        my @features = $embl_record->get_SeqFeatures();
        $embl_record->flush_SeqFeatures();

        # retrieve amosvalidate results for this contig and sort by start co-ordinate...
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
                }
            }
        }

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

######################################################################
#
# run_cgview
#
# Runs cgview to generate a genome map from the annotations
#
# Required parameters: $ (tmpdir)
#
# Returns: $ (0)
#
######################################################################

sub run_cgview {

    my $tmpdir = shift;

    my $cgview_dir  = $config->{'cgview_dir'};
    my $xml_creator = $cgview_dir . "cgview_xml_builder/cgview_xml_builder.pl";
    my $java        = $config->{'java'};

    message("Creating genome visualisaton...");

    chdir $tmpdir  or die "Error chdiring to $tmpdir: $!";
    mkdir "cgview" or die "Error creating cgview directory: $!";
    chdir "cgview" or die "Error chdiring to cgview $!";

    my ( $embl, $outfile );
    if ( -e "../scaffolds.embl" ) {
        $embl    = "../scaffolds.embl";
        $outfile = "scaffolds_cgview.png";
    }
    else {
        $embl    = "../contigs.embl";
        $outfile = "contigs_cgview.png";
    }

    my $cmd = "$xml_creator -sequence $embl -output scaffolds_cgview.xml -gc_skew T >xml_creator.log 2>&1";
    system("$cmd") == 0 or die "Error running $cmd: $!";

    $cmd = "${java} -jar ${cgview_dir}/cgview.jar -f png -i scaffolds_cgview.xml -o $outfile >cgview.log 2>&1";
    system("$cmd") == 0 or die "Error running $cmd: $!";

    chdir $tmpdir or die "Error chdiring to $tmpdir: $!";
    symlink( "cgview/$outfile", "$outfile" ) or die "Error creating symlink: $!";

    return (0);
}

######################################################################
#
# get_contig_stats
#
# Reports contig statistics on assembly. Reports on scaffolds or
# contigs depending upon 2nd argument passed - contigs gives values
# for all contigs and those >200bp
#
# required params: $ (path to contigs)
#                  $ ('scaffolds'|'contigs')
#
# returns          $ (0)
#
######################################################################

sub get_contig_stats {

    my $file = shift;
    my $type = shift;

    my $IO = Bio::SeqIO->new( -format => 'fasta', -file => $file );
    my ( %lengths, @all_lengths, $count, $tot_length, $max, $progress, $n50, $l50, $n50_tot );
    my (
         %lengths_200,  @all_lengths_200, $count_200, $tot_length_200, $max_200,
         $progress_200, $n50_200,         $l50_200,   $n50_tot_200
       );

    while ( my $seq = $IO->next_seq() ) {
        $lengths{ $seq->length() }++;
        $count++;
        $tot_length += $seq->length();
        push @all_lengths, $seq->length();
        if ( $seq->length() > 200 ) {
            $lengths_200{ $seq->length() }++;
            $count_200++;
            $tot_length_200 += $seq->length();
            push @all_lengths_200, $seq->length();
        }
    }

    my $fifty     = $tot_length / 2;
    my $fifty_200 = $tot_length_200 / 2;

    my @sorted_lengths     = sort { $b <=> $a } @all_lengths;
    my @sorted_lengths_200 = sort { $b <=> $a } @all_lengths_200;

    # l50
    foreach ( sort { $b <=> $a } keys(%lengths) ) {
        $max = $_ if ( !$max );
        $progress += $_ * ( $lengths{$_} );
        if ( $progress >= $fifty ) {
            $l50 = $_;
            last;
        }
    }
    foreach ( sort { $b <=> $a } keys(%lengths_200) ) {
        $max_200 = $_ if ( !$max_200 );
        $progress_200 += $_ * ( $lengths_200{$_} );
        if ( $progress_200 >= $fifty_200 ) {
            $l50_200 = $_;
            last;
        }
    }

    # n50
    foreach my $length (@sorted_lengths) {
        $n50_tot += $length;
        $n50++;

        # $l50 = $length;
        last if ( $n50_tot >= $fifty );
    }
    foreach my $length (@sorted_lengths_200) {
        $n50_tot_200 += $length;
        $n50_200++;

        # $l50_200 = $length;
        last if ( $n50_tot_200 >= $fifty_200 );
    }

    my $tb = Text::ASCIITable->new();
    if ( $type eq 'contigs' ) {
        $type = ucfirst($type);
        $tb->setCols( "", "All $type", "$type >200bp" );
        $tb->addRow( "$type count",   $count,      $count_200 );
        $tb->addRow( "Max Length",    $max,        $max_200 );
        $tb->addRow( "Assembly size", $tot_length, $tot_length_200 );
        $tb->addRow( "L50",           $l50,        $l50_200 );
        $tb->addRow( "N50",           $n50,        $n50_200 );
    }
    else {
        $type = ucfirst($type);
        $tb->setCols( "", "All $type" );
        $tb->addRow( "$type count",   $count );
        $tb->addRow( "Max Length",    $max );
        $tb->addRow( "Assembly size", $tot_length );
        $tb->addRow( "L50",           $l50 );
        $tb->addRow( "N50",           $n50 );
    }
    print $tb . "\n";

    return (0);
}

######################################################################
#
#  Pretty formats a status message
#
#  requried params: $ (message to display)
#
#  returns        : $ (0)
#
######################################################################

sub message {

    my $message = shift;
    print "\n\n", "*" x 80, "\n";
    print "*",    " " x 78, "*\n";
    print "* $message", " " x ( 77 - length($message) ), "*\n";
    print "*", " " x 78, "*\n";
    print "*" x 80, "\n\n";

}

######################################################################
#
# get_embl_date
#
# returns the current date in EMBL style...
#
# required params: none
#
# returns: $ (formatted date)
#
######################################################################

sub get_embl_date {

    my @months = qw( JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC );
    my @vals   = ( localtime() )[ 3 .. 5 ];
    my $date   = $vals[0] . '-' . $months[ $vals[1] ] . '-' . ( 1900 + $vals[2] );

    return ($date);

}

######################################################################
#
# parse_fasta_id
#
# attempt to parse an ID from the reference fasta file. Different fasta formats
# make this tricky, so we will specifically parse plain IDs, ENA and NCBI
# formats
#
# required params: $ (fasta file)
#
# returns: $ (hashref of parsed sequence IDs)=
#
######################################################################

sub parse_fasta_id {

    my $file = shift;
    my $IO = Bio::SeqIO->new( -format => 'fasta', -file => $file );
    my @seq_ids;
    while ( my $seq = $IO->next_seq() ) {
        my $id = $seq->display_id();
        if ( $id =~ /\|/ ) {
            my @fields = split( /\|/, $id );
            if ( $fields[0] eq 'ENA' ) {
                push( @seq_ids, $fields[1] );
            }
            elsif ( $fields[0] eq 'gi' ) {
                push( @seq_ids, $fields[3] );
            }
        }
        else {
            if ( $id =~ />([^ ])/ ) {
                push( @seq_ids, $1 );
            }
        }
    }
    return (@seq_ids);
}

######################################################################
#
# parse_ref_id
#
# Similar to parse_fasta_id above, but works directly on a passed ID
#
# required params: $ (id)
#
# returns: $ (id)
#
######################################################################

sub parse_ref_id {

    my $id = shift;
    my $parsed_id;
    if ( $id =~ /\|/ ) {
        my @fields = split( /\|/, $id );
        if ( $fields[0] eq 'ENA' ) {
            $parsed_id = $fields[1];
        }
        elsif ( $fields[0] eq 'gi' ) {
            $parsed_id = $fields[3];
        }
    }
    else {
        if ( $id =~ />([^ ])/ ) {
            $parsed_id = $1;
        }
        else {
            $parsed_id = $id;
        }
    }
    return ($parsed_id);
}

######################################################################
#
# show_tools
#
# Reports on configured assemblers, scaffolders and platforms
#
# Required params: $ (config hash)
#
# Returns:         $ ()
#
######################################################################

sub show_tools {

    my $config = shift;

    my @available_assemblers  = map { $_->{'name'} } @{ $config->{'assemblers'} };
    my @available_scaffolders = map { $_->{'name'} } @{ $config->{'scaffolders'} };
    my @available_mergers     = map { $_->{'name'} } @{ $config->{'merge_tools'} };
    my @available_finishers   = map { $_->{'name'} } @{ $config->{'finishers'} };
    my @platforms             = map { $_->{'platforms'} } @{ $config->{'assembler_categories'} };

    my %available_platforms;
    foreach my $platform (@platforms) {
        foreach my $p (@$platform) {
            $available_platforms{$p}++;
        }
    }
    my @available_platforms = sort( keys(%available_platforms) );

    print "\nWelcome to BugBuilder\n\n";
    print "Available assemblers: " . join( ", ",                @available_assemblers ),
      "\n" . "Available scaffolders: " . join( ", ",            @available_scaffolders ),
      "\n" . "Available assembly merging tools: " . join( ", ", @available_mergers ),
      "\n" . "Available finishing tools: " . join( ", ",        @available_finishers ),
      "\n" . "Configured platforms: " . join( ", ", @available_platforms ), "\n\n";

    return ();
}
