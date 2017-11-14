#!/usr/bin/env python3
"""
######################################################################
#
# $HeadURL: https://bss-srv4.bioinformatics.ic.ac.uk/svn/BugBuilder/trunk/bin/run_sis $
# $Author: jamesa $
# $Revision: 179 $
# $Date: 2016-03-10 10:32:17 +0000 (Thu, 10 Mar 2016) $
#
# Wrapper for SIS to permit use via BugBuilder scaffolding stages
#
# This file is part of BugBuilder (https://github.com/jamesabbott/BugBuilder)
# and is distributed under the Artistic License 2.0 - see accompanying LICENSE
# file for details
#
######################################################################

=pod

=head1 NAME

run_sis

=head1 SYNOPSIS

run_sis --tmpdir BugBuilder_working_directory --reference reference_genome.fasta [--help]

=head1 DESCRIPTION

Wrapper for SIS scaffolder to permit use withing BugBuilder's scaffolding configuration.
SIS makes use of MUMmer to carry out alignments, and post-processes the show_coords output.
These MUMmer stages need running prior to executing SIS itself. Following SIS execution the
generated scaffolds (which consist of ordered contigs, with one scaffold per fasta file)
are reprocessed into a multifasta file of 'N' gapped scaffold sequences.

=head1 REQUIRED ARGUMEMNTS

=over 4

=item B<tmpdir>: BugBuilder working directory

=item B<scaff_dir>: Directory within tmpdir for scaffolding these sequences...

=item B<reference>: Fasta formatted reference genome for aligning contigs against

=item B<contigs>: Fasta formatted file of contig sequences to scaffold

=back

=head1 OPTIONAL ARGUMENTS

=over 4

=item B<help>: display short help text

=item B<man>: display full documentation

=back

=head1 REPORTING BUGS

Please report any bugs/issues via github:
https://github.com/jamesabbott/BugBuilder/issues/new

=head1 AUTHOR - James Abbott

Email j.abbott@imperial.ac.uk

=cut

use warnings;
use strict;

use FindBin;
use YAML::XS qw(LoadFile);
use Getopt::Long;
use Pod::Usage;
use Carp qw(croak cluck);
use Bio::SeqIO;
use File::Basename;

"""
import glob
import subprocess
import sys
import os
from .shared_methods import make_nucmer_delta_show_cmds
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from scaffoldsis import sis, multifasta

def make_sis_etc_cmds(config, ref, contigs, scaff_dir):
    cmds = make_nucmer_delta_show_cmds(
        config, ref=ref, query=contigs,
        out_dir=scaff_dir, prefix="sis", header=True)
    # sis
    sis_cmd = "{0} {1}/sis.coords > {1}/sis.sis".format(
        config.sis, scaff_dir)
    cmds.append(sis_cmd)
    # multifasta
    multifasta_cmd = "{0} {1}/sis.sis {2} > {1}/unfilled_scaffolds.fna".format(
        config.multifasta, scaff_dir, contigs)
    cmds.append(multifasta_cmd)
    return cmds


def run(config, args, results, ref, contigs, scaff_dir, logger):
    """ we need to explicitly set the reference and the contigs files cause we may have partitioned them (using blast to see which contigs go with which contigs)
    """
    logger.warning("Using SIS to scaffold %s against %s", contigs, ref)
    cmd_list = make_sis_etc_cmds(config, ref=ref, contigs=contigs, scaff_dir=scaff_dir)
    for cmd in cmd_list:
        logger.debug(cmd)
        print(cmd)
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)

    ref_id = os.path.basename(scaff_dir)

    #$ref_id=~s/SIS_//;

    # Build a multifasta file of N-gapped scaffolds
    multi = []
    singletons = []
    scaffolds =  glob.glob(os.path.join(scaff_dir, '*.fna'))
    print(scaffolds)
    for scaf in scaffolds:
        with open(scaf, "r") as inf:
            contig_count = 0
            for rec in SeqIO.parse(inf, "fasta"):
                contig_count = contig_count + 1
            if contig_count > 1:
                multi.append(scaf)
            else:
                singletons.append(scaf)
    print("contigs: %s singltons:%s multi:%s" %\
          (contig_count, singletons, multi))
    with open(os.path.join(scaff_dir, "scaffolds.fasta"), "w") as outf:
        new_seq = ""
        for path in multi:
            with open(path, "r") as inf:
                for rec in SeqIO.parse(inf, "fasta"):
                    new_seq = new_seq + ("N" * 100) + rec.seq
        if new_seq != "":
            SeqIO.write(SeqRecord(new_seq, id="combined_from_sis"),
                        outf, "fasta")
        for spath in singletons:
            with open(spath, "r") as inf:
                for rec in SeqIO.parse(inf, "fasta"):
                    SeqIO.write(rec, outf, "fasta")
