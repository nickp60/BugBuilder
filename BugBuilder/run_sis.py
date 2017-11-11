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
from .shared_methods import make_nucmer_delta_show_cmds


def make_sis_etc_cmds(config, args, contigs, scaff_dir):
    cmds = make_nucmer_delta_show_cmds(
        config, ref=args.reference, query=contigs,
        out_dir=scaff_dir, prefix="sis", header=True)
    # sis
    sis_cmd = "{0} {1}/sis.coords > {1}/sis.sis".format(
        config.sis, scaff_dir)
    cmds.append(sis_cmd)
    # multifasta
    multifasta_cmd = "{0} {1}/sis.sis {2}".format(
        config.multifasta, scaff_dir, contigs)
    cmds.append(multifasta_cmd)
    return cmds


def run_sis(config, args, contigs, scaff_dir):
    cmd_list = make_sis_etc_cmds(config, args, contigs, scaff_dir)

    ref_id = os.path.basename(scaff_dir)

    #$ref_id=~s/SIS_//;

    # Build a multifasta file of N-gapped scaffolds
    multi = []
    singletons = []
    scaffolds =  glob.glob(scaff_dir  + '*.fna')
    for scaf in scaffolds:
        contig_count = 0
        with open(scaf, "r") as inf:
            for rec in SeqIO.parse(inf, fasta):
                contig_count = contig_count + 1
            if contig_count > 1:
                multi = multi.append(scaf)
            else:
                singletons = singletons.append(scaf)


#     my $multiFastaIO = Bio::SeqIO->new( -format => 'fasta',
#                                         -file   => ">$scaff_dir/scaffolds.fasta" );

#     my $i = 0;
#     foreach my $scaffold (@multi) {

#         my $io = Bio::SeqIO->new( -format => 'fasta',
#                                   -file   => "$scaff_dir/$scaffold" );
#         my $scaffold = Bio::Seq->new( -display_id => "scaffold_${ref_id}" . ++$i );
#         while ( my $seq = $io->next_seq ) {
#             if ( $scaffold->length() > 0 ) {
#                 $scaffold->seq( $scaffold->seq() . 'N' x 100 . $seq->seq() );
#             }
#             else {
#                 $scaffold->seq( $seq->seq() );
#             }
#         }
#         $multiFastaIO->write_seq($scaffold);
#     }

#     #  append singleton scaffolds to scaffolds.fasta
#     foreach my $scaffold (@singletons) {
#         my $io = Bio::SeqIO->new( -format => 'fasta',
#                                   -file   => "$scaff_dir/$scaffold" );
#         my $seq = $io->next_seq();
#         $seq->display_id( "scaffold_${ref_id}" . ++$i );
#         $multiFastaIO->write_seq($seq);

#     }

#     open SIS_OUTPUT, "$scaff_dir/sis.sis"         or croak "Error opening $scaff_dir/sis.sis: $!";
#     open CONTIG_IDS, ">$scaff_dir/sis.contig_ids" or croak "Error opening $scaff_dir/sis.contig_ids: $!";
#     while ( my $line = <SIS_OUTPUT> ) {
#         next if ( $line =~ /^>|^$/ );
#         my $contig = ( split( / /, $line ) )[0];
#         print CONTIG_IDS $contig, "\n";
#     }
#     close SIS_OUTPUT;
#     close CONTIG_IDS;

# }
# """
