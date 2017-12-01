#!/usr/bun/env python3

######################################################################
#
# BugBuilder wrapper for abyss-sealer
#
# This file is part of BugBuilder (https://github.com/jamesabbott/BugBuilder)
# and is distributed under the Artistic License 2.0 - see accompanying LICENSE
# file for details
#
######################################################################

def run(config, args, results, reads_ns, scaffolds, finisher_dir, logger):
    logger.info("Using abyss-sealer to finish our scaffolds: %s", scaffolds)
                             'tmpdir=s'   => \$tmpdir,
                             'encoding=s' => \$encoding,
                             'threads=s'  => \$threads,
                             'help'       => \$help,
                             'man'        => \$man,
                           );

    # Following contributed by Andrey Tovchigrechko <andreyto@gmail.com>
    filled_prefix = "gap_sealed"

    # First, try to use reads error-corrected by Spades if they exits
    reads = glob.glob(
        "{args.tmp_dir}/spades/corrected/*.cor.fastq.gz".format(locals()))
    if len(reads) == 0:
        # Otherwise, use the input reads and rely on error-correction built into abyss-sealer
        reads = [args.fastq1, args.fastq2]

    enc_opt = '--standard-quality'
    if reads_ns.encoding != 'sanger':
        enc_opt = '--illumina-quality'
    # --no-trim-masked is critical. Contrary to --help output, trimming of
    # lower case nucleotides is the default in version 1.9.0, and everything is lower case after
    # the contig orientation step, except for the Ns.
    cmd = "{config.abyss_sealer} --no-trim-masked -k95 -k90 -k80 -P 10 -j ${args.threads}threads $enc_opt -o $filled_prefix -S $tmpdir/scaffolds.fasta ".join( " ", @reads );
    for cmd in [fgap_cmd]:
        logger.debug(cmd)
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    with open(seal_scaffs, "r") as inf, open("100bpNs.fasta", "w") as outf:
        for rec in SeqIO.parse(inf, fasta):
            # expand gaps of a single N to 100 Ns as per convention for unknown gap sizes,
            # since these are probably the 1 base remaining following gap closure...
            # gap = 'N' x 100;
            # new_seq = re.sub("N*", rec.seq, gap )
            SeqIO.write(outf, rec, "fasta")
