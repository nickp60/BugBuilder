#!/usr/bin/env python3

######################################################################
#
# BugBuilder wrapper for pilon
#
# This file is part of BugBuilder (https://github.com/jamesabbott/BugBuilder)
# and is distributed under the Artistic License 2.0 - see accompanying LICENSE
# file for details
#
######################################################################
import subprocess
import os
import sys

# def make_pilon_bwa_cmds(bwa_exe, samtools_exe, args, scaffolds, varcall_dir, threads):
#     index_cmd = str("{bwa_exe} index {scaffolds} > " +
#                     "{varcall_dir}bwa_index.log 2>&1").format(**locals())
#     mem_cmd = str("{bwa_exe} mem -t {threads} {scaffolds} {args.fastq1} " +
#                   "{args.fastq2} 2> {varcall_dir}bwa_mem.log | " +
#                   "{samtools_exe} view -bS - > {varcall_dir}scaffolds.bam " +
#                   "2>{varcall_dir}samtools_view.log").format(**locals())
#     sort_cmd = str("{samtools_exe} sort {varcall_dir}scaffolds.bam -o " +
#                    "{varcall_dir}scaffolds.sorted.bam > " +
#                    "{varcall_dir}samtools_sort.log 2>&1").format(**locals())
#     idx_bam = str("{samtools_exe} index {varcall_dir}scaffolds.sorted.bam > " +
#                   "{varcall_dir}samtools_index.log 2>&1").format(**locals())
#     cmd_list = [index_cmd, mem_cmd, sort_cmd, idx_bam]
#     return (cmd_list, "{varcall_dir}scaffolds.sorted.bam".format(**locals()))


def make_pilon_varcall_cmds(pilon_exe, varcall_dir, threads, reference_bam, reference):
    """ with conda, there is a lovely runner script so we dont have
    to tango with java
    """
    return [str("{pilon_exe} --genome {reference} " +
               "--bam {reference_bam} --variant --changes --vcf " +
               " --threads {threads} --outdir {varcall_dir} > " +
               "{varcall_dir}pilon.log 2>&1").format(**locals())]
    # my $pilon_dir  = $config->{'pilon_dir'};
    # my $vcflib_dir = $config->{'vcflib_dir'};
    # my $java       = $config->{'java'};

    # my $cmd =
    #     "${java} -jar ${pilon_dir}/pilon.jar --genome ${tmpdir}/var_pilon/reference.fasta"
    #   . " --bam ${tmpdir}/var_pilon/reference.bam --variant --vcf --output var "
    #   . "--threads $threads >${tmpdir}/var_pilon/pilon.log 2>&1";
    # system($cmd) == 0 or die "Error running pilon";

    # $cmd = "cat ${tmpdir}/var_pilon/var.vcf|$vcflib_dir/vcffilter -g 'GT = 1/1'|$vcflib_dir/vcffixup - "
    #   . "|$vcflib_dir/vcffilter -f 'AC > 0'> ${tmpdir}/var_pilon/var.filtered.vcf";


def run(config, args, results, reads_ns, reference, reference_bam, varcall_dir, logger):
    logger.info("Using Pilon to call variants between  reads and %s",
                reference)
    if not reads_ns.paired:
        raise ValueError("Error: paired reads do not seem to be available....\n")

    # cmds, sorted_bam = make_pilon_bwa_cmds(
    #     bwa_exe=config.bwa, samtools_exe=config.samtools, args=args, scaffolds=scaffolds,
    #     varcall_dir=varcall_dir, threads=args.threads)
    cmds = make_pilon_varcall_cmds(pilon_exe=config.pilon,
                                   varcall_dir=varcall_dir,
                                   reference=reference,
                                   threads=args.threads,
                                   reference_bam=reference_bam)
    for cmd in cmds:
        logger.debug(cmd)
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    return os.path.join(varcall_dir, "pilon.vcf")
