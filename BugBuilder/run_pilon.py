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

def make_pilon_bwa_cmds(bwa_exe, samtools_exe, args, scaffolds, finisher_dir, threads):
    index_cmd = "{0} index {1} > {2}bwa_index.log 2>&1".format(
        bwa_exe, scaffolds, finisher_dir)
    mem_cmd = "{exe} mem -t {threads} {scaffs} {read1} {read2} 2> {outdir}bwa_mem.log | {sam_exe} view -bS - > {outdir}scaffolds.bam 2>{outdir}samtools_view.log".format(
        exe=bwa_exe, threads=args.threads,scaffs=scaffolds,
        read1=args.fastq1, read2=args.fastq2, outdir=finisher_dir,
        sam_exe=samtools_exe)
    sort_cmd = "{sam_exe} sort {outdir}scaffolds.bam -o {outdir}scaffolds.sorted.bam > {outdir}samtools_sort.log 2>&1".format(
        sam_exe=samtools_exe,  outdir=finisher_dir)
    idx_bam = "{sam_exe} index {outdir}scaffolds.sorted.bam >{outdir}samtools_index.log 2>&1".format(
        sam_exe=samtools_exe,  outdir=finisher_dir)
    cmd_list = [index_cmd, mem_cmd, sort_cmd, idx_bam]
    return (cmd_list, "{outdir}scaffolds.sorted.bam".format(outdir=finisher_dir))


def make_pilon_cmd(pilon_exe, finisher_dir, scaffolds, threads):
    """ with conda, there is a lovely runner script so we dont have
    to tango with java
    """
    return str("{pilon} --genome {scaffolds} " +
               "--bam {outdir}scaffolds.sorted.bam --changes --vcf " +
               "--tracks --threads {t} --outdir {outdir} > " +
               "{outdir}pilon.log 2>&1").format(
                   pilon=pilon_exe, outdir=finisher_dir,
                   scaffolds=scaffolds, t=threads)


def run(config, args, results, reads_ns, scaffolds, finisher_dir, logger):
    logger.info("Using Pilon to finish our scaffolds: %s", scaffolds)
    if not reads_ns.paired:
        raise ValueError("Error: paired reads do not seem to be available....\n")
    cmds, sorted_bam = make_pilon_bwa_cmds(
        bwa_exe=config.bwa, samtools_exe=config.samtools, args=args, scaffolds=scaffolds,
        finisher_dir=finisher_dir, threads=args.threads)
    cmds.append(
        make_pilon_cmd(pilon_exe=config.pilon,
                       finisher_dir=finisher_dir,
                       scaffolds=scaffolds,
                       threads=args.threads))
    for cmd in cmds:
        logger.debug(cmd)
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
