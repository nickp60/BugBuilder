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
    index_cmd = str("{bwa_exe} index {scaffolds} > " +
                    "{finisher_dir}bwa_index.log 2>&1").format(**locals())
    mem_cmd = str("{bwa_exe} mem -t {threads} {scaffolds} {args.fastq1} " +
                  "{args.fastq2} 2> {finisher_dir}bwa_mem.log | " +
                  "{samtools_exe} view -bS - > {finisher_dir}scaffolds.bam " +
                  "2>{finisher_dir}samtools_view.log").format(**locals())
    sort_cmd = str("{samtools_exe} sort {finisher_dir}scaffolds.bam -o " +
                   "{finisher_dir}scaffolds.sorted.bam > " +
                   "{finisher_dir}samtools_sort.log 2>&1)").format(**locals())
    idx_bam = str("{samtools_exe} index {finisher_dir}scaffolds.sorted.bam > " +
                  "{finisher_dir}samtools_index.log 2>&1").format(**locals())
    cmd_list = [index_cmd, mem_cmd, sort_cmd, idx_bam]
    return (cmd_list, "{finisher_dir}scaffolds.sorted.bam".format(**locals()))


def make_pilon_cmd(pilon_exe, finisher_dir, scaffolds, threads):
    """ with conda, there is a lovely runner script so we dont have
    to tango with java
    """
    return str("{pilon_exe} --genome {scaffolds} " +
               "--bam {finisher_dir}scaffolds.sorted.bam --changes --vcf " +
               "--tracks --threads {threads} --outdir {finisher_dir} > " +
               "{finisher_dir}pilon.log 2>&1").format(**locals())


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
