#!/usr/bin/env python3

######################################################################
#
# BugBuilder wrapper for FGAP
#
# This file is part of BugBuilder (https://github.com/jamesabbott/BugBuilder)
# and is distributed under the Artistic License 2.0 - see accompanying LICENSE
# file for details
#
######################################################################
import subprocess
import os
import sys

def make_fgap_cmd(fgap_exe, finisher_dir, scaffolds, threads):
    """ with conda, there is a lovely runner script so we dont have
    to tango with java
    """
    return str("{fgap_exe} --genome {scaffolds} " +
               "--bam {finisher_dir}scaffolds.sorted.bam --changes --vcf " +
               "--tracks --threads {threads} --outdir {finisher_dir} > " +
               "{finisher_dir}pilon.log 2>&1").format(**locals())


def run(config, args, results, reads_ns, scaffolds, finisher_dir, logger):
    logger.info("Using fgap to finish our scaffolds: %s", scaffolds)
    fgap_cmd = make_fgap_cmd(fgap_exe=config.pilon,
                         finisher_dir=finisher_dir,
                         scaffolds=scaffolds,
                         threads=args.threads))
    for cmd in [fgap_cmd]:
        logger.debug(cmd)
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
