#!/usr/bin/env python3

######################################################################
#
# This template can be used to wrap a new genome finisher to BugBuilder
#
# This file is part of BugBuilder (https://github.com/jamesabbott/BugBuilder)
# and is distributed under the Artistic License 2.0 - see accompanying LICENSE
# file for details
#
######################################################################
import subprocess
import os
import sys

def make_preprocessing_cmds():
    return None


def make_finisher_cmd():
    return "command used to run finisher"


def run(config, args, results, reads_ns, scaffolds, finisher_dir, logger):
    """keep these args the same"""
    logger.info("Using [FINISHER] to finish our scaffolds: %s", scaffolds)
    cmds, sorted_bam = make_preprocessing_cmds()
    cmds.append(
        make_finisher_cmd())
    for cmd in cmds:
        logger.debug(cmd)
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
