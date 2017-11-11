# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas
"""
import time
import sys
import logging
import os
import unittest
import shutil
import BugBuilder.BugBuilder as bb
import statistics

from .context import BugBuilder
from BugBuilder import run_sis as rs
from argparse import Namespace
logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class test_run_sis(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references")
        self.test_dir = os.path.join(os.path.dirname(__file__), "sis_tests")
        self.empty_config = os.path.join(self.ref_dir, "empty_config.yaml")
        self.filled_config = os.path.join(self.ref_dir, "semicomplete_config.yaml")
        self.ref_fasta = os.path.join(self.ref_dir, "AP017923.1.fasta")
        self.renaming_fq = os.path.join(self.ref_dir, "needs_renaming.fq")
        self.renamed = os.path.join(self.ref_dir, "renamed_ref.fq")
        self.fastq1 = os.path.join(self.ref_dir, "AP017923.1_reads1.fq")
        self.fastq2 = os.path.join(self.ref_dir, "AP017923.1_reads2.fq")
        self.args = Namespace()
        os.makedirs(self.test_dir, exist_ok=True)
        self.startTime = time.time() # for timing
        self.to_be_removed = []

    def test_make_sis_etc_cmds(self):
        """ test pandas import
        """
        ref_cmds = [
            "nucmer reffy contig.fa -p sisdir/sis 2>&1 > sisdir/nucmer.log",
            "delta-filter -1 sisdir/sis.delta 2> sisdir/delta-filter.log > sisdir/sis.filter",
            "show-coords sisdir/sis.filter 2> sisdir/show-coords.log > sisdir/sis.coords",
            "sis.py sisdir/sis.coords > sisdir/sis.sis",
             "multifasta.py sisdir/sis.sis contig.fa"
            ]

        test_args = Namespace(reference="reffy")
        config = bb.parse_config(self.filled_config)
        config.nucmer = "nucmer"
        config.show_coords = "show-coords"
        config.delta_filter = "delta-filter"
        config.sis = "sis.py"
        config.multifasta = "multifasta.py"
        cmds = rs.make_sis_etc_cmds(config=config, args=test_args, contigs="contig.fa", scaff_dir="sisdir")
        for i, cmd in enumerate(cmds):
            self.assertEqual(cmd, ref_cmds[i])


    def tearDown(self):
        """ delete temp files if no errors, and report elapsed time
        """
        for filename in self.to_be_removed:
            try:
                os.unlink(filename)
            except Exception as e: # could be IsADirectpry or PErmissionsError
                print (e)
                shutil.rmtree(filename)
        t = time.time() - self.startTime
        print("%s: %.3f" % (self.id(), t))


if __name__ == '__main__':
    unittest.main()
