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
from .context import BugBuilder
import BugBuilder.BugBuilder as bb
import BugBuilder.shared_methods as sm
import statistics

from argparse import Namespace
logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class test_shared_methods(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references")
        self.test_dir = os.path.join(os.path.dirname(__file__), "tmp_tests")
        self.startTime = time.time() # for timing
        self.ref_fasta = os.path.join(self.ref_dir, "AP017923.1.fasta")
        self.coords = os.path.join(self.ref_dir, "origin", "ori.coords")
        self.contigs = os.path.join(self.ref_dir, "contigs.fasta")
        self.to_be_removed = []

    def test_make_nucmer_delta_show_cmds(self):
        """ test pandas import
        """
        config = Namespace(nucmer = "nuc", delta_filter="df", show_coords="sc")
        cmds, coords_file = sm.make_nucmer_delta_show_cmds(
            config, ref="genome", query="contigs",
            out_dir="res", prefix="sis", header=False)

        ref_cmds = [
            "nuc genome contigs -p res/sis > res/nucmer.log 2>&1",
            "df -1 res/sis.delta 2> res/delta-filter.log > res/sis.filter",
            "sc -H res/sis.filter 2> res/show-coords.log > res/sis.coords"
        ]
        for idx, cmd in enumerate(cmds):
            self.assertEqual(ref_cmds[idx], cmd)
        self.assertEqual(coords_file, os.path.join("res", "sis.coords"))

    def test_run_nucmer_cmds(self):
        """
        """
        if shutil.which("nucmer") is None:
            cmds = ["echo 'I wish nucmer was here'",
                    "echo 'nucmer always loved a good unit test'"]
        else:
            config = Namespace(nucmer = shutil.which("nucmer"),
                               delta_filter=shutil.which("delta-filter"),
                               show_coords=shutil.which("show-coords"))
            cmds, coords = sm.make_nucmer_delta_show_cmds(
                config, ref=self.ref_fasta, query=self.contigs,
                out_dir=self.test_dir, prefix="test_nuc", header=False)
            self.assertEqual(coords,
                             os.path.join(self.test_dir, "test_nuc.coords"))
        sm.run_nucmer_cmds(cmds=cmds, logger=logger)

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
