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
from BugBuilder import run_fgap as rf
from argparse import Namespace
logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class test_run_fgap(unittest.TestCase):
    """
    """
    def setUp(self):
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references")
        self.test_dir = os.path.join(os.path.dirname(__file__), "fgap_tests")
        self.startTime = time.time() # for timing
        self.to_be_removed = []


    def test_make_fgap_cmd(self):
        self.assertEqual(
            rf.make_fgap_cmd(fgap_exe="FGAP",
                             finisher_dir="./fgap/",
                             refs_string="ref1.fasta,ref2.fasta",
                             scaffolds="scaffs.fasta", threads=7),
            "FGAP '--draft-file scaffs.fasta " +
            "--datasets-files ref1.fasta,ref2.fasta --threads 7 " +
            "--output-prefix ./fgap/fgap' > ./fgap/FGAP.log 2>&1"
        )

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
