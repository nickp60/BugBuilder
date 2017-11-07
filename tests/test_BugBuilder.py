# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas
"""
import sys
import logging
import os
import unittest
import shutil
import BugBuilder.BugBuilder as bb

from .context import BugBuilder
from argparse import Namespace
logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class test_BugBuilder(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references")
        self.test_dir = os.path.join(os.path.dirname(__file__), "tmp_tests")
        self.empty_config = os.path.join(self.ref_dir, "empty_config.yaml")
        self.ref_fasta = os.path.join(self.ref_dir, "AP017923.1.fasta")
        self.renaming_fq = os.path.join(self.ref_dir, "needs_renaming.fq")
        self.fastq1 = os.path.join(self.ref_dir, "AP017923.1_reads1.fq")
        self.fastq2 = os.path.join(self.ref_dir, "AP017923.1_reads2.fq")
        self.args = Namespace()
        os.makedirs(self.test_dir, exist_ok=True)

    def test_parse_config(self):
        """ test pandas import
        """
        newpath_config = os.path.join(
            self.test_dir, "empty.yaml")
        shutil.copyfile(self.empty_config, newpath_config)
        config = bb.parse_config(newpath_config)
        self.assertEqual(config.STATUS, "INCOMPLETE")

    def test_fill_in_config(self):
        """ test pandas import
        """
        newpath_config = os.path.join(
            self.test_dir, "to_be_filled.yaml")
        shutil.copyfile(self.empty_config, newpath_config)
        config = bb.return_config(newpath_config, logger=logger)
        self.assertEqual(config.STATUS, "COMPLETE")

    def test_make_fastqc_cmd(self):
        test_args = Namespace(fastq1="reads1.fastq", fastq2="reads2.fastq",
                              long_fastq=None, threads=7)
        ref_cmd = "fastqc -t 7 --extract -o ./outdir/ reads1.fastq reads2.fastq > " + \
                  "./outdir/fastqc.log 2>&1"
        cmd = bb.make_fastqc_cmd(args=test_args, outdir="./outdir/")
        self.assertEqual(ref_cmd, cmd)

    # def test_n50(self):
    #     lengths  = [2, 2, 2, 3, 3, 4, 8, 8]
    #     self.assertEqual(bb.get_L50_N50(lengths), (2, 6))

    def test_match_assembler_args_unequal(self):
        test_args = Namespace(assembler=["spades"],
                              assembler_args=["too", "many", "args"])
        with self.assertRaises(ValueError):
            bb.match_assembler_args(test_args)

    def test_match_assembler_args_no_assembler(self):
        test_args = Namespace(assembler=[],
                              assembler_args=[])
        self.assertEqual(bb.match_assembler_args(test_args),
                         [None, None])

    def test_match_assembler_args_expected(self):
        test_args = Namespace(
            assembler=["spades", "sellotape"],
            assembler_args=["--careful -k 21,33", "more sellotape" ])
        self.assertEqual(bb.match_assembler_args(test_args),
                         [["spades", "--careful -k 21,33"],
                          ["sellotape", "more sellotape"]])

    def test_check_files_present(self):
        test_args = Namespace(fastq1=self.fastq1, fastq2=self.fastq2,
                              reference=self.ref_fasta, mode="draft")
        bb.check_files_present(test_args)

    def test_check_files_present_not_exist(self):
        test_args = Namespace(fastq1=self.fastq1, fastq2=self.fastq2,
                              reference="notarealfile", mode="draft")
        with self.assertRaises(ValueError):
            bb.check_files_present(test_args)

    def test_check_files_present_same_files(self):
        test_args = Namespace(fastq1=self.fastq1, fastq2=self.fastq1,
                              reference=self.ref_fasta, mode="draft")
        with self.assertRaises(ValueError):
            bb.check_files_present(test_args)

    def test_check_files_present_bad_mode(self):
        test_args = Namespace(fastq1=self.fastq1, fastq2=self.fastq1,
                              reference=None, mode="draft")
        with self.assertRaises(ValueError):
            bb.check_files_present(test_args)

    def test_setup_tmp_dir(self):
        pass

    def test_fastq_needs_newname(self):
        test_args = Namespace(fastq1=self.renaming_fq, fastq2=self.fastq1,
                              reference=None, mode="draft")
        self.assertTrue(bb.fastq_needs_newname(test_args))

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
