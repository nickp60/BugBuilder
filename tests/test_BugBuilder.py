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
        self.encoding_dir = os.path.join(self.ref_dir, "encoding")
        self.test_dir = os.path.join(os.path.dirname(__file__), "tmp_tests")
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
        tmp_dir = os.path.join(self.test_dir, "tmp_setup")
        test_args = Namespace(fastq1=self.fastq1, fastq2=self.fastq2,
                              reference=self.ref_fasta, mode="draft",
                              tmp_dir=tmp_dir, long_fastq=None,
                              de_fere_contigs=None)
        bb.setup_tmp_dir(args=test_args, output_root=tmp_dir, logger=logger)
        intended_files = [
            os.path.join(self.test_dir, "tmp_setup",
                         os.path.basename(self.fastq1)),
            os.path.join(self.test_dir, "tmp_setup",
                         os.path.basename(self.fastq2)),
            os.path.join(self.test_dir, "tmp_setup",
                         os.path.basename(self.ref_fasta)),
        ]
        for idx, f in enumerate([self.fastq1, self.fastq2, self.ref_fasta]):
            self.assertEqual(bb.md5(intended_files[idx]), bb.md5(f))
        self.to_be_removed.append(tmp_dir)

    def test_setup_tmp_dir_existingdir(self):
        with self.assertRaises(SystemExit):
            bad_test_args = Namespace(fastq1=self.fastq1, fastq2=self.fastq2,
                                      reference=self.ref_fasta, mode="draft",
                                      tmp_dir=self.test_dir, long_fastq=None,
                                      de_fere_contigs=None)
            bb.setup_tmp_dir(args=bad_test_args,
                             output_root=bad_test_args.tmp_dir, logger=logger)


    def test_fastq_needs_newname(self):
        test_args = Namespace(fastq1=self.renaming_fq, fastq2=self.fastq2,
                              reference=None, mode="draft")
        self.assertTrue(bb.fastq_needs_newname(test_args))

    def test_rename_fastq_seqids(self):
        test_args = Namespace(fastq1=self.renaming_fq, fastq2=None)
        bb.rename_fastq_seqids(args=test_args)
        self.assertEqual(
            bb.md5(self.renamed), bb.md5(test_args.fastq1))
        self.assertEqual(
            test_args.fastq1,
            os.path.splitext(self.renaming_fq)[0] +"_renamed.fq")
        self.to_be_removed.append(os.path.splitext(self.renaming_fq)[0] +
                                  "_renamed.fq")

    def test_get_fastq_read_len(self):
        lens = bb.get_read_lens_from_fastq(fastq1=self.fastq1, logger=logger)
        self.assertEqual(150, statistics.mean(lens))


    def test_id_fastq_encoding(self):
        tests = {
            "illumina": "illumina13.fq",
            "illumina": "illumina15.fq",
            "sanger": "illumina18_sanger.fq",
            "sanger": "sanger.fq",
            "solexa": "solexa.fq"
        }
        for k, v in tests.items():
            test_args = Namespace(
                fastq1=os.path.join(self.encoding_dir, v), long_fastq=None,
                de_fere_contigs=None)
            self.assertEqual(k, bb.id_fastq_encoding(args=test_args, logger=logger))


    def test_assess_reads(self):
        test_args = Namespace(
            fastq1=self.fastq1, fastq2=self.fastq2, long_fastq=None,
            de_fere_contigs=None, reference=self.ref_fasta, genome_size=0)
        config = bb.parse_config(self.filled_config)
        reads_ns = bb.assess_reads(args=test_args, config=config,
                                platform="illumina", logger=logger)
        self.assertEqual(9.996, reads_ns.coverage)
        self.assertEqual("sanger", reads_ns.encoding)
        self.assertEqual("long_illumina", reads_ns.lib_type)

    def test_check_ref_needed(self):
        test_args = Namespace(genome_size=0, reference=None)
        with self.assertRaises(ValueError):
            bb.check_ref_needed(args=test_args, lib_type="long")

    def test_get_config_path(self):
        pass

    def test_check_assemblers_bad_assembler(self):
        test_args = Namespace(
            fastq1=self.fastq1, fastq2=self.fastq2, long_fastq=None,
            de_fere_contigs=None, reference=self.ref_fasta, genome_size=0)
        config = bb.parse_config(self.filled_config)
        reads_ns = bb.assess_reads(args=test_args, config=config,
                                platform="illumina", logger=logger)
        test_args.assemblers = ["notAnAssembler"]
        with self.assertRaises(ValueError):
            bb.check_assemblers(args=test_args, config=config,
                                reads_ns=reads_ns, logger=logger)

    def test_check_assemblers_none_provided(self):
        test_args = Namespace(
            fastq1=self.fastq1, fastq2=self.fastq2, long_fastq=None,
            de_fere_contigs=None, reference=self.ref_fasta, genome_size=0)
        config = bb.parse_config(self.filled_config)
        reads_ns = bb.assess_reads(args=test_args, config=config,
                                   platform="illumina", logger=logger)
        test_args.assemblers = []
        reads_ns.lib_type = "hybrid"
        bb.check_assemblers(args=test_args, config=config,
                            reads_ns=reads_ns, logger=logger)
        self.assertEqual(test_args.assemblers, ["masurca", "spades"])

    def test_check_assemblers_too_long_reads(self):
        test_args = Namespace(
            fastq1=self.fastq1, fastq2=self.fastq2, long_fastq=None,
            de_fere_contigs=None, reference=self.ref_fasta, genome_size=0)
        config = bb.parse_config(self.filled_config)
        reads_ns = bb.assess_reads(args=test_args, config=config,
                                   platform="illumina", logger=logger)
        test_args.assemblers = ["spades"]
        reads_ns.mean_read_length = 500
        with self.assertRaises(ValueError):
            bb.check_assemblers(args=test_args, config=config,
                                reads_ns=reads_ns, logger=logger)

    def test_check_assemblers_need_insertsize(self):
        pass

    def test_check_assemblers_too_short(self):
        pass

    def test_get_scaffolder_and_linkage(self):
        test_args = Namespace(
            scaffolder="sis", reference=self.ref_fasta)
        config = bb.parse_config(self.filled_config)
        bb.get_scaffolder_and_linkage(args=test_args, config=config,
                                      paired=True, logger=logger)


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
