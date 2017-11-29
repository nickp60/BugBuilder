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
from BugBuilder import run_pilon as rp
from argparse import Namespace
logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class test_run_pilon(unittest.TestCase):
    """
    """
    def setUp(self):
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references")
        self.pilon_dir = os.path.join(os.path.dirname(__file__), "sis_tests")
        # os.makedirs(self.test_dir, exist_ok=True)
        self.startTime = time.time() # for timing
        self.to_be_removed = []

    def test_make_pilon_cmd(self):
        """
        """
        self.assertEqual(
            "pilon --genome contigs.fasta --bam ./out/scaffolds.sorted.bam " +
            "--changes --vcf --tracks --threads 7 --outdir ./out/ > " +
            "./out/pilon.log 2>&1",
            rp.make_pilon_cmd(pilon_exe="pilon", finisher_dir="./out/",
                           scaffolds="contigs.fasta", threads=7)
        )


    def test_make_pilon_bwa_cmds(self):
        args = Namespace(fastq1="reads1.fastq", fastq2="reads2.fastq")
        ref_cmds = [
            "bwa index contigs.fasta > ./out/bwa_index.log 2>&1",
            "bwa mem -t 7 contigs.fasta reads1.fastq " +
                  "reads2.fastq 2> ./out/bwa_mem.log | " +
                  "samtools view -bS - > ./out/scaffolds.bam " +
                  "2>./out/samtools_view.log",
            "samtools sort ./out/scaffolds.bam -o ./out/scaffolds.sorted.bam" +
            " > ./out/samtools_sort.log 2>&1)",
            "samtools index ./out/scaffolds.sorted.bam > " +
            "./out/samtools_index.log 2>&1"
        ]
        cmds, outbam = rp.make_pilon_bwa_cmds(bwa_exe="bwa", samtools_exe="samtools",
                                      args=args, scaffolds="contigs.fasta",
                                      finisher_dir="./out/", threads=7)
        for i, cmd in enumerate(cmds):
            self.assertEqual(cmd, ref_cmds[i])
        self.assertEqual(outbam, "./out/scaffolds.sorted.bam")


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
