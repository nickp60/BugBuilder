# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas
"""
import time
import sys
import logging
import os
from unittest import mock
from yaml import YAMLError
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
        self.aa_dir = os.path.join(os.path.dirname(__file__),
                                   "already_assembled_dir")
        self.empty_config = os.path.join(self.ref_dir, "empty_config.yaml")
        self.broken_config = os.path.join(self.ref_dir, "broken_config.yaml")
        # the active config actually gets used to run sys commands so feel free
        # to run `mv tests/tmp_tests/to_be_filled.yaml ./tests/references/semicomplete_config.yaml`
        # if you have a error from a wrong executable
        self.active_config = os.path.join(self.ref_dir,
                                          "semicomplete_config.yaml")
        # static config should not be changed, it is just used to test parsing
        self.static_config = os.path.join(self.ref_dir, "static_config.yaml")
        self.ref_fasta = os.path.join(self.ref_dir, "2chrom.fasta")
        self.coords = os.path.join(self.ref_dir, "origin", "ori.coords")
        self.picard_stats = os.path.join(self.ref_dir, "picard_insert.txt")
        self.ref_split = os.path.join(self.ref_dir, "2chrom.fasta")
        self.ref_split_1 = os.path.join(self.ref_dir, "2chrom1.fq")
        self.ref_split_2 = os.path.join(self.ref_dir, "2chrom2.fq")
        self.sickle_log = os.path.join(self.ref_dir, "sickle.log")
        self.sickle_bad_log = os.path.join(self.ref_dir, "sickle_bad.log")
        self.contigs = os.path.join(self.ref_dir, "contigs.fasta")
        self.scaffold = os.path.join(self.ref_dir, "scaffs.fasta")
        self.contigs_split_ori = os.path.join(self.ref_dir,
                                              "contigs_split_ori.fasta")
        self.mapped_bam = os.path.join(self.ref_dir, "mapped.bam")
        self.contigs_to_scaf = os.path.join(self.ref_dir,
                                            "contigs_to_scaffold.fasta")
        self.distant_contigs = os.path.join(self.ref_dir,
                                            "distant_contigs.fasta")
        self.renaming_fq = os.path.join(self.ref_dir, "needs_renaming.fq")
        self.renamed = os.path.join(self.ref_dir, "renamed_ref.fq")
        self.fastq1 = os.path.join(self.ref_dir, "AP017923.1_reads1.fq")
        self.fastq2 = os.path.join(self.ref_dir, "AP017923.1_reads2.fq")
        self.args = Namespace()
        os.makedirs(self.test_dir, exist_ok=True)
        self.startTime = time.time() # for timing
        self.to_be_removed = []
        # #  I wish there was a way around this
        # sis_path = os.path.join(self.test_dir, "sis_1")
        # if os.path.exists(sis_path):
        #     shutil.rmtree(sis_path)

    def test_parse_config(self):
        """
        """
        config = bb.parse_config(self.empty_config)
        self.assertEqual(config.STATUS, "INCOMPLETE")

    def test_parse_config_broken(self):
        """
        """
        with self.assertRaises(YAMLError):
            config = bb.parse_config(self.broken_config)


    def test_fill_in_config(self):
        """ test pandas import
        """
        newpath_config = os.path.join(
            self.test_dir, "to_be_filled.yaml")
        shutil.copyfile(self.empty_config, newpath_config)
        config = bb.return_config(newpath_config,
                                  hardfail=False, logger=logger)
        self.assertEqual(config.STATUS, "COMPLETE")
        self.to_be_removed.append(newpath_config)

    def test_make_fastqc_cmd(self):
        test_args = Namespace(fastq1="reads1.fastq", fastq2="reads2.fastq",
                              long_fastq=None, threads=7)
        ref_cmd = "fastqc -t 7 --extract -o ./outdir/ reads1.fastq " + \
                  "reads2.fastq > ./outdir/fastqc.log 2>&1"
        cmd = bb.make_fastqc_cmd(exe="fastqc", args=test_args,
                                 outdir="./outdir/")
        self.assertEqual(ref_cmd, cmd)

    # def test_n50(self):
    #     lengths  = [2, 2, 2, 3, 3, 4, 8, 8]
    #     self.assertEqual(bb.get_L50_N50(lengths), (2, 6))

    def test_match_assembler_args_unequal(self):
        test_args = Namespace(assemblers=["spades"],
                              assembler_args=["too", "many", "args"])
        with self.assertRaises(ValueError):
            bb.match_assembler_args(test_args)

    def test_match_assembler_args_no_assembler(self):
        test_args = Namespace(assemblers=[],
                              assembler_args=[])
        self.assertEqual(bb.match_assembler_args(test_args),
                         [None, None])

    def test_match_assembler_args_expected(self):
        test_args = Namespace(
            assemblers=["spades", "sellotape"],
            assembler_args=["--careful -k 21,33", "more sellotape" ])
        self.assertEqual(bb.match_assembler_args(test_args),
                         [["spades", "--careful -k 21,33"],
                          ["sellotape", "more sellotape"]])

    def test_check_files_present(self):
        test_args = Namespace(fastq1=self.fastq1, fastq2=self.fastq2,
                              references=[self.ref_fasta], mode="draft")
        bb.check_files_present(test_args)

    def test_check_files_present_not_exist(self):
        test_args = Namespace(fastq1=self.fastq1, fastq2=self.fastq2,
                              references=["notarealfile"], mode="draft")
        with self.assertRaises(ValueError):
            bb.check_files_present(test_args)

    def test_check_files_present_same_files(self):
        test_args = Namespace(fastq1=self.fastq1, fastq2=self.fastq1,
                              references=[self.ref_fasta], mode="draft")
        with self.assertRaises(ValueError):
            bb.check_files_present(test_args)

    def test_check_files_present_bad_mode(self):
        test_args = Namespace(fastq1=self.fastq1, fastq2=self.fastq1,
                              references=[], mode="draft")
        with self.assertRaises(ValueError):
            bb.check_files_present(test_args)

    def test_setup_tmp_dir(self):
        tmp_dir = os.path.join(self.test_dir, "tmp_setup")
        test_args = Namespace(fastq1=self.fastq1, fastq2=self.fastq2,
                              references=[self.ref_fasta], mode="draft",
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
            bad_test_args = Namespace(
                fastq1=self.fastq1, fastq2=self.fastq2,
                references=[self.ref_fasta], mode="draft",
                tmp_dir=self.test_dir, long_fastq=None,
                de_fere_contigs=None)
            bb.setup_tmp_dir(args=bad_test_args,
                             output_root=bad_test_args.tmp_dir, logger=logger)


    def test_fastq_needs_newname(self):
        test_args = Namespace(fastq1=self.renaming_fq, fastq2=self.fastq2,
                              references=[], mode="draft")
        self.assertTrue(bb.fastq_needs_newname(test_args))

    def test_fastq_dont_needs_newname(self):
        test_args = Namespace(fastq1=self.fastq1, fastq2=self.fastq2,
                              references=[], mode="draft")
        self.assertFalse(bb.fastq_needs_newname(test_args))

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
            self.assertEqual(
                k,
                bb.id_fastq_encoding(args=test_args, logger=logger))

    def test_assess_reads(self):
        test_args = Namespace(
            fastq1=self.fastq1, fastq2=self.fastq2, long_fastq=None,
            de_fere_contigs=None, references=[self.ref_fasta], genome_size=0)
        config = bb.parse_config(self.active_config)
        reads_ns = bb.assess_reads(args=test_args, config=config,
                                platform="illumina", logger=logger)
        self.assertEqual(4.998, reads_ns.coverage)
        self.assertEqual("sanger", reads_ns.encoding)
        self.assertEqual("long_illumina", reads_ns.lib_type)

    # def test_assess_reads_nanopore(self):
    #     test_args = Namespace(
    #         fastq1=self.fastq1, fastq2=self.fastq2, long_fastq=None,
    #         de_fere_contigs=None, references=[self.ref_fasta], genome_size=0)
    #     config = bb.parse_config(self.active_config)
    #     reads_ns = bb.assess_reads(args=test_args, config=config,
    #                             platform="illumina", logger=logger)
    #     self.assertEqual(9.996, reads_ns.coverage)
    #     self.assertEqual("sanger", reads_ns.encoding)
    #     self.assertEqual("long_illumina", reads_ns.lib_type)


    def test_check_ref_needed(self):
        test_args = Namespace(genome_size=0, references=[])
        with self.assertRaises(ValueError):
            bb.check_ref_required(args=test_args, lib_type="long")

    def test_get_config_path(self):
        pass

    def test_check_assemblers_too_many(self):
        test_args = Namespace(
            fastq1=self.fastq1, fastq2=self.fastq2, long_fastq=None,
            de_fere_contigs=None, references=[self.ref_fasta], genome_size=0)
        config = bb.parse_config(self.active_config)
        reads_ns = bb.assess_reads(args=test_args, config=config,
                                platform="illumina", logger=logger)
        test_args.assemblers = ["abyss", "spades", "too_many"]
        with self.assertRaises(ValueError):
            bb.check_and_get_assemblers(args=test_args, config=config,
                                reads_ns=reads_ns, logger=logger)

    # def test_check_assemblers_bad_assembler(self):
    #     test_args = Namespace(
    #         fastq1=self.fastq1, fastq2=self.fastq2, long_fastq=None,
    #         de_fere_contigs=None, references=[self.ref_fasta], genome_size=0)
    #     config = bb.parse_config(self.active_config)
    #     reads_ns = bb.assess_reads(args=test_args, config=config,
    #                             platform="illumina", logger=logger)
    #     test_args.assemblers = ["notAnAssembler"]
    #     with self.assertRaises(ValueError):
    #         bb.check_and_get_assemblers(args=test_args, config=config,
    #                             reads_ns=reads_ns, logger=logger)

    def test_check_assemblers_none_provided(self):
        test_args = Namespace(
            fastq1=self.fastq1, fastq2=self.fastq2, long_fastq=None,
            de_fere_contigs=None, references=[self.ref_fasta], genome_size=0)
        config = bb.parse_config(self.active_config)
        reads_ns = bb.assess_reads(args=test_args, config=config,
                                   platform="illumina", logger=logger)
        test_args.assemblers = []
        reads_ns.lib_type = "hybrid"
        bb.check_and_get_assemblers(args=test_args, config=config,
                            reads_ns=reads_ns, logger=logger)
        self.assertEqual(test_args.assemblers, ["masurca", "spades"])

    def test_check_assemblers_too_long_reads(self):
        test_args = Namespace(
            fastq1=self.fastq1, fastq2=self.fastq2, long_fastq=None,
            de_fere_contigs=None, references=[self.ref_fasta], genome_size=0)
        config = bb.parse_config(self.active_config)
        reads_ns = bb.assess_reads(args=test_args, config=config,
                                   platform="illumina", logger=logger)
        test_args.assemblers = ["spades"]
        reads_ns.read_length_mean = 500
        with self.assertRaises(ValueError):
            bb.check_and_get_assemblers(args=test_args, config=config,
                                reads_ns=reads_ns, logger=logger)

    def test_check_assemblers_need_insertsize(self):
        test_args = Namespace(assemblers=['mascura'],
            fastq1=self.fastq1, fastq2=self.fastq2, long_fastq=None,
            de_fere_contigs=None, references=[], genome_size=0)
        config = bb.parse_config(self.active_config)
        reads_ns = bb.assess_reads(args=test_args, config=config,
                                   platform="illumina", logger=logger)
        test_args.assemblers = ["spades"]
        reads_ns.read_length_mean = 500
        with self.assertRaises(ValueError):
            bb.check_and_get_assemblers(args=test_args, config=config,
                                reads_ns=reads_ns, logger=logger)

    def test_check_assemblers_too_short(self):
        test_args = Namespace(assemblers=['PBcR'],
            fastq1=self.fastq1, fastq2=self.fastq2, long_fastq=None,
            de_fere_contigs=None, references=[], genome_size=0)
        config = bb.parse_config(self.active_config)
        reads_ns = bb.assess_reads(args=test_args, config=config,
                                   platform="illumina", logger=logger)
        test_args.assemblers = ["spades"]
        reads_ns.read_length_mean = 400
        with self.assertRaises(ValueError):
            bb.check_and_get_assemblers(args=test_args, config=config,
                                reads_ns=reads_ns, logger=logger)

    def test_get_scaffolder_and_linkage(self):
        test_args = Namespace(
            scaffolder="sis", references=[self.ref_fasta])
        config = bb.parse_config(self.active_config)
        bb.get_scaffolder_and_linkage(args=test_args, config=config,
                                      paired=True, logger=logger)

    def test_check_already_assembled_bad_path(self):
        test_args = Namespace(
            already_assembled_dirs=["notAnActualDir"], assemblers=["spades"])
        config = bb.parse_config(self.active_config)
        with self.assertRaises(FileNotFoundError):
            bb.check_already_assembled_dirs(
                args=test_args, config=config, logger=logger)

    def test_check_already_assembled_missing_file(self):
        test_args = Namespace(
            already_assembled_dirs=[self.ref_dir], assemblers=["spades"])
        config = bb.parse_config(self.active_config)
        with self.assertRaises(FileNotFoundError):
            bb.check_already_assembled_dirs(
                args=test_args, config=config, logger=logger)

    def test_check_already_assembled_uneven_args(self):
        test_args = Namespace(
            already_assembled_dirs=[self.ref_dir],
            assemblers=["spades", "masurca"])
        config = bb.parse_config(self.active_config)
        with self.assertRaises(ValueError):
            bb.check_already_assembled_dirs(
                args=test_args, config=config, logger=logger)

    def test_check_already_assembled(self):
        already_dir = os.path.join(self.ref_dir,"already_assembled_test")
        already_ctg = os.path.join(already_dir,"contigs.fasta")
        already_scf = os.path.join(already_dir,"scaffolds.fasta")
        test_args = Namespace(
            already_assembled_dirs=[
                already_dir],
            assemblers=["spades" ])
        config = bb.parse_config(self.active_config)
        self.assertEqual(
            bb.check_already_assembled_dirs(
                args=test_args, config=config, logger=logger),
            [{"name": "spades",
             "contigs": already_ctg,
              "scaffolds": already_scf}])

    # def test_check_already_assembled_dirs(self):
    #     check_already_assembled_dirs(args, config, logger)

    def test_check_args_fail_assembler(self):
        config = bb.parse_config(self.active_config)
        test_args = Namespace(
            merger=None, assemblers=["spades", "shovels"])
        with self.assertRaises(ValueError):
            bb.check_args(test_args, config)

    def test_check_args_fail_finisher(self):
        config = bb.parse_config(self.active_config)
        test_args = Namespace(
            stages="f",
            scaffolder=None,
            downsample=0,
            merger=None,
            finisher=None,assemblers=[])
        with self.assertRaises(ValueError):
            bb.check_args(test_args, config)

    def test_parse_origin_from_coords(self):
        ref_dict = {'AP017923.1': 'NODE_4_length_5704_cov_4.881287:5704',
                    'AP017923.2': 'NODE_2_length_10885_cov_4.771743:869'}
        origin_dict = bb.parse_origin_from_coords(
            coords=self.coords, flex=5,
            reference=self.ref_split, logger=logger)
        self.assertEqual(origin_dict, ref_dict)

    def test_make_empty_results_object(self):
        key_list = [
            "organism",
            "assemblers_list",
            "assemblers_results_dict_list",
            "current_contigs",
            "current_contigs_source",
            "current_scaffolds",
            "current_scaffolds_source",
            "current_reference",
            "current_embl",
            "current_embl_source",
            "ID_OK",
            "reference_percent_sim",
            "old_contigs",
            "old_embl",
            "old_scaffolds",
            "old_references"]
        results = bb.make_empty_results_object()
        for k, v in sorted(vars(results).items()):
            if k not in key_list:
                raise KeyError("%s not a valid results key!" %k)

    def test_check_trim_length(self):
        self.assertEqual(
            25,
            bb.check_and_set_trim_length(args=Namespace(trim_length=50),
                              reads_ns=Namespace(read_length_mean=45),
                              logger=logger))

    def test_report_trim_quality_good(self):
        with self.assertLogs('', level='DEBUG') as cm:
            bb.report_trim_quality(trim_log=self.sickle_log, logger=logger)
            self.assertEqual(
                cm.output[-1],
                "DEBUG:root:Discarded 0% of the reads during trimming")

    def test_report_trim_quality_bad(self):
        with self.assertLogs('', level='DEBUG') as cm:
            bb.report_trim_quality(trim_log=self.sickle_bad_log, logger=logger)
            self.assertEqual(
                cm.output[-1],
                "WARNING:root:>10% of reads (26%) discarded during read " +
                "trimming. Please examine the FastQC outputs...")



    def test_replace_placeholders(self):
        test_args = Namespace(memory=1)
        string = "executable __MEMORY__"
        self.assertEqual(
            "executable 1",
            bb.replace_placeholders(string, args=test_args))

    def test_replace_placeholders_fail(self):
        string = "executable __MEMORY__ __CONTIGS__"
        test_args = Namespace(memory=1)
        with self.assertRaises(ValueError):
            bb.replace_placeholders(string, args=test_args)


    # def test_parse_available_noconfig(self):
    #     self.assertEqual(
    #         [],
    #         bb.parse_available("assemblers", None)
    #     )
    # @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
    #                  "Skipping this test on Travis CI. Too hard to debug")
    def test_parse_available_assemblers(self):
        self.assertEqual(
            ["abyss", "spades", "ribo"],
            bb.parse_available("assemblers", self.static_config)
        )

    def test_get_scaffolder_and_linkage(self):
        config = bb.parse_config(self.active_config)
        test_args = Namespace(scaffolder="sis",
                              references=[self.ref_fasta])
        scaffolder, linkage = bb.get_scaffolder_and_linkage(
            args=test_args, config=config, paired=True, logger=logger)
        self.assertEqual(linkage, "align_genus")
        self.assertEqual(scaffolder['name'], "SIS")

    def test_get_scaffolder_and_linkage_none(self):
        config = bb.parse_config(self.active_config)
        test_args = Namespace(scaffolder=None,
                              references=[self.ref_fasta])
        scaffolder, linkage = bb.get_scaffolder_and_linkage(
            args=test_args, config=config, paired=True, logger=logger)
        self.assertEqual(linkage, None)
        self.assertEqual(scaffolder, None)

    def test_get_merger_tool(self):
        config = bb.parse_config(self.active_config)
        test_args = Namespace(merger='gfinisher',
                              scaffolder=None,
                              references=[self.ref_fasta])
        merger = bb.get_merger_tool(
            args=test_args, config=config, paired=False)
        self.assertEqual(merger['name'], 'gfinisher')

    def test_get_finisher(self):
        config = bb.parse_config(self.active_config)
        test_args = Namespace(finisher='gapfiller',
                              scaffolder=None,
                              references=[self.ref_fasta])
        finisher = bb.get_finisher(args=test_args, config=config, paired=True)
        self.assertEqual(finisher['name'], 'gapfiller')

    def test_get_finisher_fail_needs_pair(self):
        config = bb.parse_config(self.active_config)
        test_args = Namespace(finisher='gapfiller',
                              scaffolder=None,
                              references=[self.ref_fasta])
        with self.assertRaises(ValueError):
            bb.get_finisher(args=test_args, config=config, paired=False)

    def test_get_finisher_fail_needs_refernce(self):
        config = bb.parse_config(self.active_config)
        test_args = Namespace(finisher='pilon',
                              scaffolder=None,
                              references=[])
        with self.assertRaises(ValueError):
            bb.get_finisher(args=test_args, config=config, paired=False)

    def test_get_varcaller_fail_needs_reference(self):
        config = bb.parse_config(self.active_config)
        test_args = Namespace(varcaller="pilon",
                              scaffolder=None,
                              references=[])
        with self.assertRaises(ValueError):
            bb.get_varcaller(args=test_args, config=config, paired=False)

    def test_get_varcaller(self):
        config = bb.parse_config(self.active_config)
        test_args = Namespace(varcaller="pilon",
                              scaffolder=None,
                              references=[self.ref_fasta])
        vc = bb.get_varcaller(args=test_args, config=config, paired=False)
        self.assertEqual("pilon", vc['name'])

    # @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
    #                  "Skipping this test on Travis CI. Too hard to debug")
    def test_select_tools(self):
        config = bb.parse_config(self.static_config)
        test_args = Namespace(scaffolder="sis",
                              assemblers=["spades"],
                              merger=None,
                              finisher="pilon",
                              varcaller="pilon",
                              references=[self.ref_fasta])
        tools = bb.select_tools(args=test_args, config=config,
                             reads_ns=Namespace(read_length_mean=45,
                                                paired=True),
                             logger=logger)
        self.assertEqual("SIS", tools.scaffolder['name'])

    def test_assembler_needs_downsampling(self):
        assemblers = [
            x for x in bb.parse_config(self.active_config).assemblers if \
            x['name' ]== "spades"]
        tools = Namespace(assemblers=assemblers)
        self.assertEqual(True, bb.assembler_needs_downsampling(tools))

    def run_fastqc(reads_ns, args, logger=None):
        pass

    def test_make_sickle_cmd(self):
        test_args = Namespace(
            fastq1=self.fastq1,
            trim_length=50,
            trim_qv=10,
            references=[self.ref_fasta])
        test_cmd = str("sickle se -f {0} -t sanger -q 10 -l 50 -o " +
                       "sickle/read1.fastq > sickle/sickle.log").format(
            self.fastq1)
        self.assertEqual(
            bb.make_sickle_cmd(
                args=test_args,
                reads_ns=Namespace(encoding="sanger", paired=False),
                out_dir="sickle/",
                paired=False),
            test_cmd)

    def quality_trim_reads(args, config, reads_ns, logger):
        pass

    def test_make_seqtk_ds_cmd(test):
        cmds = bb.make_seqtk_ds_cmd(
            args=Namespace(fastq1="reads_1.fq", fastq2="reads_2.fq"),
            reads_ns=Namespace(coverage=1000),
            new_coverage=100, outdir="./seqtk/",
            config=Namespace(seqtk="seqtk"), logger=logger)
        ref_cmds = [
            "seqtk sample -s 100 .1 reads_1.fq > ./seqtk/reads1.fq",
            "seqtk sample -s 100 .1 reads_2.fq > ./seqtk/reads2.fq"]
        for i, cmd in enumerate(cmds):
            self.assertEqual(cmd, ref_cmds[i])

    def downsample_reads(args, reads_ns, config, new_cov=100):
        pass

    def test_make_bwa_cmds_aln_se(self):
        config = bb.parse_config(self.active_config)
        config.bwa = "bWA"
        cmds, mapping_sam = bb.make_bwa_cmds(
            ref="ref.fasta",
            fastq1="reads_1.fq",
            fastq2=None,
            reads_ns=Namespace(read_length_mean=75),
            args=Namespace(threads=13),
            outdir="./bwa/",
            config=config)
        ref_cmds = [
            "bWA index ref.fasta  > ./bwa/bwa_index.log 2>&1",
            "bWA aln -t 13 ref.fasta reads_1.fq > ./bwa/read1.sai",
            "bWA samse ref.fasta ./bwa/read1.sai reads_1.fq " +
            "2> ./bwa/samse.log > ./bwa/mapping.sam"]
        for i, cmd in enumerate(cmds):
            self.assertEqual(cmd, ref_cmds[i])
        self.assertEqual(mapping_sam, "./bwa/mapping.sam")

    def test_make_bwa_cmds_aln_pe(self):
        config = bb.parse_config(self.active_config)
        config.bwa = "bWA"
        cmds, mapping_sam = bb.make_bwa_cmds(
            ref="ref.fasta",
            fastq1="reads_1.fq",
            fastq2="reads_2.fq",
            reads_ns=Namespace(read_length_mean=75),
            args=Namespace(threads=13),
            outdir="./bwa/",
            config=config)
        ref_cmds = [
            "bWA index ref.fasta  > ./bwa/bwa_index.log 2>&1",
            "bWA aln -t 13 ref.fasta reads_1.fq > ./bwa/read1.sai",
            "bWA aln -t 13 ref.fasta reads_2.fq > ./bwa/read2.sai",
            "bWA sampe ref.fasta ./bwa/read1.sai ./bwa/read2.sai reads_1.fq " +
            "reads_2.fq 2> ./bwa/sampe.log > ./bwa/mapping.sam"]
        for i, cmd in enumerate(cmds):
            self.assertEqual(cmd, ref_cmds[i])
        self.assertEqual(mapping_sam, "./bwa/mapping.sam")

    def test_make_bwa_cmds_mem_pe(self):
        config = bb.parse_config(self.active_config)
        config.bwa = "bWA"
        cmds, mapping_sam = bb.make_bwa_cmds(
            ref="ref.fasta",
            fastq1="reads_1.fq",
            fastq2="reads_2.fq",
            reads_ns=Namespace(read_length_mean=101),
            args=Namespace(threads=13),
            outdir="./bwa/",
            config=config)
        ref_cmds = [
            "bWA index ref.fasta  > ./bwa/bwa_index.log 2>&1",
            "bWA mem -t 13 -M ref.fasta reads_1.fq reads_2.fq > " +
            "./bwa/mapping.sam 2> ./bwa/bwa_mem.log"]
        for i, cmd in enumerate(cmds):
            self.assertEqual(cmd, ref_cmds[i])
        self.assertEqual(mapping_sam, "./bwa/mapping.sam")

    def test_make_bwa_cmds_mem_se(self):
        config = bb.parse_config(self.active_config)
        config.bwa = "bWA"
        cmds, mapping_sam = bb.make_bwa_cmds(
            ref="ref.fasta",
            fastq1="reads_1.fq",
            fastq2=None,
            reads_ns=Namespace(read_length_mean=101),
            args=Namespace(threads=13),
            outdir="./bwa/",
            config=config)
        ref_cmds = [
            "bWA index ref.fasta  > ./bwa/bwa_index.log 2>&1",
            "bWA mem -t 13 -M ref.fasta reads_1.fq > " +
            "./bwa/mapping.sam 2> ./bwa/bwa_mem.log"]
        for i, cmd in enumerate(cmds):
            self.assertEqual(cmd, ref_cmds[i])
        self.assertEqual(mapping_sam, "./bwa/mapping.sam")

    def test_make_samtools_cmds(self):
        self.assertEqual(
            ["sam view  -Shb in.sam | sam sort - > ./out/out.bam",
             "sam index ./out/out.bam 2> ./out/samtools_index.log"],
            bb.make_samtools_cmds(exe="sam", sam="in.sam", outdir="./out/",
                                out_bam="./out/out.bam")
        )

    def align_reads(dirname, reads_ns,  downsample, args, config, logger):
        pass

    def test_make_picard_stats_command(self):
        self.assertEqual(
            (str("picard CollectInsertSizeMetrics " +
              "INPUT=file.bam HISTOGRAM_FILE=./test/insert_histogram.pdf " +
              "OUTPUT=./test/insert_stats.txt QUIET=true VERBOSITY=ERROR " +
              "ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT > " +
              "./test/CollectInsertMetrics.log 2>&1 "),
             "./test/insert_stats.txt"),
            bb.make_picard_stats_command(bam="file.bam",
                                      config=Namespace(picard="picard"),
                                      picard_outdir="./test/")
            )


    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "Skipping this test on Travis CI. Too hard to debug")
    def test_get_insert_stats(self):
        test_args = Namespace(
            tmp_dir=self.test_dir)
        integration_config = bb.return_config(
            os.path.join(self.test_dir, "integration_config.yaml"), force=True,
            hardfail=False, logger=None)
        self.assertEqual(
            ('199.816601', '9.857524'),
            bb.get_insert_stats(bam=self.mapped_bam, config=integration_config,
                                args=test_args, logger=logger)
        )
        self.to_be_removed.append(
            os.path.join(self.test_dir, "insert_stats", ""))

    def test_parse_picard_insert_stats(self):
        self.assertEqual(('207.5',  '3.535534'),
                         bb.parse_picard_insert_stats(self.picard_stats))

    def replace_placeholders(string, config=None, reads_ns=None,
                             args=None, results=None):
        pass
    def get_assembler_cmds(assembler, assembler_args, args, config, reads_ns):
        pass
    def standardize_fasta_output(infile, outfile, ctype):
        pass
    def get_L50_N50(lengths):
        pass

    def test_update_results(self):
        results = bb.make_empty_results_object()
        self.assertTrue(results.current_scaffolds is None)
        path1 = "path/to/me/scaffs.fa"
        path2 = "path/to/me/other/scaffs.fa"
        source1 = "arbitrary point number 1"
        source2 = "look, more abitrary string!"
        bb.update_results(results=results,
                       thing="scaffolds",
                       path=path1,
                       source=source1)
        self.assertEqual(results.current_scaffolds, path1)
        self.assertEqual(results.current_scaffolds_source, source1)
        self.assertEqual(results.old_scaffolds, [(None, "init")])
        bb.update_results(results=results,
                       thing="scaffolds",
                       path=path2,
                       source=source2)
        self.assertEqual(results.current_scaffolds, path2)
        self.assertEqual(results.current_scaffolds_source, source2)
        self.assertEqual(results.old_scaffolds, [(None, "init"),
                                                 (path1, source1)])

    def test_get_contig_info(self):
        info_dict = bb.get_contig_info(self.contigs, x=200)
        for metric in [["count", 4],
                       ["all_lengths", [5704, 7626, 10885, 13713]]]:
            self.assertEqual(info_dict[metric[0]], metric[1])

    def test_get_contig_stats(self):
        sample_results = [
            ['', 'Old seqs', 'Old seqs > 200', 'New seqs', 'New seqs > 200'],
            ['count', 7, 7, 4, 4],
            ['Max Length', 13713, 13713, 13713, 13713],
            ['Assembly size', 38720, 38720, 37928, 37928],
            ['L50', 7626, 7626, 10885, 10885],
            ['N50', 6, 6, 3, 3],
            ["N's", 0, 0, 0, 0],
            ['path', 'references/contigs_to_scaffold.fasta', '', 'references/contigs.fasta', '']]
        results = bb.get_contig_stats(
            contigs=self.contigs, old_contigs=self.contigs_to_scaf)
        for i, line in enumerate(sample_results):
            self.assertEqual(results[i], line)


    def run_assembler(assembler, assembler_args, args, reads_ns, config, logger):
        pass
    def merge_assemblies(args, config, reads_ns, logger):
        pass

    @unittest.skipIf(shutil.which("blastn") is None,
                     "Cannot test check_id without blastn")
    def test_check_id_close(self):
        test_args = Namespace(tmp_dir=self.test_dir,
                              references=[self.ref_fasta],
                              )
        config = Namespace(blastn=shutil.which("blastn"))
        ID_OK, percent = bb.check_id(ref=self.ref_fasta,
                                     results=Namespace(current_scaffolds=self.contigs),
                                     args=test_args, contigs=self.contigs,
                                     config=config, logger=logger)
        self.assertAlmostEqual(percent,  98.088, 2)
        self.assertTrue(ID_OK)
        self.to_be_removed.append(
            os.path.join(self.test_dir, "id_check_AP017923.1", ""))

    @unittest.skipIf(shutil.which("blastn") is None,
                     "Cannot test check_id without blastn")
    def test_check_id_distant(self):
        test_args = Namespace(tmp_dir=self.test_dir,
                              references=[self.ref_fasta],
                              )
        config = Namespace(blastn=shutil.which("blastn"))
        self.assertEqual(
            (False, 0.0),
            bb.check_id(ref=self.ref_fasta, args=test_args, contigs=self.distant_contigs,
                        results=Namespace(current_scaffolds=self.contigs),
                        config=config, logger=logger))
        self.to_be_removed.append(
            os.path.join(self.test_dir, "id_check_AP017923.1", ""))

    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "Skipping this test on Travis CI. Too hard to debug")
    def test_find_origin(self):
        test_args = Namespace(tmp_dir=self.test_dir,
                              references=[self.ref_split],
                              )
        results = Namespace(current_scaffolds=self.contigs,
                            current_reference=self.ref_split)
        config = bb.parse_config(self.active_config)
        ref_dict = {'AP017923.1': 'NODE_4_length_5704_cov_4.88129:1',
                    'AP017923.2': 'NODE_2_length_10885_cov_4.76832:869'}
        ori = bb.find_origin(
            args=test_args, config=config, results=results,
            ori_dir=self.test_dir, flex=10,
            logger=logger)
        for k, v in sorted(ori.items()):
            self.assertEqual(ref_dict[k], v)


    def test_split_origin(self):
        ref_dict = {'AP017923.1': 'NODE_4_length_5704_cov_4.88129:1',
                    'AP017923.2': 'NODE_2_length_10885_cov_4.76832:869'}
        ori_dir=os.path.join(self.test_dir, "origin", "")
        os.makedirs(ori_dir, exist_ok=True)
        results_path = bb.split_origin(origin_dict=ref_dict,
                                       ori_dir=ori_dir, min_len=5,
                                       scaffolds=self.contigs, logger=logger)
        self.assertEqual(
            os.path.join(ori_dir, "scaffolds_ori_split.fasta"),
            results_path)
        self.assertEqual(
            bb.md5(self.contigs_split_ori),
            bb.md5(results_path))
        self.to_be_removed.append(
            os.path.join(self.test_dir, "origin", ""))

    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "Skipping this test on Travis CI. Too hard to debug")
    def test_find_and_split_origin(self):
        test_args = Namespace(tmp_dir=self.test_dir,
                              references=[self.ref_split],
                              )
        results = Namespace(current_scaffolds=self.contigs,
                            current_reference=self.ref_split)
        config = bb.parse_config(self.active_config)
        ref_dict = {'AP017923.1': 'NODE_4_length_5704_cov_4.88129:1',
                    'AP017923.2': 'NODE_2_length_10885_cov_4.76832:869'}
        bb.find_and_split_origin(args=test_args, config=config, results=results,
                                       tools=None, reads_ns=None,
                                       logger=logger)
        ###################
        #  write a propper test here
        ###################
        self.to_be_removed.append(
            os.path.join(self.test_dir, "origin", ""))

    def test_use_already_assembled(self):
        self.assertTrue(bb.use_already_assembled(
            args=Namespace(already_assembled_dirs=["notempty"])))
        self.assertFalse(bb.use_already_assembled(
            args=Namespace(already_assembled_dirs=[],
                           already_assembled_scaffolds=[],
                           already_assembled_contigs=[])))

    def test_parse_stages_error(self):
        with self.assertRaises(ValueError):
            bb.parse_stages(args=Namespace(stages="za"),
                            reads_ns=Namespace(TRIM=None),
                            logger=logger,
                            all_stages=bb.all_stages)

    def test_parse_stages(self):
        reads_ns = Namespace(TRIM=None)
        bb.parse_stages(args=Namespace(stages="t"),
                        reads_ns=reads_ns, logger=logger,
                        all_stages=bb.all_stages)
        self.assertTrue(reads_ns.TRIM)

    def order_scaffolds():
        pass

    def store_orientation(orientations, contig, or_count):
        pass

    def finish_assembly(args, config, reads_ns, results, logger):
        pass


    def test_get_prokka_cmd(self):
        args = Namespace(genus="A", species="b", strain="c", locustag="LOC",
                         centre="NCBI", threads=7, tmp_dir="./")
        self.assertEqual(
            "prokk --addgenes --outdir pro --prefix prokka --genus A " +
            "--species b --strain c --locustag LOC --centre NCBI " +
              "--cpus 7 seqs.fasta > ./prokka.log 2>&1",
            bb.get_prokka_cmd(exe="prokk", outdir="pro",
                           args=args, seqs="seqs.fasta")
        )


    def test_check_spades_kmers(self):
        cmd = "spades.py -k 22,33,55,77 -o someassembly"
        self.assertEqual(
            "spades.py -k 33,55 -o someassembly",
            bb.check_spades_kmers(assembler="spades", cmd=cmd, readlen=78, min_diff=2, logger=logger))

    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "Skipping this test on Travis CI. Too hard to debug")
    def test_run_scaffolder_sis(self):
        # need all these args because the replace_placeholders function
        test_args = Namespace(
            fastq1=self.fastq1, fastq2=self.fastq2, long_fastq=None,
            # de_fere_contigs=None,
            references=[self.ref_split],
            genome_size=0,
            tmp_dir=self.test_dir,
            scaffolder="sis",
            scaffolder_args=None,
            # memory=8, untrimmed_fastq1=None, untrimmed_fastq2=None, platform="illumina", threads=1
        )
        config = bb.return_config(
            os.path.join(self.test_dir, "integration_config.yaml"), force=True,
            hardfail=False, logger=None)
        tools = Namespace(
            scaffolder=[x for x in config.scaffolders if x['name'].lower() == "sis"][0])
        results = bb.make_empty_results_object()
        results.current_reference = self.ref_split
        results.current_contigs = self.contigs_to_scaf
        results.current_scaffolds = "scaffs"
        reads_ns = bb.assess_reads(args=test_args, config=config,
                                   platform="illumina", logger=logger)
        bb.run_ref_scaffolder(
            args=test_args, config=config, tools=tools, results=results, reads_ns=reads_ns,
            run_id=1, logger=logger)
        self.to_be_removed.append(os.path.join(self.test_dir, "SIS_1"))

    def test_build_agp(self):
        test_args = Namespace(
            fastq1=self.fastq1, fastq2=self.fastq2, long_fastq=None,
            # de_fere_contigs=None,
            references=[self.ref_split],
            genome_size=0,
            tmp_dir=self.test_dir,
            scaffolder="sis",
            scaffolder_args=None,
            mode="draft"
            # memory=8, untrimmed_fastq1=None, untrimmed_fastq2=None, platform="illumina", threads=1
        )
        config = bb.parse_config(self.active_config)
        results = bb.make_empty_results_object()
        results.current_scaffolds=self.scaffold

        reads_ns = bb.assess_reads(args=test_args, config=config,
                                platform="illumina", logger=logger)
        bb.build_agp(
            args=test_args, results=results, reads_ns=reads_ns,
            evidence="paired-ends", logger=logger)
        self.to_be_removed.append(
            os.path.join(self.test_dir, "agp", ""))

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
