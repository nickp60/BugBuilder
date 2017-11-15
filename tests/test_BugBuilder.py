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
        self.empty_config = os.path.join(self.ref_dir, "empty_config.yaml")
        self.broken_config = os.path.join(self.ref_dir, "broken_config.yaml")
        self.filled_config = os.path.join(self.ref_dir, "semicomplete_config.yaml")
        self.ref_fasta = os.path.join(self.ref_dir, "AP017923.1.fasta")
        self.ref_split = os.path.join(self.ref_dir, "split_AP017923.1.fasta")
        self.contigs = os.path.join(self.ref_dir, "contigs.fasta")
        self.contigs_to_scaf = os.path.join(self.ref_dir, "contigs_to_scaffold.fasta")
        self.distant_contigs = os.path.join(self.ref_dir, "distant_contigs.fasta")
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
        config = bb.return_config(newpath_config, hardfail=False, logger=logger)
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
            bb.check_and_get_assemblers(args=test_args, config=config,
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
        bb.check_and_get_assemblers(args=test_args, config=config,
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
        reads_ns.read_length_mean = 500
        with self.assertRaises(ValueError):
            bb.check_and_get_assemblers(args=test_args, config=config,
                                reads_ns=reads_ns, logger=logger)

    def test_check_assemblers_need_insertsize(self):
        test_args = Namespace(assemblers=['mascura'],
            fastq1=self.fastq1, fastq2=self.fastq2, long_fastq=None,
            de_fere_contigs=None, reference=None, genome_size=0)
        config = bb.parse_config(self.filled_config)
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
            de_fere_contigs=None, reference=None, genome_size=0)
        config = bb.parse_config(self.filled_config)
        reads_ns = bb.assess_reads(args=test_args, config=config,
                                   platform="illumina", logger=logger)
        test_args.assemblers = ["spades"]
        reads_ns.read_length_mean = 400
        with self.assertRaises(ValueError):
            bb.check_and_get_assemblers(args=test_args, config=config,
                                reads_ns=reads_ns, logger=logger)

    def test_get_scaffolder_and_linkage(self):
        test_args = Namespace(
            scaffolder="sis", reference=self.ref_fasta)
        config = bb.parse_config(self.filled_config)
        bb.get_scaffolder_and_linkage(args=test_args, config=config,
                                      paired=True, logger=logger)

    def test_check_already_assembled_bad_path(self):
        test_args = Namespace(
            already_assembled_dirs=["notAnActualDir"], assemblers=["spades"])
        config = bb.parse_config(self.filled_config)
        with self.assertRaises(FileNotFoundError):
            bb.check_already_assembled_dirs(
                args=test_args, config=config, logger=logger)

    def test_check_already_assembled_missing_file(self):
        test_args = Namespace(
            already_assembled_dirs=[self.ref_dir], assemblers=["spades"])
        config = bb.parse_config(self.filled_config)
        with self.assertRaises(FileNotFoundError):
            bb.check_already_assembled_dirs(
                args=test_args, config=config, logger=logger)

    def test_check_already_assembled_uneven_args(self):
        test_args = Namespace(
            already_assembled_dirs=[self.ref_dir],
            assemblers=["spades", "masurca"])
        config = bb.parse_config(self.filled_config)
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
        config = bb.parse_config(self.filled_config)
        self.assertEqual(
            bb.check_already_assembled_dirs(
                args=test_args, config=config, logger=logger),
            [{"name": "spades",
             "contigs": already_ctg,
              "scaffolds": already_scf}])


    def test_check_args(self):
        config = bb.parse_config(self.filled_config)
        test_args = Namespace(
            merge_method=None, assemblers=["spades", "shovels"])
        with self.assertRaises(ValueError):
            bb.check_args(test_args, config)

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
            "ID_OK",
            "old_contigs",
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

    def check_get_merger_and_linkage(self):
        pass
    def parse_available_platforms():
        pass

    @mock.patch('shutil.which')
    def test_parse_available_assemblers(self, shumock):
        shumock.which = "exectuable"
        self.assertEqual(
            ["abyss", "spades"],
            bb.parse_available("assemblers", self.filled_config)
        )

    def test_parse_available_scaffolders(self):
        pass
    def parse_available_mergers():
        pass
    def parse_available_finisher():
        pass
    def parse_available_varcaller():
        pass
    def get_scaffolder_and_linkage(args, config, paired, logger):
        pass
    def get_merger_and_linkage(args, config, paired):
        pass
    def get_finisher(args, config, paired):
        pass
    def get_varcaller(args, config, paired):
        pass
    def select_tools(args, config, reads_ns, logger):
        pass
    def assembler_needs_downsampling(args, config):
        pass
    def make_fastqc_cmd(args, outdir):
        pass
    def run_fastqc(reads_ns, args, logger=None):
        pass
    def make_sickle_cmd(args, reads_ns, out_dir):
        pass
    def quality_trim_reads(args, config, reads_ns, logger):
        pass
    def make_seqtk_ds_cmd(args, reads_ns, new_coverage, outdir, config, logger):
        pass
    def downsample_reads(args, reads_ns, config, new_cov=100):
        pass
    def make_bwa_cmds(args, config, outdir, ref, reads_ns, fastq1, fastq2):
        pass
    def make_samtools_cmds(config, sam, outdir, sorted_bam):
        pass
    def align_reads(dirname, reads_ns,  downsample, args, config, logger):
        pass
    def make_picard_stats_command(bam, config, picard_outdir):
        pass
    def get_insert_stats(bam, reads_ns, args, config, logger):
        pass
    def replace_placeholders(string, config=None, reads_ns=None, args=None, results=None):
        pass
    def get_assembler_cmds(assembler, assembler_args, args, config, reads_ns):
        pass
    def standardize_fasta_output(infile, outfile, ctype):
        pass
    def get_L50_N50(lengths):
        pass
    def get_contig_stats(contigs, ctype):
        pass
    def run_assembler(assembler, assembler_args, args, reads_ns, config, logger):
        pass
    def merge_assemblies(args, config, reads_ns, logger):
        pass

    @unittest.skipIf(shutil.which("blastn") is None,
                     "Cannot test check_id without blastn")
    def test_check_id_close(self):
        test_args = Namespace(tmp_dir=self.test_dir,
                              reference=self.ref_fasta,
                              )
        config = Namespace(blastn=shutil.which("blastn"))
        self.assertEqual(
            True,
            bb.check_id(args=test_args, contigs=self.contigs,
                        config=config, logger=logger))
        self.to_be_removed.append(
            os.path.join(self.test_dir, "id_check", ""))

    @unittest.skipIf(shutil.which("blastn") is None,
                     "Cannot test check_id without blastn")
    def test_check_id_distant(self):
        test_args = Namespace(tmp_dir=self.test_dir,
                              reference=self.ref_fasta,
                              )
        config = Namespace(blastn=shutil.which("blastn"))
        self.assertEqual(
            False,
            bb.check_id(args=test_args, contigs=self.distant_contigs,
                        config=config, logger=logger))
        self.to_be_removed.append(
            os.path.join(self.test_dir, "id_check", ""))

    def check_already_assembled_dirs(args, config, logger):
        pass
    def check_already_assembled_args(args, config, logger):
        pass
    def make_empty_results_object():
        pass
    def log_read_and_run_data(reads_ns, args, results, logger):  # pragma nocover
        pass
    def run_scaffolder(args, config, reads_ns, results, run_id, logger=None):
        pass
    def find_origin(args, config, results, logger):
        pass
    def check_and_set_trim_length(reads_ns, args, logger):
        pass
    def use_already_assembled(args):
        pass
    def order_scaffolds():
        pass
    def store_orientation(orientations, contig, or_count):
        pass
    def finish_assembly(args, config, reads_ns, results, logger):
        pass
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
            reference=self.ref_split,
            genome_size=0,
            tmp_dir=self.test_dir,
            scaffolder="sis",
            scaffolder_args=None,
            # memory=8, untrimmed_fastq1=None, untrimmed_fastq2=None, platform="illumina", threads=1
        )
        config = bb.parse_config(self.filled_config)
        tools = Namespace(
            scaffolder=[x for x in config.scaffolders if x['name'].lower() == "sis"][0])
        results = bb.make_empty_results_object()
        results.current_contigs = self.contigs_to_scaf
        results.current_scaffolds = "scaffs"
        reads_ns = bb.assess_reads(args=test_args, config=config,
                                   platform="illumina", logger=logger)
        bb.run_ref_scaffolder(
            args=test_args, config=config, tools=tools, results=results, reads_ns=reads_ns,
            run_id=1, logger=logger)
        # self.to_be_removed.append(os.path.join(self.test_dir, "SIS_1"))


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
