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
        config = bb.return_config(newpath_config)
        self.assertEqual(config.STATUS, "COMPLETE")

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
