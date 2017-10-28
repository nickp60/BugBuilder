# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:57:31 2016
@author: nicholas
"""
import sys
import logging
import os
import unittest

from .context import BugBuilder
from BugBuilder.BugBuilder import parse_config

logger = logging


@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class test_BugBuilder(unittest.TestCase):
    """ tests for riboSeed.py
    """
    def setUp(self):
        self.ref_dir = os.path.join(os.path.dirname(__file__), "references")

        self.empty_config = os.path.join(self.ref_dir, "empty_config.yaml")
        pass

    def test_parse_config(self):
        """ test pandas import
        """
        print(parse_config(self.empty_config))


    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
