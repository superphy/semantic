#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import unittest

from superphy.uploader.metadata_validation import MetadataValidator

__author__ = 'Stephen Kan'
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = 'Stephen Kan'
__email__ = 'stebokan@gmail.com'


class MetadataValidatorTestCase(unittest.TestCase):
    def setUp(self):
        self.case = MetadataValidator()

    def tearDown(self):
        pass


    def test_something(self):
        self.assertEqual(True, True)


if __name__ == '__main__':
    unittest.main()
