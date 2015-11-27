#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import unittest
from superphy.uploader.sequence_validation import SequenceValidator

__author__ = 'Stephen Kan'
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = 'Stephen Kan'
__email__ = 'stebokan@gmail.com'

class SequenceValidatorTestCase(unittest.TestCase):

    def setUp(self):
        self.validSequence = SequenceValidator("ABCD00000000",
                                               ["AAAAGCGCCTTTAGGGCGCTTTTTTACATTGGTGGGTCGTGCAGGATTCGAACCTGCGACCAATTGATTA",
                                                "AAAGTCAACTGCTCTACCAACTGAGCTAACGACCCGAAGTGGTGGGTGATGACGGGATCGAACCGCCGAC"],
                                               140,
                                               2,
                                               "asdfasdf")
        self.invalidSequence = SequenceValidator("ABCD00000000",
                                               ["AAAAGCGCCTTTAGGGCGC**&&$$ACATTGGTGGGTCGTGCAGGATTCGAACCTGCGACCAATTGATTA",
                                                "AAAGTCAACTGCTCTACCAACTGAGCTAACGACCCGAAGTGGTGGGTGATGACGGGATCGAACCGCCGAC"],
                                               140,
                                               2,
                                               "asdfasdf")
    def tearDown(self):
        del self.validSequence
        del self.invalidSequence

    def test_validate(self):
        pass

    def test_check_chars(self):
        self.assertTrue(self.validSequence.check_chars())
        self.assertFalse(self.invalidSequence.check_chars())

    def test_filter_passing_hits(self):
        pass

    def test_blastn_commandline(self):
        pass

    def test_create_fasta(self):
        pass

if __name__ == '__main__':
    unittest.main()
