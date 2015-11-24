__author__ = 'Stephen Kan'

import unittest
from superphy.uploader.sequence_validation import SequenceValidator

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
