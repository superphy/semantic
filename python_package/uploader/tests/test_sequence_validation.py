#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import unittest
import mock

from Bio.Blast import Record

from superphy.uploader.sequence_upload import SequenceMetadata
from superphy.uploader.sequence_validation import SequenceValidator

__author__ = 'Stephen Kan'
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = 'Stephen Kan'
__email__ = 'stebokan@gmail.com'


class SequenceValidatorTestCase(unittest.TestCase):
    def setUp(self):
        self.sequence = SequenceMetadata("ABCD00000000", "ABCD00000000")

    def tearDown(self):
        del self.sequence

    def test_validate(self):
        pass

    def test_check_hits(self):
        test_set = [2, 11]

        for x in test_set:
            self.sequence.hits = [n for n in range(x)]
            self.assertFalse(SequenceValidator(self.sequence).check_hits())

        test_set = [3, 4, 9, 10]

        for x in test_set:
            self.sequence.hits = [n for n in range(x)]
            self.assertTrue(SequenceValidator(self.sequence).check_hits())

    def test_check_bp(self):
        test_set = [0, 1000, 3499999, 7500001]

        for x in test_set:
            self.sequence.dict["bp"] = x
            self.assertFalse(SequenceValidator(self.sequence).check_bp())

        test_set = [3500000, 3500001, 5000000, 7490000, 7500000]

        for x in test_set:
            self.sequence.dict["bp"] = x
            self.assertTrue(SequenceValidator(self.sequence).check_bp())

    def test_check_contigs(self):
        test_set = [-1, 0, 10001]

        for x in test_set:
            self.sequence.dict["contigs"] = x
            self.assertFalse(SequenceValidator(self.sequence).check_contigs())

        test_set = [1, 2, 5000, 9999, 10000]

        for x in test_set:
            self.sequence.dict["contigs"] = x
            self.assertTrue(SequenceValidator(self.sequence).check_contigs())

    def test_check_chars(self):
        self.sequence.dict["sequences"] = ["AAAAGCGCCTTTAGGGCGCTTTTTTACATTGGTGGGTCGTGCAGGATTCGAACCTGCGACCAATTGATTA",
                                           "AAAGTCAACTGCTCTACCAACTGAGCTAACGACCCGAAGTGGTGGGTGATGACGGGATCGAACCGCCGAC"]
        self.assertTrue(SequenceValidator(self.sequence).check_chars())

        self.sequence.dict["sequences"] = ["AAAAGCGCCTTTAGGGCGC**&&$$ACATTGGTGGGTCGTGCAGGATTCGAACCTGCGACCAATTGATTA",
                                           "AAAGTCAACTGCTCTACCAACTGAGCTAACGACCCGAAGTGGTGGGTGATGACGGGATCGAACCGCCGAC"]
        self.assertFalse(SequenceValidator(self.sequence).check_chars())

    @mock.patch('superphy.uploader.sequence_validation.NCBIXML.parse', autospec=True)
    @mock.patch('superphy.uploader.sequence_validation.open')
    def test_filter_passing_hits(self, mock_open, mock_parse):
        mock_open.return_value = mock.MagicMock(spec=file)

        tests = [x * 50  for x in range(15,18)]

        for n in tests:
            mock_parse.return_value = [self.create_test_record("Test Region A", n, 1000),
                                       self.create_test_record("Test Region A", n, 1000),
                                       self.create_test_record("Test Region B", n, 1000),
                                       self.create_test_record("Test Region B", n, 1000),
                                       self.create_test_record("Test Region C", n, 1000)]
            SequenceValidator(self.sequence).filter_passing_hits()
            self.assertFalse(self.sequence.hits)

        tests = [x * 50  for x in range(18,21)]
        for n in tests:
            mock_parse.return_value = [self.create_test_record("Test Region A", n, 1000),
                                       self.create_test_record("Test Region A", n, 1000),
                                       self.create_test_record("Test Region B", n, 1000),
                                       self.create_test_record("Test Region B", n, 1000),
                                       self.create_test_record("Test Region C", n, 1000)]
            SequenceValidator(self.sequence).filter_passing_hits()
            self.assertTrue(self.sequence.hits)

        tests = [x * 50  for x in range(21, 23)]
        for n in tests:
            mock_parse.return_value = [self.create_test_record("Test Region A", n, 1000),
                                       self.create_test_record("Test Region A", n, 1000),
                                       self.create_test_record("Test Region B", n, 1000),
                                       self.create_test_record("Test Region B", n, 1000),
                                       self.create_test_record("Test Region C", n, 1000)]
            SequenceValidator(self.sequence).filter_passing_hits()
            self.assertFalse(self.sequence.hits)


    def create_test_record(self, hit_def, positives, length):
        record = mock.MagicMock(spec=Record)
        entry = mock.MagicMock(spec=Record.Alignment)
        hsp = mock.MagicMock(spec=Record.HSP)

        hsp.positives = positives

        entry.hit_def = hit_def
        entry.length = length
        entry.hsps = [hsp]

        record.alignments = [entry]
        return record


if __name__ == '__main__':
    unittest.main()
