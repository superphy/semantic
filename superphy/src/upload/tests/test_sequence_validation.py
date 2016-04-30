#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import unittest
import mock

from Bio.Blast import Record

from superphy.upload.contig_upload import ContigsWrapper
from superphy.upload.sequence_validation import SequenceValidator

__author__ = 'Stephen Kan'
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = 'Stephen Kan'
__email__ = 'stebokan@gmail.com'


class SequenceValidatorTestCase(unittest.TestCase):
    def setUp(self):
        self.seqdata = ContigsWrapper("ABCD00000000", "ABCD00000000")

    def tearDown(self):
        del self.seqdata

    @mock.patch("superphy.upload.sequence_validation.open")
    @mock.patch('superphy.upload.sequence_validation.check_checksum')
    @mock.patch('superphy.upload.sequence_validation.SequenceValidator.check_chars')
    @mock.patch('superphy.upload.sequence_validation.SequenceValidator.check_contigs')
    @mock.patch('superphy.upload.sequence_validation.SequenceValidator.check_bp')
    @mock.patch('superphy.upload.sequence_validation.SequenceValidator.check_hits')
    @mock.patch('superphy.upload.sequence_validation.SequenceValidator.filter_passing_hits')
    def test_validate(self, mock_filter, mock_hits, mock_bp, mock_contigs, mock_chars, mock_checksum, mock_open):
        mock_open.return_value = mock.MagicMock(spec=file)
        mock_hits.return_value = True
        mock_bp.return_value = True
        mock_contigs.return_value = True
        mock_chars.return_value = True
        mock_checksum.return_value = False

        SequenceValidator(self.seqdata).validate()
        self.assertTrue(self.seqdata.valid)

        mock_open.reset_mock()
        mock_contigs.return_value = False
        mock_chars.return_value = False
        SequenceValidator(self.seqdata).validate()
        self.assertFalse(self.seqdata.valid)
        self.assertEqual(mock_open.call_count, 2)

    def test_check_hits(self):
        test_set = [2, 11]

        for x in test_set:
            self.seqdata.hits = [n for n in range(x)]
            self.assertFalse(SequenceValidator(self.seqdata).check_hits())

        test_set = [3, 4, 9, 10]

        for x in test_set:
            self.seqdata.hits = [n for n in range(x)]
            self.assertTrue(SequenceValidator(self.seqdata).check_hits())

    def test_check_bp(self):
        test_set = [0, 1000, 3499999, 7500001]

        for x in test_set:
            self.seqdata.basepairs = x
            self.assertFalse(SequenceValidator(self.seqdata).check_bp())

        test_set = [3500000, 3500001, 5000000, 7490000, 7500000]

        for x in test_set:
            self.seqdata.basepairs = x
            self.assertTrue(SequenceValidator(self.seqdata).check_bp())

    def test_check_contigs(self):
        test_set = [-1, 0, 10001]

        for x in test_set:
            self.seqdata.numcontigs = x
            self.assertFalse(SequenceValidator(self.seqdata).check_contigs())

        test_set = [1, 2, 5000, 9999, 10000]

        for x in test_set:
            self.seqdata.numcontigs = x
            self.assertTrue(SequenceValidator(self.seqdata).check_contigs())

    def test_check_chars(self):
        self.seqdata.contigs = [("ANVW", "AAAAGCGCCTTTAGGGCGCTTTTTTACATTGGTGGGTCGTGCAGGATTCGAACCTGCGACCAATTGATTA"),
                                 ("JHNV", "AAAGTCAACTGCTCTACCAACTGAGCTAACGACCCGAAGTGGTGGGTGATGACGGGATCGAACCGCCGAC")]
        self.assertTrue(SequenceValidator(self.seqdata).check_chars())

        self.seqdata.contigs = [("ANVW", "AAAAGCGCCTTTAGGGCGCTTTTTTACATTGGTGGGTCGTGCAGGATTCGAACCTGCGACCAATTGATTA"),
                                ("JHNV", "AAAGTCAACTGCTCTACCAACTGAGCTAACGACCCGAAGTGGTGGGTGATGACGGGATCGAACCGC3@#CGAC")]
        self.assertFalse(SequenceValidator(self.seqdata).check_chars())

    @mock.patch('superphy.upload.sequence_validation.NCBIXML.parse', autospec=True)
    @mock.patch('superphy.upload.sequence_validation.open')
    @mock.patch('superphy.upload.sequence_validation.SequenceValidator.blastn_commandline')
    @mock.patch('superphy.upload.sequence_validation.SequenceValidator.create_fasta')
    def test_filter_passing_hits(self, mock_fasta, mock_blast, mock_open, mock_parse):
        mock_open.return_value = mock.MagicMock(spec=file)

        tests = [x * 50  for x in range(15,18)]

        for n in tests:
            mock_parse.return_value = [self.create_sample_record("Test Region A", n, 1000),
                                       self.create_sample_record("Test Region A", n, 1000),
                                       self.create_sample_record("Test Region B", n, 1000),
                                       self.create_sample_record("Test Region B", n, 1000),
                                       self.create_sample_record("Test Region C", n, 1000)]
            SequenceValidator(self.seqdata).filter_passing_hits()
            self.assertFalse(self.seqdata.hits)

        tests = [x * 50  for x in range(18,21)]
        for n in tests:
            mock_parse.return_value = [self.create_sample_record("Test Region A", n, 1000),
                                       self.create_sample_record("Test Region A", n, 1000),
                                       self.create_sample_record("Test Region B", n, 1000),
                                       self.create_sample_record("Test Region B", n, 1000),
                                       self.create_sample_record("Test Region C", n, 1000)]
            SequenceValidator(self.seqdata).filter_passing_hits()
            self.assertTrue(self.seqdata.hits)

        tests = [x * 50  for x in range(21, 23)]
        for n in tests:
            mock_parse.return_value = [self.create_sample_record("Test Region A", n, 1000),
                                       self.create_sample_record("Test Region A", n, 1000),
                                       self.create_sample_record("Test Region B", n, 1000),
                                       self.create_sample_record("Test Region B", n, 1000),
                                       self.create_sample_record("Test Region C", n, 1000)]
            SequenceValidator(self.seqdata).filter_passing_hits()
            self.assertFalse(self.seqdata.hits)

    def create_sample_record(self, hit_def, positives, length):
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
