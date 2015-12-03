#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import mock
import unittest

from Bio.SeqRecord import SeqRecord

from superphy.uploader.classes import Sequence
from superphy.uploader.sequence_upload import SequenceUploader, SequenceMetadata


__author__ = 'Stephen Kan'
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = 'Stephen Kan'
__email__ = 'stebokan@gmail.com'

class SequenceUploaderTestCase(unittest.TestCase):
    def setUp(self):
        self.case = SequenceUploader()
        self.mock_seqdata = mock.MagicMock(spec=SequenceMetadata)

    def tearDown(self):
        del self.case

    def test_upload_missing_sequences(self):
        pass

    def test_load_sequence(self):
        pass

    def test_upload(self):
        pass

    def test_plasmid_rdf(self):
        seq_ref = mock.MagicMock(spec=Sequence, create=True)

        self.mock_seqdata.accession = "AAAA"
        self.mock_seqdata.genome = "AAAA"
        self.case.plasmid_rdf(self.mock_seqdata, seq_ref)
        self.assertEqual(len(seq_ref.mock_calls), 2)

        seq_ref.reset_mock()
        self.mock_seqdata.accession = "AAAA"
        self.mock_seqdata.genome = "BBBB"
        self.case.plasmid_rdf(self.mock_seqdata, seq_ref)
        self.assertEqual(len(seq_ref.mock_calls), 1)

    def test_nonplasmid_rdf(self):
        pass

    def test_error_logging(self):
        pass

    @mock.patch('superphy.uploader.sequence_upload.SequenceUploader.from_ftp')
    @mock.patch('superphy.uploader.sequence_upload.SequenceUploader.from_nuccore')
    def test_get_seqdata(self, mock_nuccore, mock_ftp):
        mock_nuccore.side_effect = [ValueError, None]
        mock_ftp.side_effect = None

        self.case.get_seqdata(self.mock_seqdata)
        mock_nuccore.assert_called_once_with(mock.ANY)
        mock_ftp.assert_called_once_with(mock.ANY)

        mock_nuccore.reset_mock()
        mock_ftp.reset_mock()
        self.case.get_seqdata(self.mock_seqdata)
        mock_nuccore.assert_called_once_with(mock.ANY)
        mock_ftp.assert_not_called()

    @mock.patch('superphy.uploader.sequence_upload.SequenceUploader.read_fasta')
    @mock.patch('superphy.uploader.sequence_upload.Entrez', autospec=True)
    def test_from_nuccore(self, mock_entrez, mock_read):
        mock_entrez.efetch.return_value = mock.MagicMock(spec=file)
        self.mock_seqdata.accession = None
        self.mock_seqdata.dict = {}

        self.mock_seqdata.dict["sequences"] = ['']
        with self.assertRaises(ValueError):
            self.case.from_nuccore(self.mock_seqdata)
        self.assertEqual(mock_entrez.email, "superphy.info@gmail.com")
        self.assertEqual(mock_entrez.efetch.call_count, 1)
        self.assertEqual(mock_read.call_count, 1)

        mock_entrez.reset_mock()
        mock_read.reset_mock()
        self.mock_seqdata.dict["sequences"] = ['asdf']
        self.mock_seqdata.dict["is_from"] = "PLASMID"
        self.case.from_nuccore(self.mock_seqdata)
        self.assertEqual(self.mock_seqdata.dict["is_from"], "PLASMID")

        mock_entrez.reset_mock()
        mock_read.reset_mock()
        self.mock_seqdata.dict["is_from"] = None
        self.mock_seqdata.dict["sequences"] = ['asdf']
        self.case.from_nuccore(self.mock_seqdata)
        self.assertEqual(self.mock_seqdata.dict["is_from"], "CORE")

    @mock.patch('superphy.uploader.sequence_upload.SequenceUploader.read_fasta')
    @mock.patch('superphy.uploader.sequence_upload.open')
    @mock.patch('superphy.uploader.sequence_upload.SequenceUploader.download_file')
    def test_from_ftp(self, mock_download, mock_open, mock_read):
        mock_open.return_value = mock.MagicMock(spec=file)
        self.mock_seqdata.accession = None
        self.mock_seqdata.dict = {}

        self.case.from_ftp(self.mock_seqdata)

        mock_download.assert_called_once_with(mock.ANY, mock.ANY)
        mock_open.assert_called_once_with(mock.ANY, mock.ANY)
        mock_read.assert_called_once_with(mock.ANY, mock.ANY)
        self.assertEqual(self.mock_seqdata.dict["is_from"], "WGS")

    @mock.patch('superphy.uploader.sequence_upload.SeqIO.parse')
    def test_read_fasta(self, mock_SeqIO):
        record1 = SeqRecord(seq="abcd")
        record2 = SeqRecord(seq="efghi")
        record3 = SeqRecord(seq="jklmno")
        plasmid = SeqRecord(seq="abcd", description="PLASMID")

        mock_handle = mock.MagicMock(spec=file)
        self.mock_seqdata.dict = {}

        mock_SeqIO.return_value = [record3, record1, record2]

        self.case.read_fasta(mock_handle, self.mock_seqdata)
        print self.mock_seqdata.mock_calls
        (name, args, keywords) = self.mock_seqdata.mock_calls

        print name, args, keywords

        mock_SeqIO.return_value = [plasmid]
        self.case.read_fasta(mock_handle, self.mock_seqdata)

        self.assertEqual(self.mock_seqdata.dict["sequences"],["abcd"])
        self.assertEqual(self.mock_seqdata.dict["is_from"], "PLASMID")

    @mock.patch('superphy.uploader.sequence_upload.open')
    @mock.patch('superphy.uploader.sequence_upload.FTP', autospec=True)
    def test_download_file(self, mock_FTP, mock_open):
        mock_FTP().nlst.return_value = ["wgs.JHNI.1.fsa_nt.gz", "wgs.ABCD.1.fsa_nt.gz", "wgs.IDLE.1.fsa_nt.gz"]
        mock_open.return_value = mock.MagicMock(spec=file)

        self.case.download_file("JHNI", "fsa_nt.gz")
        self.assertEqual(mock_open.call_count, 2)

        mock_open.reset_mock()
        with self.assertRaises(TypeError):
            self.case.download_file("NYAN", "fsa_nt.gz")
        mock_open.assert_not_called()


class SequenceMetadataTestCase(unittest.TestCase):

    def setUp(self):
        self.case = SequenceMetadata("ABCD", "EFGH")

    def tearDown(self):
        del self.case

    def test_add_sequence(self):
        sequences = ["abcd", "efghi", "jklmno", "pqrstuv"]
        self.case.add_sequences(sequences)
        self.assertEqual(self.case.dict["sequences"], sequences)
        self.assertEqual(self.case.dict["contigs"], 4)
        self.assertEqual(self.case.dict["bp"], 22)
        self.assertIsNotNone(self.case.dict["checksum"])

    def test_generate_checksum(self):
        self.case.dict["sequences"] = ["abcd", "efghi", "jklmno", "pqrstuv"]
        self.case.generate_checksum()
        result = "44a66044834cbe55040089cabfc102d5"
        self.assertEqual(self.case.dict["checksum"],result)

    def test_generate_kwargs(self):
        with self.assertRaises(TypeError):
            self.case.generate_kwargs()

        sequences = ["abcd", "efghi", "jklmno", "pqrstuv"]
        self.case.add_sequences(sequences)
        self.case.dict["is_from"] = "CORE"
        results = {"name":"EFGH_seq",
                   "genome":"ABCD",
                   "sequences":sequences,
                   "bp":22,
                   "contigs":4,
                   "checksum":"44a66044834cbe55040089cabfc102d5",
                   "is_from":"CORE"}
        self.assertEqual(self.case.generate_kwargs(), results)


if __name__ == '__main__':
    unittest.main()
