#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import mock
import unittest

from Bio.SeqRecord import SeqRecord
from rdflib import Graph, Namespace, XSD

from db_integration import BlazegraphIntegration
from superphy.uploader.sequence_upload import SequenceUploader, SequenceMetadata


__author__ = 'Stephen Kan'
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = 'Stephen Kan'
__email__ = 'stebokan@gmail.com'

n = Namespace("https://github.com/superphy#")
owl = Namespace("http://www.w3.org/2002/07/owl#")
rdf = Namespace("http://www.w3.org/1999/02/22-rdf-syntax-ns#")
rdfs = Namespace("http://www.w3.org/2000/01/rdf-schema#")
gfvo = Namespace("http://www.biointerchange.org/gfvo#")

class SequenceUploaderTestCase(unittest.TestCase):

    """
    @classmethod
    def setUpClass(cls):
        super(SequenceUploaderTestCase, cls).setUpClass()
        cls.setupBlazegraph()

        cls.case = SequenceUploader()

    @classmethod
    def setupBlazegraph(cls):
        g = Graph()
    """

    def setUp(self):
        self.case = SequenceUploader()

    def tearDown(self):
        del self.case

    def seqdata(self):
        return SequenceMetadata("AAAA","AAAA")

    def test_load_sequence(self):
        pass

    def test_validated_upload(self):
        pass

    def test_error_logging(self):
        pass

    @mock.patch('superphy.uploader.sequence_upload.SequenceUploader.from_ftp')
    @mock.patch('superphy.uploader.sequence_upload.SequenceUploader.from_nuccore')
    def test_get_seqdata(self, mock_nuccore, mock_ftp):
        mock_nuccore.side_effect = [ValueError, None]
        mock_ftp.side_effect = None

        self.case.get_seqdata(self.seqdata())
        mock_nuccore.assert_called_once_with(mock.ANY)
        mock_ftp.assert_called_once_with(mock.ANY)

        mock_nuccore.reset_mock()
        mock_ftp.reset_mock()
        self.case.get_seqdata(self.seqdata())
        mock_nuccore.assert_called_once_with(mock.ANY)
        mock_ftp.assert_not_called()

    @mock.patch('superphy.uploader.sequence_upload.SequenceUploader.read_fasta')
    @mock.patch('superphy.uploader.sequence_upload.Entrez', autospec=True)
    def test_from_nuccore(self, mock_entrez, mock_read):
        mock_entrez.efetch.return_value = mock.MagicMock(spec=file)

        seqdata = self.seqdata()
        seqdata.dict["sequences"] = ['']

        with self.assertRaises(ValueError):
            self.case.from_nuccore(seqdata)
        self.assertEqual(mock_entrez.email, "superphy.info@gmail.com")
        self.assertEqual(mock_entrez.efetch.call_count, 1)
        self.assertEqual(mock_read.call_count, 1)

        mock_entrez.reset_mock()
        mock_read.reset_mock()
        seqdata = self.seqdata()
        seqdata.dict["sequences"] = ['asdf']
        seqdata.dict["is_from"] = "PLASMID"
        self.case.from_nuccore(seqdata)
        self.assertEqual(seqdata.dict["is_from"], "PLASMID")

        mock_entrez.reset_mock()
        mock_read.reset_mock()
        seqdata = self.seqdata()
        seqdata.dict["sequences"] = ['asdf']
        self.case.from_nuccore(seqdata)
        self.assertEqual(seqdata.dict["is_from"], "CORE")

    @mock.patch('superphy.uploader.sequence_upload.SequenceUploader.read_fasta')
    @mock.patch('superphy.uploader.sequence_upload.open')
    @mock.patch('superphy.uploader.sequence_upload.SequenceUploader.download_file')
    def test_from_ftp(self, mock_download, mock_open, mock_read):
        mock_open.return_value = mock.MagicMock(spec=file)
        seqdata = self.seqdata()

        self.case.from_ftp(seqdata)

        mock_download.assert_called_once_with(mock.ANY, mock.ANY)
        mock_open.assert_called_once_with(mock.ANY, mock.ANY)
        mock_read.assert_called_once_with(mock.ANY, mock.ANY)
        self.assertEqual(seqdata.dict["is_from"], "WGS")

    @mock.patch('superphy.uploader.sequence_upload.SeqIO.parse')
    def test_read_fasta(self, mock_SeqIO):
        record1 = SeqRecord(seq="abcd")
        record2 = SeqRecord(seq="efghi")
        record3 = SeqRecord(seq="jklmno")
        plasmid = SeqRecord(seq="abcd", description="PLASMID")
        mock_handle = mock.MagicMock(spec=file)

        mock_SeqIO.return_value = [record3, record1, record2]

        seqdata = self.seqdata()
        self.case.read_fasta(mock_handle, seqdata)

        self.assertEqual(seqdata.dict["sequences"],["abcd", "efghi", 'jklmno'])
        self.assertIsNone(seqdata.dict["is_from"])

        mock_SeqIO.return_value = [plasmid]
        seqdata = self.seqdata()
        self.case.read_fasta(mock_handle, seqdata)

        self.assertEqual(seqdata.dict["sequences"],["abcd"])
        self.assertEqual(seqdata.dict["is_from"], "PLASMID")

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
        pass

    def test_generate_checksum(self):
        pass

    def test_generate_kwargs(self):
        pass


if __name__ == '__main__':
    unittest.main()
