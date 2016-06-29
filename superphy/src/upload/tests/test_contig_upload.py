#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import mock
import unittest

from Bio.SeqRecord import SeqRecord

from SuperPhy.models.upload.classes import Contig
from SuperPhy.models.upload.contig_upload import ContigUploader, ContigsWrapper


__author__ = 'Clarice Ng'
__copyright__ = "© Copyright Government of Canada 2012-2016. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = 'Clarice Ng'
__email__ = 'c32ng@uwaterloo.ca'

class ContigUploaderTestCase(unittest.TestCase):
    def setUp(self):
        self.case = ContigUploader()
        self.mock_contigswrapper = mock.MagicMock(spec=ContigsWrapper)

    def tearDown(self):
        del self.case


    @mock.patch('SuperPhy.models.upload.contig_upload.ContigUploader.error_logging')
    @mock.patch('SuperPhy.models.upload.contig_upload.SequenceValidator')
    @mock.patch('SuperPhy.models.upload.contig_upload.ContigUploader.upload')
    @mock.patch('SuperPhy.models.upload.contig_upload.ContigUploader.get_seqdata')
    @mock.patch('SuperPhy.models.upload.contig_upload.ContigsWrapper')
    @mock.patch('SuperPhy.models.upload.contig_upload.find_missing_sequences')
    def test_upload_missing_contigs(self, mock_find, mock_wrapper, mock_data, mock_upload, mock_valid, mock_error):
        self.mock_contigswrapper.dict = {}

        # Test: No genomes with missing contig data
        mock_find.return_value = []
        self.case.upload_missing_contigs()

        self.assertEqual(len(mock_find.mock_calls), 1)
        mock_wrapper.assert_not_called()
        mock_data.assert_not_called()
        mock_upload.assert_not_called()
        mock_valid.assert_not_called()
        mock_error.assert_not_called()

        # Test: One genome returned from find_missing_contig (from Plasmid)
        mock_find.return_value = [("AAAA", "AAAA")]
        self.mock_contigswrapper.dict["is_from"] = "PLASMID"
        mock_wrapper.return_value = self.mock_contigswrapper
        self.case.upload_missing_contigs()

        self.assertEqual(len(mock_find.mock_calls), 2)
        self.assertEqual(len(mock_wrapper.mock_calls), 1)
        self.assertEqual(len(mock_data.mock_calls), 1)
        self.assertEqual(len(mock_upload.mock_calls), 1)
        mock_valid.assert_not_called()
        mock_error.assert_not_called()

        # Test: One genome returned from find_missing_contig (not from Plasmid)
        mock_find.return_value = [("AAAA", "AAAA")]
        self.mock_contigswrapper.dict["is_from"] = "CORE"
        mock_wrapper.return_value = self.mock_contigswrapper
        self.case.upload_missing_contigs()

        self.assertEqual(len(mock_find.mock_calls), 3)
        self.assertEqual(len(mock_wrapper.mock_calls), 2)
        self.assertEqual(len(mock_data.mock_calls), 2)
        self.assertEqual(len(mock_upload.mock_calls), 2)
        self.assertEqual(len(mock_valid.mock_calls), 2)
        mock_error.assert_not_called()

        # Test: One genome returned from find_missing_contig (not from Plasmid), raise type_error in try/except [ERROR]
        mock_find.return_value = [("AAAA", "AAAA")]
        self.mock_contigswrapper.dict["is_from"] = "WGS"
        mock_wrapper.return_value = self.mock_contigswrapper
        mock_data.side_effect = TypeError
        self.case.upload_missing_contigs()

        self.assertEqual(len(mock_find.mock_calls), 4)
        self.assertEqual(len(mock_data.mock_calls), 3)
        self.assertEqual(len(mock_wrapper.mock_calls), 3)
        self.assertEqual(len(mock_upload.mock_calls), 2)
        self.assertEqual(len(mock_valid.mock_calls), 2)
        self.assertEqual(len(mock_error.mock_calls), 1)



    @mock.patch('SuperPhy.models.upload.contig_upload.BlazegraphUploader', autospec=True)
    @mock.patch('SuperPhy.models.upload.contig_upload.Contig')
    def test_upload(self, mock_rdf, mock_bg):
        contig_ref = mock.MagicMock(spec=Contig, create=True)
        mock_rdf.return_value = contig_ref
        sample_func = mock.Mock()
        self.mock_contigswrapper.contigs = [("ANVW01000001", "ACGTTGCA")]

        self.case.upload(self.mock_contigswrapper, sample_func)
        self.assertEqual(self.mock_contigswrapper.generate_kwargs.call_count, 1)
        sample_func.assert_called_once_with(mock.ANY, mock.ANY)
        mock_rdf.assert_called_once_with(mock.ANY)
        mock_bg.assert_called_once_with()

    def test_plasmid_rdf(self):
        contig_ref = mock.MagicMock(spec=Contig, create=True)

        self.mock_contigswrapper.accession = "AAAA"
        self.mock_contigswrapper.genome = "AAAA"
        self.case.plasmid_rdf(self.mock_contigswrapper, contig_ref)
        self.assertEqual(len(contig_ref.mock_calls), 2)

        contig_ref.reset_mock()
        self.mock_contigswrapper.accession = "AAAA"
        self.mock_contigswrapper.genome = "BBBB"
        self.case.plasmid_rdf(self.mock_contigswrapper, contig_ref)
        self.assertEqual(len(contig_ref.mock_calls), 1)

    def test_nonplasmid_rdf(self):
        """
        Mainly tests the if case in the nonplasmid_rdf method
        """
        contig_ref = mock.MagicMock(spec=Contig, create=True)

        self.mock_contigswrapper.valid = True
        self.mock_contigswrapper.hits = []
        self.case.nonplasmid_rdf(self.mock_contigswrapper, contig_ref)
        self.assertEqual(len(contig_ref.mock_calls), 2)

        contig_ref.reset_mock()
        self.mock_contigswrapper.valid = False
        self.mock_contigswrapper.hits = []
        self.case.nonplasmid_rdf(self.mock_contigswrapper, contig_ref)
        self.assertEqual(len(contig_ref.mock_calls), 1)


    def test_error_logging(self):
        pass

    @mock.patch('SuperPhy.models.upload.contig_upload.ContigUploader.download_file')
    @mock.patch('SuperPhy.models.upload.contig_upload.open')
    @mock.patch('SuperPhy.models.upload.contig_upload.ContigUploader.load_contigs')
    @mock.patch('SuperPhy.models.upload.contig_upload.SeqIO.parse')
    @mock.patch('SuperPhy.models.upload.contig_upload.Entrez', autospec=True)
    def test_get_seqdata(self, mock_entrez, mock_SeqIO, mock_load, mock_open, mock_download):
        mock_entrez.efetch.return_value = mock.MagicMock(spec=file)
        mock_handle = mock.MagicMock(spec=file)
        self.mock_contigswrapper.accession = None
        self.mock_contigswrapper.dict = {}
        self.mock_contigswrapper.genome = "BA000007"
        
        #complete genome
        record1 = SeqRecord(seq="abcd", description="complete genome")
        mock_SeqIO.return_value = [record1]
        self.case.get_seqdata(self.mock_contigswrapper)
        
        self.assertEqual(len(mock_load.mock_calls), 1)
        self.assertEqual(self.mock_contigswrapper.dict["is_from"], "CORE")
        

        #non-complete genome
        self.mock_contigswrapper.reset_mock()
        mock_SeqIO.reset_mock()
        record2 = SeqRecord(seq="actg", description="some genome")
        mock_SeqIO.return_value = [record2]
        self.case.get_seqdata(self.mock_contigswrapper)

        mock_open.return_value = mock.MagicMock(spec=file)
        mock_download.assert_called_once_with(mock.ANY, mock.ANY)
        mock_open.assert_called_once_with(mock.ANY, mock.ANY)
        self.assertEqual(self.mock_contigswrapper.dict["is_from"], "WGS")


    @mock.patch('SuperPhy.models.upload.contig_upload.open')
    @mock.patch('SuperPhy.models.upload.contig_upload.FTP', autospec=True)
    def test_download_file(self, mock_FTP, mock_open):
        mock_FTP().nlst.return_value = ["wgs.JHNI.1.fsa_nt.gz", "wgs.ABCD.1.fsa_nt.gz", "wgs.IDLE.1.fsa_nt.gz"]
        mock_open.return_value = mock.MagicMock(spec=file)

        self.case.download_file("JHNI", "fsa_nt.gz")
        self.assertEqual(mock_open.call_count, 2)

        mock_open.reset_mock()
        with self.assertRaises(TypeError):
            self.case.download_file("NYAN", "fsa_nt.gz")
        mock_open.assert_not_called()


class ContigsWrapperTestCase(unittest.TestCase):

    def setUp(self):
        self.case = ContigsWrapper("ANVW00000000", "ANVW01000000")

    def tearDown(self):
        del self.case

    def test_add_contigs(self):
        contigs = [('ANVW01000001', 'ACTGATGGCTGATCGTACGGATTGGCATCGTACTATCGATCAGCTGAGCTACG'),
                   ('ANVW01000002', 'TTTCGAGCTGGACTGATGGCATCGGCAGGACGAGCTTCAGCTACGATCGAAA'),
                   ('ANVW01000003', 'CCTAAGCTGGGATCGATTGAGCTAGCATTGACACACCACACCCATTTTTGT'),
                   ('ANVW01000004', 'GAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTCCATTTTTACTACTAGC')]
        self.case.add_contigs(contigs)
        self.assertEqual(self.case.contigs, contigs)
        self.assertEqual(self.case.numcontigs, 4)
        self.assertEqual(self.case.basepairs, 208)
        self.assertIsNotNone(self.case.checksum)

    def test_generate_checksum(self):
        self.case.contigs = [("a", "atcg"), ("b", "t"), ("c", "actgatgtg")]
        self.case.generate_checksum()
        result = "517e0a5acb941ce868da1e285a8c1372"
        self.assertEqual(self.case.checksum,result)

    def test_generate_kwargs(self):
        with self.assertRaises(TypeError):
            self.case.generate_kwargs()

        sequences = [("a", "atcg"), ("b", "t"), ("c", "actgatgtg")]
        self.case.add_contigs(sequences)
        self.case.dict["is_from"] = "CORE"
        results = {"name":"ANVW01000001",
                   "genome":"ANVW00000000",
                   "sequence":"actactactact",
                   "is_from":"CORE"}
        self.assertEqual(self.case.generate_kwargs("ANVW01000001", "actactactact"), results)


if __name__ == '__main__':
    unittest.main()
