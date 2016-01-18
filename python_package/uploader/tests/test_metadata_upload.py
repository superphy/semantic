#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import unittest
import mock

from ijson.backends import YAJLImportError
try:
    import ijson.backends.yajl2 as ijson
except YAJLImportError:
    import ijson.backends.yajl as ijson

from superphy.uploader.metadata_upload import MetadataUploader, GenomeMetadata
from superphy.uploader._utils import generate_path

__author__ = "Stephen Kan"
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Stephen Kan"
__email__ = "stebokan@gmail.com"

class MetadataUploaderTestCase(unittest.TestCase):
    def setUp(self):
        self.case = MetadataUploader('samples/2_genome.json', 'ecoli')
        self.metadata = GenomeMetadata("ANVW00000000", "ecoli")

    def tearDown(self):
        del self.case
        del self.metadata

    @mock.patch('superphy.uploader.metadata_upload.ijson.parse')
    @mock.patch('superphy.uploader.metadata_upload.MetadataUploader.parse_metadata')
    def test_upload(self, mock_parse, mock_ijson):
        self.case.upload()
        mock_ijson.assert_called_with(mock.ANY)
        mock_parse.assert_called_with(mock.ANY)

    @mock.patch('superphy.uploader.metadata_upload.MetadataUploader.add_to_graph')
    def test_parse_metadata(self, mock_add):

        def side_effect(metadata):
            if metadata:
                self.assertTrue("accession" in metadata.dict.keys())

        mock_add.side_effect = side_effect

        with open(generate_path("../" + self.case.filename), "r") as fd:
            data = ijson.parse(fd)
            self.case.parse_metadata(data)

    @mock.patch('superphy.uploader.metadata_upload.MetadataUploader.error_logging')
    @mock.patch('superphy.uploader.metadata_upload.MetadataUploader.create_pending_genome')
    def test_add_to_graph(self, mock_create, mock_error):
        mock_create.side_effect = TypeError
        self.case.add_to_graph(self.metadata)
        self.case.add_to_graph(self.metadata)
        mock_create.assert_called_with(mock.ANY)
        mock_error.assert_called_with(mock.ANY)


    @mock.patch('superphy.uploader.metadata_upload.open')
    def test_error_logging(self, mock_open):
        mock_open.return_value = mock.MagicMock(spec=file)

        self.case.error_logging(self.metadata)
        self.assertEqual(self.case.error, 1)


    @mock.patch('superphy.uploader.metadata_upload.BlazegraphUploader.upload_data')
    @mock.patch('superphy.uploader.metadata_upload.MetadataUploader.get_ncbi_ids')
    @mock.patch('superphy.uploader.metadata_upload.check_NamedIndividual')
    def test_create_pending_genome(self, mock_check, mock_ncbi, mock_upload):
        mock_check.side_effect = [False, True]
        mock_ncbi.return_value = {}
        mock_upload.side_effect = ValueError('End of function')

        try:
            self.case.create_pending_genome(self.metadata)
            self.fail("Did not execute mock_upload and throw a ValueError")
        except ValueError:
            print "Uploading!"

        try:
            self.case.create_pending_genome(self.metadata)
        except ValueError:
            self.fail("Did execute mock_upload and throw a ValueError")

    @mock.patch('superphy.uploader.metadata_upload.return_elink_uid')
    @mock.patch('superphy.uploader.metadata_upload.return_esearch_uid')
    def test_get_ncbi_ids(self, mock_esearch, mock_elink):
        def elink_side_effect(*args, **kwargs):
            if args[1] is "bioproject":
                return {"164875"}
            elif args[1] is "biosample":
                return {"854618"}

        mock_esearch.return_value = "431333453"
        mock_elink.side_effect = elink_side_effect

        self.case.get_ncbi_ids(self.metadata)
        self.assertEqual(self.metadata.dict["bioproject"], {"164875"})
        self.assertEqual(self.metadata.dict["biosample"], {"854618"})


class GenomeMetadataTestCase(unittest.TestCase):
    def setUp(self):
        self.case = GenomeMetadata('JHNI00000000', 'ecoli')

    def tearDown(self):
        del self.case

    def test_add_genome_parameter(self):
        self.case.add_genome_parameter("isolation_host", "btaurus")
        self.assertEqual(self.case.dict["isolation_host"], {"btaurus"})

    def test_build_genome_kwargs(self):
        self.case.add_genome_parameter("serotype","O157:H7")
        self.case.add_genome_parameter("strain", "CS03")
        self.case.add_genome_parameter("syndrome", "Urinary tract infection (cystitis)")

        kwargs = self.case.build_genome_kwargs()
        expected = {'strain': {'CS03'}, 'Otype': '157', 'name': 'JHNI00000000', 'Htype': '7', 'organism': 'ecoli',
                    'accession': {'JHNI00000000'}, 'syndrome': {'Urinary tract infection (cystitis)'}}

        self.assertEqual(kwargs, expected)

    def test_get_serotypes(self):
        self.assertEqual(self.case.get_serotypes({"ONT:NM"}), {"Otype": "NT", "Htype": "NM"})
        self.assertEqual(self.case.get_serotypes({"O157:NA"}), {"Otype": "157", "Htype": None})
        self.assertEqual(self.case.get_serotypes({"O157:H7"}), {"Otype": "157", "Htype": "7"})

if __name__ == '__main__':
    unittest.main()
