#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import unittest
import mock

from ijson.backends import YAJLImportError
try:
    import ijson.backends.yajl2 as ijson
except YAJLImportError:
    import ijson.backends.yajl as ijson

from superphy.upload.metadata_upload import MetadataUploader, Metadata, GenomeMetadataUploader, GenomeMetadata, GeneMetadataUploader, GeneMetadata
from superphy.upload._utils import generate_path

__author__ = "Stephen Kan"
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Stephen Kan"
__email__ = "stebokan@gmail.com"


class MetadataUploaderTestCase(unittest.TestCase):
    def setUp(self):
        self.case = MetadataUploader('samples/2_genome.json')
        self.metadata = Metadata("JHNI00000000")

    def tearDown(self):
        del self.case

    @mock.patch('superphy.upload.metadata_upload.open')
    def test_error_logging(self, mock_open):
        mock_open.return_value = mock.MagicMock(spec=file)

        self.case.error_logging(self.metadata)
        self.assertEqual(self.case.error, 1)


class MetadataTestCase(unittest.TestCase):
    def setUp(self):
        self.case = Metadata("JHNI00000000")

    def test_add_parameter(self):
        self.case.add_parameter("isolation_host", "btaurus")
        self.assertEqual(self.case.dict["isolation_host"], {"btaurus"})

    def test_build_kwargs(self):
        self.case.add_parameter("serotype","O157:H7")
        self.case.add_parameter("strain", "CS03")
        self.case.add_parameter("syndrome", "Urinary tract infection (cystitis)")

        kwargs = self.case.build_kwargs()
        expected = {'strain': {'CS03'}, 'name': 'JHNI00000000', 'serotype': set(['O157:H7']), 'syndrome': {'Urinary tract infection (cystitis)'}}

        print kwargs
        self.assertEqual(kwargs, expected)


class GenomeMetadataUploaderTestCase(unittest.TestCase):
    def setUp(self):
        self.case = GenomeMetadataUploader('samples/2_genome.json', 'ecoli')
        self.metadata = GenomeMetadata("ANVW00000000", "ecoli")

    def tearDown(self):
        del self.case
        del self.metadata

    @mock.patch('superphy.upload.metadata_upload.ijson.parse')
    @mock.patch('superphy.upload.metadata_upload.GenomeMetadataUploader.parse_metadata')

    def test_upload(self, mock_parse, mock_ijson):
        self.case.upload()
        mock_ijson.assert_called_with(mock.ANY)
        mock_parse.assert_called_with(mock.ANY)

    @mock.patch('superphy.upload.metadata_upload.GenomeMetadataUploader.add_to_graph')
    def test_parse_metadata(self, mock_add):

        def side_effect(metadata):
            if metadata:
                self.assertTrue("accession" in metadata.dict.keys())

        mock_add.side_effect = side_effect

        with open(generate_path("../python/" + self.case.filename), "r") as fd:
            data = ijson.parse(fd)
            self.case.parse_metadata(data)

    @mock.patch('superphy.upload.metadata_upload.GenomeMetadataUploader.error_logging')
    @mock.patch('superphy.upload.metadata_upload.GenomeMetadataUploader.create_pending_genome')
    def test_add_to_graph(self, mock_create, mock_error):
        mock_create.side_effect = TypeError
        self.case.add_to_graph(self.metadata)
        self.case.add_to_graph(self.metadata)
        mock_create.assert_called_with(mock.ANY)
        mock_error.assert_called_with(mock.ANY)

    @mock.patch('superphy.upload.metadata_upload.BlazegraphUploader.upload_data')
    @mock.patch('superphy.upload.metadata_upload.GenomeMetadataUploader.get_ncbi_ids')
    @mock.patch('superphy.upload.metadata_upload.check_NamedIndividual')
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

    @mock.patch('superphy.upload.metadata_upload.return_elink_uid')
    @mock.patch('superphy.upload.metadata_upload.return_esearch_uid')
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
    
    def test_get_serotypes(self):
        self.assertEqual(self.case.get_serotypes({"ONT:NM"}), {"Otype": "NT", "Htype": "NM"})
        self.assertEqual(self.case.get_serotypes({"O157:NA"}), {"Otype": "157", "Htype": None})
        self.assertEqual(self.case.get_serotypes({"O157:H7"}), {"Otype": "157", "Htype": "7"})


class GeneMetadataUploaderTestCase(unittest.TestCase):
    def setUp(self):
        self.case = GeneMetadataUploader('data/superphy_vf.json', "virulence_factor")
        self.metadata = GeneMetadata("hlyA")

    def tearDown(self):
        del self.case
        del self.metadata

    @mock.patch('superphy.upload.metadata_upload.GeneMetadataUploader.parse_amr')
    @mock.patch('superphy.upload.metadata_upload.GeneMetadataUploader.parse_vf')
    @mock.patch('superphy.upload.metadata_upload.json.load')
    @mock.patch('superphy.upload.metadata_upload.open')
    def test_upload_genes(self, mock_open, mock_load, mock_vf, mock_amr):
        mock_open.return_value = mock.MagicMock(spec=file)

        self.case.upload_genes()
        mock_load.assert_called_with(mock.ANY)
        self.assertEqual(mock_open.call_count, 1)
        self.assertEqual(mock_vf.call_count, 1)
        self.assertEqual(mock_amr.call_count, 0)

    
    @mock.patch('superphy.upload.metadata_upload.GeneMetadataUploader.add_to_graph')
    @mock.patch('superphy.upload.metadata_upload.GeneMetadata')
    def test_parse_vf(self, mock_gene, mock_add):
        data = {"hlyA": {"vfo_id": "301222", "category": "adherence"}}
        self.case.parse_vf(data)

        mock_gene.assert_called_with("hlyA")
        self.assertEqual(len(mock_add.mock_calls), 1)

    @mock.patch('superphy.upload.metadata_upload.GeneMetadataUploader.add_to_graph')
    @mock.patch('superphy.upload.metadata_upload.GeneMetadata')
    def test_parse_amr(self, mock_gene, mock_add):
        data = {"1": {"ARO_name": "CblA-1", "ARO_category": {
                  "36268": {
                    "category_aro_accession": "3000129",
                    "category_aro_cvterm_id": "36268",
                    "category_aro_name": "beta-lactam resistance gene",
                    "category_aro_description": "Genes conferring resistance to beta-lactams."
                  },
                  "36696": {
                    "category_aro_accession": "3000557",
                    "category_aro_cvterm_id": "36696",
                    "category_aro_name": "antibiotic inactivation enzyme",
                    "category_aro_description": "Enzyme that catalyzes the inactivation of an antibiotic resulting in resistance.  Inactivation includes chemical modification, destruction, etc."
                  }
                }}}
        self.case.parse_amr(data)

        mock_gene.assert_called_with("CblA-1")
        self.assertEqual(len(mock_add.mock_calls), 1)


    def test_remove_bad_chars(self):
        s = "hlyA/yexT"
        s2 = self.case.remove_bad_chars(s)
        self.assertEqual(s2, "hlyA_yexT")

        s = "some gene name like AAC(6'')"
        s2 = self.case.remove_bad_chars(s)
        self.assertEqual(s2, "some_gene_name_like_AAC_6_prime")

    @mock.patch('superphy.upload.metadata_upload.GeneMetadataUploader.error_logging')
    @mock.patch('superphy.upload.metadata_upload.GeneMetadataUploader.create_gene')
    def test_add_to_graph(self, mock_create, mock_error):
        mock_create.side_effect = TypeError
        self.case.add_to_graph(self.metadata)
        self.case.add_to_graph(self.metadata)
        mock_create.assert_called_with(mock.ANY)
        mock_error.assert_called_with(mock.ANY)


    @mock.patch('superphy.upload.metadata_upload.BlazegraphUploader.upload_data')
    @mock.patch('superphy.upload.metadata_upload.check_NamedIndividual')
    def test_create_gene(self, mock_check, mock_upload):
        mock_check.side_effect = [False, True]
        mock_upload.side_effect = ValueError('End of function')

        try:
            self.case.create_gene(self.metadata)
            self.fail("Did not execute mock_upload and throw a ValueError")
        except ValueError:
            print "Uploading!"

        try:
            self.case.create_gene(self.metadata)
        except ValueError:
            self.fail("Did execute mock_upload and throw a ValueError")


class GeneMetadataTestCase(unittest.TestCase):
    def setUp(self):
        self.case = GeneMetadata("hlyA")

    def tearDown(self):
        del self.case

if __name__ == '__main__':
    unittest.main()
