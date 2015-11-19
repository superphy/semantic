__author__ = 'Stephen Kan'

import unittest
import mock
from superphy.uploader.miner_ijson import MinerDataUploader


class MinerIJSONTestCase(unittest.TestCase):
    def setUp(self):
        self.case = MinerDataUploader('samples/2_genome.json', 'ecoli')

    def tearDown(self):
        self.case = None

    def test_loadJSON(self):
        self.case.load_JSON()
        self.assertTrue(self.case.data)

    def test_start_new_genome(self):
        self.case.start_new_genome("ANVW00000000")
        self.assertEqual(self.case.dict["name"], "ANVW00000000")
        self.assertEqual(self.case.dict["accession"], {"ANVW00000000"})
        self.assertEqual(self.case.dict["organism"], "ecoli")

    def test_add_genome_parameter(self):
        self.case.add_genome_parameter("ANVW00000000.strain.item.displayname", "KTE59")
        self.assertTrue("strain" in self.case.dict)
        self.assertEqual(self.case.dict["strain"], {"KTE59"})

    @mock.patch('superphy.uploader.miner_ijson.open')
    def test_error_logging(self, mock_open):
        mock_open.return_value = mock.MagicMock(spec=file)
        self.case.dict["name"] = "ANVW00000000"
        self.case.error_logging()
        self.assertEqual(self.case.error, 1)

    @mock.patch('superphy.uploader.miner_ijson.upload_data')
    @mock.patch('superphy.uploader.miner_ijson.MinerDataUploader.get_ncbi_ids')
    @mock.patch('superphy.uploader.miner_ijson.check_NamedIndividual')
    def test_create_pending_genome(self, mock_check, mock_ncbi, mock_upload):
        mock_check.side_effect = [False, True]
        mock_ncbi.return_value = {}
        mock_upload.side_effect = ValueError('End of function')

        self.case.dict = {"name":"ANVW00000000"}

        try:
            self.case.create_pending_genome()
            self.fail("Did not execute mock_upload and throw a ValueError")
        except ValueError:
            print "Uploading!"

        try:
            self.case.create_pending_genome()
        except ValueError:
            self.fail("Did execute mock_upload and throw a ValueError")

    @mock.patch('superphy.uploader.miner_ijson.MinerDataUploader.get_ncbi_ids')
    def test_setup_genome_kwargs(self, mock_ncbi):
        mock_ncbi.return_value = {"bioproject":{"164875"},
                                  "biosample":{"854618"}}
        self.case.dict = {"name":"ANVW00000000",
                          "isolation_host":{"Homo sapiens (human)"},
                          "serotype":{"O157:H7"},
                          "strain":{"KTE59"},
                          "isolation_date":{"2010-01-01"},
                          "isolation_location":{"Denmark"},
                          "isolation_source":{"Feces"}}
        kwargs = self.case.setup_genome_kwargs()
        self.assertEqual(kwargs["name"], "ANVW00000000")
        self.assertEqual(kwargs["date"], {"2010-01-01"})
        self.assertEqual(kwargs["host"], {"Homo sapiens (human)"})
        self.assertEqual(kwargs["Otype"], "157")
        self.assertEqual(kwargs["Htype"], "7")

    @mock.patch('superphy.uploader.miner_ijson.return_elink_uid')
    @mock.patch('superphy.uploader.miner_ijson.return_esearch_uid')
    def test_get_ncbi_ids(self, mock_esearch, mock_elink):
        def elink_side_effect(*args, **kwargs):
            if args[1] is "bioproject":
                return {"164875"}
            elif args[1] is "biosample":
                return {"854618"}

        mock_esearch.return_value = "431333453"
        mock_elink.side_effect = elink_side_effect

        result = self.case.get_ncbi_ids("ANVW00000000")
        expected = {"bioproject": {"164875"}, "biosample": {"854618"}}
        self.assertEqual(result, expected)

    def test_get_serotypes(self):
        self.assertEqual(self.case.get_serotypes({"ONT:NM"}), {"Otype": None, "Htype": "-"})
        self.assertEqual(self.case.get_serotypes({"O157:NA"}), {"Otype": "157", "Htype": None})
        self.assertEqual(self.case.get_serotypes({"O157:H7"}), {"Otype": "157", "Htype": "7"})
