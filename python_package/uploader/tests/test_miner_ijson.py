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

    def test_error_logging(self):
        self.case.dict["name"] = "ANVW00000000"
        self.case.error_logging()
        self.assertEqual(self.case.error, 1)

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
