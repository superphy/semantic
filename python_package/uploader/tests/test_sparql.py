#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import unittest
import os
import subprocess
from superphy.uploader._utils import generate_path
from superphy.uploader import _sparql
from superphy.uploader.metadata_upload import MetadataUploader

__author__ = 'Stephen Kan'
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = 'Stephen Kan'
__email__ = 'stebokan@gmail.com'


class sparqlTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        top_dir = generate_path("../../../")
        os.chdir(top_dir)
        print os.getcwd()
        src = os.path.join(os.getcwd(),"db/bigdata.jnl")
        dst = os.path.join(os.getcwd(),"db/bigdata.jnl.bk")
        subprocess.call("bash bash/kill_port_9999", shell=True)
        subprocess.call("cp %s %s" %(src, dst), shell=True)
        subprocess.call("bash bash/start_blazegraph", shell=True)
        MetadataUploader("python_package/uploader/samples/3_sequence.json", "ecoli").upload()

    @classmethod
    def tearDownClass(cls):
        top_dir = generate_path("../../../")
        os.chdir(top_dir)
        print os.getcwd()
        src = os.path.join(os.getcwd(),"db/bigdata.jnl.bk")
        dst = os.path.join(os.getcwd(),"db/bigdata.jnl")
        subprocess.call("bash bash/kill_port_9999", shell=True)
        subprocess.call("cp %s %s" %(src, dst), shell=True)
        subprocess.call("rm -f %s" % src, shell=True)
        subprocess.call("bash bash/start_blazegraph", shell=True)

    def test_find_from_host(self):
        self.assertEqual(_sparql.find_from_host("Bos taurus (cow)"), "from_btaurus")
        self.assertIsNone(_sparql.find_from_host("Asymptomatic"))
        self.assertIsNone(_sparql.find_from_host("Enteral feeding tube"))
        self.assertIsNone(_sparql.find_from_host("asdlkjasdfjklsdf"))

    def test_find_syndrome(self):
        self.assertEqual(_sparql.find_syndrome("Asymptomatic"), "asymptomatic")
        self.assertIsNone(_sparql.find_syndrome("Enteral feeding tube"))
        self.assertIsNone(_sparql.find_syndrome("Bos taurus (cow)"))
        self.assertIsNone(_sparql.find_syndrome("ADSFKJDF"))

    def test_find_source(self):
        self.assertEqual(_sparql.find_source("Enteral feeding tube"), "enteral_feeding_tube")
        self.assertIsNone(_sparql.find_source("Asymptomatic"))
        self.assertIsNone(_sparql.find_source("Bos taurus (cow)"))
        self.assertIsNone(_sparql.find_source("ADSFKJDF"))

    def test_check_NamedIndividual(self):
        self.assertTrue(_sparql.check_NamedIndividual("hsapiens"))
        self.assertFalse(_sparql.check_NamedIndividual("asdfasdf"))

    def test_find_missing_sequences(self):
        self.assertEqual(len(list(_sparql.find_missing_sequences())), 3)

    def test_find_duplicate_biosamples(self):
        for (biosample, sequences) in _sparql.find_duplicate_biosamples():
            expected = ["CP001165", "CP001164", "CP001163"]
            for sequence in sequences:
                self.assertTrue(str(sequence) in expected)

    def test_find_core_genome(self):
        pass

    def test_delete_instance(self):
        pass

    def test_insert_accession_sequence(self):
        pass

    def test_check_checksum(self):
        pass

if __name__ == '__main__':
    unittest.main()

