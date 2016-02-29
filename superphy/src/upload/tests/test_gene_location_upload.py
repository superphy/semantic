#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import unittest
import mock

from collections import defaultdict
import sys
import traceback
import subprocess
import json
import re
import os

from rdflib import Graph
from Bio.Blast import NCBIXML

from superphy.upload._sparql import check_NamedIndividual, has_ref_gene, _sparql_query
from superphy.upload._utils import generate_output, generate_path
from superphy.upload.classes import GeneLocation
from superphy.upload.blazegraph_upload import BlazegraphUploader
from superphy.upload.contig_upload import ContigUploader
from superphy.upload.gene_location_upload import GeneLocationUploader

class GeneLocationUploaderTestCase(unittest.TestCase):
    def setUp(self):
        self.case = GeneLocationUploader()

    def tearDown(self):
        del self.case

    def test_accession_name(self):
        contig = "ANVW00000000"
        desc = "a description"

        returned_string = self.case.accession_name(contig, desc)
        self.assertEqual(returned_string, "ANVW00000000")

        contig = "CP102000"
        desc = "a complete genome"
        returned_string = self.case.accession_name(contig, desc)
        self.assertEqual(returned_string, "CP102000_closed")

    @mock.patch('superphy.upload.gene_location_upload.GeneLocationUploader.get_num_gene_copies')
    def test_add_contig(self, mock_copies):
        #Adding a location to an existing gene with existing contig
        self.case.dict = {"hlyA":{"ANVW01000001":0}}
        self.case.add_contig("hlyA", "ANVW01000001")
        self.assertEqual(self.case.dict, {"hlyA":{"ANVW01000001":1}})

        #Adding location to an existing gene with non-existing contig
        mock_copies.return_value = 0
        self.case.add_contig("hlyA", "JPQG01000001")
        self.assertEqual(self.case.dict, {"hlyA":{"ANVW01000001":1, "JPQG01000001": 0}})

        #Adding location w/ non-existing genome and contig
        self.case.add_contig("espF", "JPQG01000002")
        self.assertEqual(self.case.dict, {"hlyA":{"ANVW01000001":1, "JPQG01000001": 0}, 
                                          "espF":{"JPQG01000002": 0}})

        #Adding location w/ exisitng genome and contig in Blazegraph (but not the dict)
        mock_copies.reset_mock()
        mock_copies.return_value = 4
        self.case.add_contig("aafA", "JPQG01000001")
        self.assertEqual(self.case.dict, {"hlyA":{"ANVW01000001":1, "JPQG01000001": 0}, 
                                          "espF":{"JPQG01000002": 0},
                                          "aafA":{"JPQG01000001": 4}})

    def test_get_gene_name(self):
        string = "AF332513.1.gene1 blaTEM-63. blaTEM-63 encodes extended spectrum beta-lactamase TEM-63. \
                  ARO:1000001 process or component of antibiotic biology or chemistry. ARO:3000931 TEM-63. \
                  ARO:3000873 TEM-1. [Escherichia coli]"
        gene_name = self.case.get_gene_name(string)
        self.assertEqual(gene_name, "blaTEM-63")


    def test_get_num_gene_copies(self):
        pass

    def test_ncbixml_parse(self):
        pass

    def test_is_complete_genome(self):
        pass

    def test_add_to_graph(self):
        pass

    # @mock.patch('superphy.upload.classes.GeneLocation.rdf')
    # @mock.patch('superphy.upload.contig_upload.BlazegraphUploader', autospec=True)
    # @mock.patch('superphy.upload.gene_location_upload.GeneLocationUploader.check_gene_copy')
    # def test_create_gene_location(self, mock_check, mock_location, mock_bg):
    #     mock_check.return_value = False
    #     self.case.create_gene_location("name", "gene", "contig", "begin", "end", "seq", "ref_gene")

    #     #self.assertEqual(len(mock_genelocation.mock_calls), 1)
    #     self.assertEqual(len(mock_check.mock_calls), 1)
    #     self.assertEqual(len(mock_location.mock_calls), 1)
    #     mock_bg.assert_called_once_with()


    def test_check_gene_copy(self):
        pass

    def test_get_reference_genes(self):
        pass

    def test_parse_result(self):
        pass