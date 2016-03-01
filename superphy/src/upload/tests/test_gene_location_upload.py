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

from rdflib import Graph, Namespace, Literal, XSD
from Bio.Blast import NCBIXML

from superphy.upload._sparql import check_NamedIndividual, has_ref_gene, _sparql_query
from superphy.upload._utils import generate_output, generate_path
from superphy.upload.classes import GeneLocation
from superphy.upload.blazegraph_upload import BlazegraphUploader
from superphy.upload.contig_upload import ContigUploader
from superphy.upload.gene_location_upload import GeneLocationUploader

n = Namespace("https://github.com/superphy#")
owl = Namespace("http://www.w3.org/2002/07/owl#")
rdf = Namespace("http://www.w3.org/1999/02/22-rdf-syntax-ns#")
rdfs = Namespace("http://www.w3.org/2000/01/rdf-schema#")
gfvo = Namespace("http://www.biointerchange.org/gfvo#")
faldo = Namespace("http://biohackathon.org/resource/faldo#")

class GeneLocationUploaderTestCase(unittest.TestCase):
    def setUp(self):
        self.case = GeneLocationUploader()

    @classmethod
    def setUpClass(cls):
        super(GeneLocationUploaderTestCase, cls).setUpClass()
        cls.setupBlazegraph()

    @classmethod
    def setupBlazegraph(cls):
        """
        Assumes that Blazegraph is initially empty when testing.
        """
        g = Graph()
        g.add((n.agn43_AP009048_closed_0, rdf.type, faldo.Region))
        g.add((n.agn43_AP009048_closed_0, rdf.type, n.reference_gene))
        g.add((n.agn43_AP009048_closed_0, n.is_gene_of, n.AP009048))
        g.add((n.agn43_AP009048_closed_0, n.has_sequence, Literal("ACGTTGCA", datatype=XSD.string)))
        g.add((n.agn43, n.has_copy, n.agn43_AP009048_closed_0))
        g.add((n.agn43_AP009048_closed_0_begin, faldo.position, Literal("2073676", datatype=XSD.string)))
        g.add((n.agn43_AP009048_closed_0_end, faldo.position, Literal("2076795", datatype=XSD.string)))
        g.add((n.agn43_AP009048_closed_0, faldo.begin, n.agn43_AP009048_closed_0_begin))
        g.add((n.agn43_AP009048_closed_0, faldo.end, n.agn43_AP009048_closed_0_end))
        g.add((n.agn43, n.has_copy, n.agn43_AP009048_closed_0))
        g.add((n.AP009048, n.has_gene, n.agn43_AP009048_closed_0))

        g.add((n.ecpC_CP002729_closed_0, rdf.type, faldo.Region))
        g.add((n.ecpC_CP002729_closed_0, n.has_sequence, Literal("ACGTTGCA", datatype=XSD.string)))
        g.add((n.ecpC, n.has_copy, n.ecpC_CP002729_closed_0))
        g.add((n.CP002729, n.has_gene, n.ecpC_CP002729_closed_0))

        BlazegraphUploader().upload_data(generate_output(g))

        del g

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
        # AMR
        string = "AF332513.1.gene1 blaTEM-63. blaTEM-63 encodes extended spectrum beta-lactamase TEM-63. \
                  ARO:1000001 process or component of antibiotic biology or chemistry. ARO:3000931 TEM-63. \
                  ARO:3000873 TEM-1. [Escherichia coli]"
        gene_name = self.case.get_gene_name(string)
        self.assertEqual(gene_name, "blaTEM-63")

        #VF
        string = "agn43|VFO:3000001| - (b2000) - CP4-44 prophage; antigen 43 (Ag43) phase-variable biofilm \
                  formation autotransporter [Escherichia coli str. MG1655 (K12)]"
        gene_name = self.case.get_gene_name(string)
        self.assertEqual(gene_name, "agn43")

    
    def test_get_num_gene_copies(self):
        self.assertEqual(self.case.get_num_gene_copies("agn43", "AP009048"), 1)
        self.assertEqual(self.case.get_num_gene_copies("espF", "ANVW01000001"), 0)
        
    @mock.patch('superphy.upload.gene_location_upload.NCBIXML.parse')
    @mock.patch('superphy.upload.gene_location_upload.open')
    def test_ncbixml_parse(self, mock_open, mock_parse):
        mock_open.return_value = mock.MagicMock(spec=file)
        mock_handle = mock.MagicMock(spec=file)

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
        self.assertTrue(self.case.check_gene_copy("agn43", "AP009048", "2073676", "2076795"))

    def test_get_reference_genes(self):
        self.assertEqual(len(list(self.case.get_reference_genes())), 1)

    @mock.patch('superphy.upload.gene_location_upload.NCBIXML.parse')
    @mock.patch('superphy.upload.gene_location_upload.open')
    def test_parse_result(self, mock_open, mock_parse):
        pass