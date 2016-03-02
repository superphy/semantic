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
from Bio.Blast import NCBIXML, Record

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
        Sets up some data for testing methods in this class that queries Blazegraph.
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
        ## Non-complete genome
        contig = "ANVW00000000"
        desc = "a description"

        returned_string = self.case.accession_name(contig, desc)
        self.assertEqual(returned_string, "ANVW00000000")

        ## Complete genome
        contig = "CP102000"
        desc = "a complete genome"
        returned_string = self.case.accession_name(contig, desc)
        self.assertEqual(returned_string, "CP102000_closed")


    @mock.patch('superphy.upload.gene_location_upload.GeneLocationUploader.get_num_gene_copies')
    def test_add_contig(self, mock_copies):
        ## Adding a location to an existing gene with existing contig
        self.case.dict = {"hlyA":{"ANVW01000001":0}}
        self.case.add_contig("hlyA", "ANVW01000001")
        self.assertEqual(self.case.dict, {"hlyA":{"ANVW01000001":1}})

        ## Adding location to an existing gene with non-existing contig
        mock_copies.return_value = 0
        self.case.add_contig("hlyA", "JPQG01000001")
        self.assertEqual(self.case.dict, {"hlyA":{"ANVW01000001":1, "JPQG01000001": 0}})

        ## Adding location w/ non-existing genome and contig
        self.case.add_contig("espF", "JPQG01000002")
        self.assertEqual(self.case.dict, {"hlyA":{"ANVW01000001":1, "JPQG01000001": 0}, 
                                          "espF":{"JPQG01000002": 0}})

        ## Adding location w/ exisitng genome and contig in Blazegraph (but not the dict)
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

        # VF
        string = "agn43|VFO:3000001| - (b2000) - CP4-44 prophage; antigen 43 (Ag43) phase-variable biofilm \
                  formation autotransporter [Escherichia coli str. MG1655 (K12)]"
        gene_name = self.case.get_gene_name(string)
        self.assertEqual(gene_name, "agn43")

    
    def test_get_num_gene_copies(self):
        self.assertEqual(self.case.get_num_gene_copies("agn43", "AP009048"), 1)
        self.assertEqual(self.case.get_num_gene_copies("espF", "ANVW01000001"), 0)
    

    @mock.patch('superphy.upload.gene_location_upload.GeneLocationUploader.create_gene_location')
    @mock.patch('superphy.upload.gene_location_upload.GeneLocationUploader.check_gene_copy')
    @mock.patch('superphy.upload.gene_location_upload.GeneLocationUploader.get_num_gene_copies')
    @mock.patch('superphy.upload.gene_location_upload.NCBIXML.parse')
    @mock.patch('superphy.upload.gene_location_upload.open')
    def test_ncbixml_parse(self, mock_open, mock_parse, mock_copies, mock_check, mock_create):
        mock_handle = mock.MagicMock(spec=file)
        mock_open.return_value = mock.MagicMock(spec=file)

        ## Gene Location w/ 100% query and identity, incomplete genome
        mock_check.return_value = False
        mock_copies.return_value = 0
        mock_parse.return_value = [self.create_sample_record("aafA", "gnl|BL_ORD_ID|56 gi|606962173|gb|JHNV01000057.1| \
                Escherichia coli O119:H4 str. 03-3458 contig57, whole genome shotgun sequence", 0, 1146, 123, 123, "ATGC")]

        self.case.parse_result()
        mock_create.assert_called_once_with("aafA_JHNV01000057_0", "aafA", "JHNV01000057", '1146', '1268', "ATGC", False)

        ## Gene Location w/ 100% query and identity, complete genome
        mock_check.reset_mock()
        mock_parse.reset_mock()
        mock_create.reset_mock()
        mock_check.return_value = False
        mock_parse.return_value = [self.create_sample_record("bapF", "gnl|BL_ORD_ID|56 gi|606962173|gb|CP002729.1| \
                complete genome", 0, 1146, 1230, 1230, "CAT")]
        self.case.parse_result()
        mock_create.assert_called_once_with("bapF_CP002729_closed_0", "bapF", "CP002729_closed", '1146', '2375', "CAT", False)


    @mock.patch('superphy.upload.gene_location_upload.BlazegraphUploader', autospec=True)
    @mock.patch('superphy.upload.gene_location_upload.Graph')
    @mock.patch('superphy.upload.gene_location_upload.GeneLocation')
    @mock.patch('superphy.upload.gene_location_upload.GeneLocationUploader.check_gene_copy')
    def test_create_gene_location(self, mock_check, mock_rdf, mock_graph, mock_bg):
        mock_check.return_value = False
        mock_genelocation = mock.MagicMock(spec=GeneLocation, create=True)
        mock_g = mock.MagicMock(spec=Graph, create=True)
        mock_graph.return_value = mock_g
        mock_rdf.return_value = mock_genelocation
        self.case.create_gene_location("name", "gene", "contig", "begin", "end", "seq", "ref_gene")

        self.assertEqual(len(mock_check.mock_calls), 1)
        mock_rdf.assert_called_once_with(mock_g, "name", "gene", "contig", "begin", "end", "seq", "ref_gene")
        mock_bg.assert_called_once_with()


    def test_check_gene_copy(self):
        """
        Assumes that there is no data uploaded to Blazegraph before executing these tests.
        """
        self.assertTrue(self.case.check_gene_copy("agn43", "AP009048", "2073676", "2076795"))


    def test_get_reference_genes(self):
        """
        Assumes that there is no data uploaded to Blazegraph before executing these tests.
        """
        self.assertEqual(len(list(self.case.get_reference_genes())), 1)


    @mock.patch('superphy.upload.gene_location_upload.GeneLocationUploader.create_gene_location')
    @mock.patch('superphy.upload.gene_location_upload.NCBIXML.parse', autospec=True)
    @mock.patch('superphy.upload.gene_location_upload.open')
    def test_parse_result(self, mock_open, mock_parse, mock_create):
        mock_open.return_value = mock.MagicMock(spec=file)
        mock_parse.return_value = [self.create_sample_record("gaa", "gnl|BL_ORD_ID|56 gi|606962173|gb|JHNV01000056.1| \
                Escherichia coli O119:H4 str. 03-3458 contig56, whole genome shotgun sequence", 
                                            0, 1146, 123, 123, "ATGC")]

        self.case.parse_result()
        mock_create.assert_called_once_with("gaa_JHNV01000056_0", "gaa", "JHNV01000056", '1146', '1268', "ATGC", False)


    def create_sample_record(self, query, title, expect, start, score, ident, seq):
        """
        Helper function that creates Blast record handles for testing NCBI parse-related methods.
        """
        record = mock.MagicMock(spec=Record)
        entry = mock.MagicMock(spec=Record.Alignment)
        hsp = mock.MagicMock(spec=Record.HSP)

        record.query = query

        entry.title = title
        entry.hsps = [hsp]

        hsp.expect = expect
        hsp.sbjct_start = start
        hsp.score = score
        hsp.identities = ident
        hsp.sbjct = seq

        record.alignments = [entry]
        return record







