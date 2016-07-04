#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import unittest

from rdflib import Graph, Namespace, Literal, XSD, BNode

from db_integration import BlazegraphIntegration
from SuperPhy.models.upload import _sparql
from SuperPhy.models.upload._utils import generate_output
from SuperPhy.models.upload.blazegraph_upload import BlazegraphUploader
from SuperPhy.models.sparql.endpoint import Endpoint

N = Namespace("https://github.com/superphy#")
OWL = Namespace("http://www.w3.org/2002/07/owl#")
RDF = Namespace("http://www.w3.org/1999/02/22-rdf-syntax-ns#")
RDFS = Namespace("http://www.w3.org/2000/01/rdf-schema#")
GFVO = Namespace("http://www.biointerchange.org/gfvo#")

class SPARQLTestCase(BlazegraphIntegration, unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(SPARQLTestCase, cls).setUpClass()
        cls.setupBlazegraph()

    @classmethod
    def setupBlazegraph(cls):
        graph = Graph()

        graph.add((N.btaurus, RDF.type, OWL.NamedIndividual))
        graph.add((N.btaurus, RDFS.label, Literal("Bos taurus (cow)", datatype=XSD.string)))
        graph.add((N.btaurus, N.is_object_of, N.from_btaurus))
        graph.add((N.from_btaurus, RDF.type, N.isolation_from_host))

        graph.add((N.asymptomatic, RDF.type, N.isolation_syndrome))
        graph.add((N.asymptomatic, RDFS.label, Literal("Asymptomatic", datatype=XSD.string)))

        graph.add((N.enteral_feeding_tube, RDF.type, N.isolation_from_source))
        graph.add((N.enteral_feeding_tube, RDFS.label, Literal("Enteral feeding tube", datatype=XSD.string)))

        graph.add((N.CP001165, RDF.type, GFVO.Genome))
        graph.add((N.CP001165, N.has_biosample, Literal("2603441", datatype=XSD.string)))
        graph.add((N.CP001165, N.has_accession, Literal("CP001165", datatype=XSD.string)))

        graph.add((N.CP001164, RDF.type, GFVO.Genome))
        graph.add((N.CP001164, N.has_valid_sequence, Literal("True", datatype=XSD.string)))
        graph.add((N.CP001164, N.has_biosample, Literal("2603441", datatype=XSD.string)))
        graph.add((N.CP001164, N.has_accession, Literal("CP001164", datatype=XSD.string)))
        graph.add((N.CP001164, N.has_sequence, N.CP001164_seq))
        graph.add((N.CP001164_seq, N.is_from, Literal("CORE", datatype=XSD.string)))

        graph.add((N.CP001163, RDF.type, GFVO.Genome))
        graph.add((N.CP001163, N.has_biosample, Literal("2603441", datatype=XSD.string)))
        graph.add((N.CP001163, N.has_accession, Literal("CP001163", datatype=XSD.string)))

        graph.add((N.fakeGenome, RDF.type, GFVO.Genome))
        graph.add((N.fakeGenome, N.has_valid_sequence, Literal("False", datatype=XSD.string)))

        graph.add((N.testEntity, RDF.type, OWL.NamedIndividual))


        graph.add((BNode(), RDFS.label, Literal("Test node", datatype=XSD.string)))

        graph.add((N.test_object, N.has_checksum, Literal("asdfghjkl", datatype=XSD.string)))

        BlazegraphUploader().upload_data(generate_output(graph))

        del graph

    def test_find_from_host(self):
        self.assertEqual(_sparql.find_from_host("Bos taurus (cow)"), "from_btaurus")

        with self.assertRaises(IndexError):
            self.assertIsNone(_sparql.find_from_host("Asymptomatic"))
        with self.assertRaises(IndexError):
            self.assertIsNone(_sparql.find_from_host("Enteral feeding tube"))
        with self.assertRaises(IndexError):
            self.assertIsNone(_sparql.find_from_host("foo"))

    def test_find_syndrome(self):
        self.assertEqual(_sparql.find_syndrome("Asymptomatic"), "asymptomatic")

        with self.assertRaises(IndexError):
            self.assertIsNone(_sparql.find_syndrome("Enteral feeding tube"))
        with self.assertRaises(IndexError):
            self.assertIsNone(_sparql.find_syndrome("Bos taurus (cow)"))
        with self.assertRaises(IndexError):
            self.assertIsNone(_sparql.find_syndrome("foo"))

    def test_find_source(self):
        self.assertEqual(_sparql.find_source("Enteral feeding tube"), "enteral_feeding_tube")

        with self.assertRaises(IndexError):
            self.assertIsNone(_sparql.find_source("Asymptomatic"))
        with self.assertRaises(IndexError):
            self.assertIsNone(_sparql.find_source("Bos taurus (cow)"))
        with self.assertRaises(IndexError):
            self.assertIsNone(_sparql.find_source("foo"))

    def test_check_named_individual(self):
        self.assertTrue(_sparql.check_named_individual("btaurus"))
        self.assertFalse(_sparql.check_named_individual("foo"))

    def test_find_missing_sequences(self):
        self.assertEqual(len(list(_sparql.find_missing_sequences())), 2)

    def test_find_duplicate_biosamples(self):
        for (biosample, sequences) in _sparql.find_duplicate_biosamples():
            expected = ["CP001165", "CP001164", "CP001163"]
            for sequence in sequences:
                self.assertTrue(str(sequence) in expected)
                self.assertFalse(str(sequence) == "fakeGenome")

    def test_find_core_genome(self):
        self.assertEqual(_sparql.find_core_genome("2603441")[0], "CP001164")
        self.assertEqual(_sparql.find_core_genome("456123"), [])

    def test_delete_instance(self):
        self.assertTrue(_sparql.check_named_individual("testEntity"))
        _sparql.delete_instance("testEntity")
        self.assertFalse(_sparql.check_named_individual("testEntity"))

    def test_insert_accession_sequence(self):
        _sparql.insert_accession_sequence("fakeGenome", "fakePlasmid", "fakePlasmidSeq")

        self.assertTrue(Endpoint.ask(
            'PREFIX : <https://github.com/superphy#>\n'
            'ASK { :fakeGenome :has_accession "fakePlasmid"^^xsd:string}'
        ))

        self.assertTrue(Endpoint.ask(
            'PREFIX : <https://github.com/superphy#>\n'
            'ASK { :fakeGenome :has_sequence :fakePlasmidSeq}'
        ))

    def test_blank_nodes(self):
        self.assertTrue(_sparql.check_blank_nodes())
        _sparql.delete_blank_nodes()
        self.assertFalse(_sparql.check_blank_nodes())

    def test_check_checksum(self):
        self.assertTrue(_sparql.check_checksum("asdfghjkl"))
        self.assertFalse(_sparql.check_checksum("123456789"))

if __name__ == '__main__':
    unittest.main()

