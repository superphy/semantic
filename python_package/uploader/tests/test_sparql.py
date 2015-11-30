#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import unittest
import os
import subprocess
from rdflib import Graph, Namespace, Literal, XSD, BNode
from superphy.uploader._utils import generate_path, generate_output
from superphy.uploader import _sparql
from superphy.uploader.blazegraph_upload import BlazegraphUploader

__author__ = 'Stephen Kan'
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = 'Stephen Kan'
__email__ = 'stebokan@gmail.com'

n = Namespace("https://github.com/superphy#")
owl = Namespace("http://www.w3.org/2002/07/owl#")
rdf = Namespace("http://www.w3.org/1999/02/22-rdf-syntax-ns#")
rdfs = Namespace("http://www.w3.org/2000/01/rdf-schema#")
gfvo = Namespace("http://www.biointerchange.org/gfvo#")

class sparqlTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        top_dir = generate_path("../../../")
        os.chdir(top_dir)
        src = os.path.join(os.getcwd(),"db/bigdata.jnl")
        dst = os.path.join(os.getcwd(),"db/bigdata.jnl.bk")
        subprocess.call("bash bash/kill_port_9999", shell=True)
        print "Killing existing Blazegraph process"
        subprocess.call("cp %s %s" %(src, dst), shell=True)
        subprocess.call("bash bash/start_blazegraph", shell=True)
        cls.setupBlazegraph()


    @classmethod
    def tearDownClass(cls):
        top_dir = generate_path("../../../")
        os.chdir(top_dir)
        src = os.path.join(os.getcwd(),"db/bigdata.jnl.bk")
        dst = os.path.join(os.getcwd(),"db/bigdata.jnl")
        subprocess.call("bash bash/kill_port_9999", shell=True)
        print "Killing existing Blazegraph process"
        subprocess.call("cp %s %s" %(src, dst), shell=True)
        subprocess.call("rm -f %s" % src, shell=True)
        subprocess.call("bash bash/start_blazegraph", shell=True)

    @classmethod
    def setupBlazegraph(cls):
        g = Graph()

        g.add((n.btaurus, rdf.type, owl.NamedIndividual))
        g.add((n.btaurus, rdfs.label, Literal("Bos taurus (cow)", datatype=XSD.string)))
        g.add((n.btaurus, n.is_object_of, n.from_btaurus))
        g.add((n.from_btaurus, rdf.type, n.isolation_from_host))

        g.add((n.asymptomatic, rdf.type, n.isolation_syndrome))
        g.add((n.asymptomatic, rdfs.label, Literal("Asymptomatic", datatype=XSD.string)))

        g.add((n.enteral_feeding_tube, rdf.type, n.isolation_from_source))
        g.add((n.enteral_feeding_tube, rdfs.label, Literal("Enteral feeding tube", datatype=XSD.string)))

        g.add((n.CP001165, rdf.type, gfvo.Genome))
        g.add((n.CP001165, n.has_biosample, Literal("2603441", datatype=XSD.string)))
        g.add((n.CP001165, n.has_accession, Literal("CP001165", datatype=XSD.string)))

        g.add((n.CP001164, rdf.type, gfvo.Genome))
        g.add((n.CP001164, n.has_valid_sequence, Literal("True", datatype=XSD.string)))
        g.add((n.CP001164, n.has_biosample, Literal("2603441", datatype=XSD.string)))
        g.add((n.CP001164, n.has_accession, Literal("CP001164", datatype=XSD.string)))
        g.add((n.CP001164, n.has_sequence, n.CP001164_seq))
        g.add((n.CP001164_seq, n.is_from, Literal("CORE", datatype=XSD.string)))

        g.add((n.CP001163, rdf.type, gfvo.Genome))
        g.add((n.CP001163, n.has_biosample, Literal("2603441", datatype=XSD.string)))
        g.add((n.CP001163, n.has_accession, Literal("CP001163", datatype=XSD.string)))

        g.add((n.fakeGenome, rdf.type, gfvo.Genome))
        g.add((n.fakeGenome, n.has_valid_sequence, Literal("False", datatype=XSD.string)))

        g.add((n.testEntity, rdf.type, owl.NamedIndividual))


        g.add((BNode(), rdfs.label, Literal("Test node", datatype=XSD.string)))

        g.add((n.test_object, n.has_checksum, Literal("asdfghjkl", datatype=XSD.string)))

        BlazegraphUploader().upload_data(generate_output(g))

        del g

    def test_find_from_host(self):
        self.assertEqual(_sparql.find_from_host("Bos taurus (cow)"), "from_btaurus")
        self.assertIsNone(_sparql.find_from_host("Asymptomatic"))
        self.assertIsNone(_sparql.find_from_host("Enteral feeding tube"))
        self.assertIsNone(_sparql.find_from_host("foo"))

    def test_find_syndrome(self):
        self.assertEqual(_sparql.find_syndrome("Asymptomatic"), "asymptomatic")
        self.assertIsNone(_sparql.find_syndrome("Enteral feeding tube"))
        self.assertIsNone(_sparql.find_syndrome("Bos taurus (cow)"))
        self.assertIsNone(_sparql.find_syndrome("foo"))

    def test_find_source(self):
        self.assertEqual(_sparql.find_source("Enteral feeding tube"), "enteral_feeding_tube")
        self.assertIsNone(_sparql.find_source("Asymptomatic"))
        self.assertIsNone(_sparql.find_source("Bos taurus (cow)"))
        self.assertIsNone(_sparql.find_source("foo"))

    def test_check_NamedIndividual(self):
        self.assertTrue(_sparql.check_NamedIndividual("btaurus"))
        self.assertFalse(_sparql.check_NamedIndividual("foo"))

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
        self.assertTrue(_sparql.check_NamedIndividual("testEntity"))
        _sparql.delete_instance("testEntity")
        self.assertFalse(_sparql.check_NamedIndividual("testEntity"))

    def test_insert_accession_sequence(self):
        self.assertFalse(_sparql._sparql_query(
            'PREFIX : <https://github.com/superphy#>\n'
            'ASK { :fakeGenome :has_accession ?o}'
        )["boolean"])

        self.assertFalse(_sparql._sparql_query(
            'PREFIX : <https://github.com/superphy#>\n'
            'ASK { :fakeGenome :has_sequence ?o}'
        )["boolean"])

        _sparql.insert_accession_sequence("fakeGenome", "fakePlasmid", "fakePlasmidSeq")

        self.assertTrue(_sparql._sparql_query(
            'PREFIX : <https://github.com/superphy#>\n'
            'ASK { :fakeGenome :has_accession "fakePlasmid"^^xsd:string}'
        )["boolean"])

        self.assertTrue(_sparql._sparql_query(
            'PREFIX : <https://github.com/superphy#>\n'
            'ASK { :fakeGenome :has_sequence :fakePlasmidSeq}'
        )["boolean"])


    def test_blank_nodes(self):
        self.assertTrue(_sparql.check_blank_nodes())
        _sparql.delete_blank_nodes()
        self.assertFalse(_sparql.check_blank_nodes())

    def test_check_checksum(self):
        self.assertTrue(_sparql.check_checksum("asdfghjkl"))
        self.assertFalse(_sparql.check_checksum("123456789"))

if __name__ == '__main__':
    unittest.main()

