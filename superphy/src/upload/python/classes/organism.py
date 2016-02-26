#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Organism class converts inputted data into RDF triples in accordance with the Superphy ontology

"""

from .named_individual import *
from attribute import FromHost

class Organism(NamedIndividual):
    """
    Umbrella class for a biological organism with associated metadata.
    """

    def __init__(self, graph, name, label, scientific_name, common_name, taxonomy_id=None):
        """
        Create an Organism with associated metadata

        Args:
            graph (rdflib.Graph): container object to store RDF triples
            name (str): the name of the Organism
            label (str): a label used by Meta::Miner to refer to the Organism
            scientific_name (str): the scientific name (genus species) of the Organism
            common_name (str): the common name of the Organism
            taxonomy_id: the taxonomy id for the Organism as given by NCBI
        """

        super(Organism, self).__init__(graph, name)
        self.label = label
        self.scientific_name = scientific_name
        self.common_name = common_name
        self.taxonomy_id = taxonomy_id

    def rdf(self):
        """
        Convert Organism metadata into RDF and adds them to the graph
        """

        super(Organism, self).rdf()
        self.graph.add((n[self.name], rdf.type, n.Organism))
        self.graph.add((n[self.name], rdfs.label, Literal(str(self.label), datatype=XSD.string)))
        self.graph.add((n[self.name], n.scientific_name, Literal(str(self.scientific_name), datatype=XSD.string)))
        self.graph.add((n[self.name], n.common_name, Literal(str(self.common_name), datatype=XSD.string)))
        self.graph.add((n[self.name], n.has_taxonomy_id, Literal(str(self.taxonomy_id), datatype=XSD.string)))


class Host(Organism):
    """
    A host organism for another organism with associated metadata
    """

    def __init__(self, graph, name, label, scientific_name, common_name, host_category, taxonomy_id=None):
        """
        Create a Host instance with associated metadata

        Args:
            - refer to Organism
            host_category (str): the host category that Host belongs to
        """

        super(Host, self).__init__(graph, name, label, scientific_name, common_name, taxonomy_id)
        self.host_category = host_category

    def rdf(self):
        """
        Convert Host metadata into RDF and adds them to the graph
        """

        super(Host, self).rdf()
        self.graph.add((n[self.name], rdf.type, n.Host))
        self.graph.add((n[self.name], n.has_host_category, n[self.host_category]))
        self.graph.add((n[self.host_category], n.is_host_category_of, n[self.name]))
        FromHost(self.graph, self.name, self.host_category).rdf()


class Microbe(Organism):
    """
    A microbe with associated metadata
    """

    def rdf(self):
        """
        Convert Microbe metadata into RDF and adds them to the graph
        """

        super(Microbe, self).rdf()
        self.graph.add((n[self.name], rdf.type, n.Microbe))