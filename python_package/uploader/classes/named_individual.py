#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from rdflib import Literal, Namespace, XSD
from namespaces import *

class NamedIndividual(object):
    """
    The superclass for all RDF Object instances (excluding Blank Nodes as they are not named)
    """

    def __init__(self, graph, name):
        """
        Create a NamedIndividual with associated metadata

        Args:
            graph (rdflib.Graph): container object to store RDF triples
            name (str): name of the instance
        """
        self.graph = graph
        self.name = name

    def rdf(self):
        """
        Convert NamedIndividual metadata into RDF and adds them to the graph
        """

        self.graph.add((n[self.name], rdf.type, owl.NamedIndividual))