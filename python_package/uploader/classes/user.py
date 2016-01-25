#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from named_individual import *
from namespaces import *

class User(NamedIndividual):
    """
    This class is created when a new user is registered.
    """

    def __init__(self, graph, email):
        """
        Create a new User with associated metadata

        Args:
            graph (rdflib.Graph): container object to store RDF triples
            email (str): is both the individual name, and the only field literal
        """

        super(User, self).__init__(graph, email)
        self.email = email

    def rdf(self):
        """
        Convert User metadata into RDF and adds them to the graph
        """

        super(User, self).rdf()
        self.graph.add((n[self.name], rdf.type, n.User))
        self.graph.add((n[self.name], n.email, Literal(str(self.email), datatype=XSD.string)))