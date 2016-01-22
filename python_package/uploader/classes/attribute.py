#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Classes:
    Attribute: attributes of objects in the ontology (abstract class)
    HostCategory: host category for an object
    IsolationAttribute: details associated with isolation of a genome (abstract class)
    FromHost: the host from which the organism was extracted from
    FromSource: the biological source of the sample
    IsolationSyndrome: the syndromes associated with the sample
    Serotype: the serotype of the genome (abstract class)
    Otype: the id of the O-antigen associated with the sample
    Htype: the id of the H-antigen associated with the sample
"""

from named_individual import *
from namespaces import *

class Attribute(NamedIndividual):
    """
    An attribute class used to describe objects in the ontology
    """

    def rdf(self):
        """
        Convert Attribute metadata into RDF and adds them to the graph
        """

        super(Attribute, self).rdf()
        self.graph.add((n[self.name], rdf.type, n.Attribute))


class HostCategory(Attribute):
    """
    An host category attribute with associated metadata

    An attribute describing valid host categories for other attributes
    """

    def __init__(self, graph, name, label):
        """
        Create a HostCategory with associated metadata

        Args:
            graph (rdflib.Graph): container object to store RDF triples
            name (str): name of the HostCategory
            label (str): a label used by Meta::Miner to refer to the HostCategory
        """

        super(HostCategory, self).__init__(graph, name)
        self.label = label

    def rdf(self):
        """
        Convert HostCategory metadata into RDF and adds them to the graph
        """

        super(Attribute, self).rdf()
        self.graph.add((n[self.name], rdf.type, n.host_category))
        self.graph.add((n[self.name], rdfs.label, Literal(str(self.label), datatype=XSD.string)))


class IsolationAttribute(Attribute):
    """
    An isolation attribute with associated metadata

    Attributes that are associated or describes how a genome sample was isolated
    """

    def rdf(self):
        """
        Convert IsolationAttribute metadata into RDF and adds them to the graph
        """

        super(IsolationAttribute, self).rdf()
        self.graph.add((n[self.name], rdf.type, n.isolation_attribute))


class FromHost(IsolationAttribute):
    """
    A "from host" attribute with associated metadata
    """

    def __init__(self, graph, host, host_category):
        """
        Create a FromHost with associated metadata

        Args:
            graph (rdflib.Graph): container object to store RDF triples
            host (str): the host that FromHost references
            host_category (str): the host category that FromHost belongs to
        """

        self.name = "from_" + host
        super(FromHost, self).__init__(graph, self.name)
        self.host = host
        self.host_category = host_category

    def rdf(self):
        """
        Convert FromHost metadata into RDF and adds them to the graph
        """

        super(FromHost, self).rdf()
        self.graph.add((n[self.name], rdf.type, n.isolation_from_host))
        self.graph.add((n[self.name], n.has_object, n[self.host]))
        self.graph.add((n[self.host], n.is_object_of, n[self.name]))
        self.graph.add((n[self.name], n.has_host_category, n[self.host_category]))
        self.graph.add((n[self.host_category], n.is_host_category_of, n[self.name]))


class FromSource(IsolationAttribute):
    """
    A biological source with associated metadata
    """

    def __init__(self, graph, name, label, host_category):
        """
        Create a FromSource with associated metadata

        Args:
            graph (rdflib.Graph): container object to store RDF triples
            name (str): the name of the FromSource
            label (str): a label used by Meta::Miner to refer to the FromSource
            host_category (str): the host category of the FromSource
        """

        super(FromSource, self).__init__(graph, name)
        self.label = label
        self.host_category = host_category

    def rdf(self):
        """
        Convert FromSource metadata into RDF and adds them to the graph
        """

        super(FromSource, self).rdf()
        self.graph.add((n[self.name], rdf.type, n.isolation_from_source))
        self.graph.add((n[self.name], rdfs.label, Literal(str(self.label), datatype=XSD.string)))
        self.graph.add((n[self.name], n.has_host_category, n[self.host_category]))
        self.graph.add((n[self.host_category], n.is_host_category_of, n[self.name]))


class IsolationSyndrome(IsolationAttribute):
    """
    A syndrome with associated metadata
    """

    def __init__(self, graph, name, label, host_category):
        """
        Create an IsolationSyndrome with associated metadata

        Args:
            graph (rdflib.Graph): container object to store RDF triples
            name (str): the name of the IsolationSyndrome
            label (str): a label used by Meta::Miner to refer to the IsolationSyndrome
            host_category (str): the host category of the IsolationSyndrome
        """
        super(IsolationSyndrome, self).__init__(graph, name)
        self.label = label
        self.host_category = host_category

    def rdf(self):
        """
        Convert IsolationSyndrome metadata into RDF and adds them to the graph
        """

        super(IsolationSyndrome, self).rdf()
        self.graph.add((n[self.name], rdf.type, n.isolation_syndrome))
        self.graph.add((n[self.name], rdfs.label, Literal(str(self.label), datatype=XSD.string)))
        self.graph.add((n[self.name], n.has_host_category, n[self.host_category]))
        self.graph.add((n[self.host_category], n.is_host_category_of, n[self.name]))


class Serotype(Attribute):
    """
    A serotype with associated metadata
    """

    def rdf(self):
        """
        Convert Serotype metadata into RDF and adds them to the graph
        """

        super(Serotype, self).rdf()
        self.graph.add((n[self.name], rdf.type, n.serotype))


class Otype(Serotype):
    """
    A O-type antigen with associated metadata
    """

    def __init__(self, graph, id):
        """
        Create a Otype with associated metadata

        Args:
            graph (rdflib.Graph): container object to store RDF triples
            id (str): the id of the O antigen
        """
        self.id = str(id)
        self.name = "O" + str(id)
        super(Otype, self).__init__(graph, self.name)

    def rdf(self):
        """
        Convert Otype metadata into RDF and adds them to the graph
        """
        super(Otype, self).rdf()
        self.graph.add((n[self.name], rdf.type, n.Otype))
        literal = Literal(self.id, datatype=XSD.string)
        self.graph.add((n[self.name], rdfs.label, literal))


class Htype(Serotype):
    """
    A H-type antigen with associated metadata
    """

    def __init__(self, graph, id):
        """
        Create a Htype with associated metadata

        Args:
            graph (rdflib.Graph): container object to store RDF triples
            id (str): the id of the H antigen
        """
        self.id = str(id)
        self.name = "H" + str(id)
        super(Htype, self).__init__(graph, self.name)

    def rdf(self):
        """
        Convert Htype metadata into RDF and adds them to the graph
        """

        super(Htype, self).rdf()
        self.graph.add((n[self.name], rdf.type, n.Htype))
        literal = Literal(self.id, datatype=XSD.string)
        self.graph.add((n[self.name], rdfs.label, literal))