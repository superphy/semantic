#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Classes for a Genome:
    Genome: a genome sequenced from a sample
    PendingGenome: a genome that has not finished processing
    CompletedGenome: a genome that has finished processing
"""

from named_individual import *
from Superphy.models.upload import _sparql

class Genome(NamedIndividual):
    """
    A genome with associated metadata

    """

    def __init__(self, graph, name, **kwargs):
        """
        Create a Genome with associated metadata

        Since most of the arguments for genome are optional, **kwargs is used to pass them in, after
        filtering them via SearchParam's keywords to prevent spurious and potentially malicious keys
        from being kept.

        Attributes:
            searchparam: list of keywords to filter kwargs by

        Args:
            graph (rdflib.Graph): container object to store RDF triples
            name (str): name of the Genome
            **kwargs (dict): optional arguments for Genome, that will be filtered for
                           spurious entries

        """

        searchparam = ["isolation_date", "isolation_location", "accession", "bioproject", "biosample", "strain",
                       "organism", "isolation_host", "isolation_source", "syndrome", "Htype", "Otype", "User"]

        super(Genome, self).__init__(graph, name)
        self.kwargs = {key: value for key, value in kwargs.items() if key in searchparam}

    def rdf(self):
        """
        Convert all Genome metadata to RDF and adds it to the graph

        To ensure that H and O types are assigned Unknowns, if kwargs does not include at least
        one of them, it would call the appropriate RDF method.

        As methods are called based on their keys, there is some coupling with minerJSON.py to
        ensure that the keys are properly named.
        """

        super(Genome, self).rdf()

        if "Htype" not in self.kwargs:
            self.Htype()
        if "Otype" not in self.kwargs:
            self.Otype()

        for key, value in self.kwargs.iteritems():
            getattr(self, key)(value)

        self.graph.add((n[self.name], rdf.type, gfvo.Genome))

    def isolation_date(self, date):
        """
        Convert all date entries into RDF and adds them to the graph

        Args:
            date: a collection of sampling dates for the Genome in the XSD date format (YYYY-MM-DD)
        """

        for item in date:
            literal = Literal(item, datatype=XSD.date)
            self.graph.add((n[self.name], n.has_isolation_date, literal))

    def isolation_location(self, location):
        """
        Convert all location entries into RDF and adds them to the graph

        Args:
            location: a collection of sampling locations for the Genome
        """

        for item in location:
            literal = Literal(item, datatype=XSD.string)
            self.graph.add((n[self.name], n.has_geographic_location, literal))

    def accession(self, accession):
        """
        Convert all NCBI Genbank/Nucleotide accession numbers into RDF and adds then to the graph

        Args:
            accession: a collection of Nucleotide accession ids associated with the Genome
        """

        for item in accession:
            literal = Literal(item, datatype=XSD.string)
            self.graph.add((n[self.name], n.has_accession, literal))

    def bioproject(self, bioproject):
        """
        Convert all BioProject ids into RDF and adds them to the graph

        Args:
            bioproject: a collection of BioProject ids associated with the Genome
        """

        for item in bioproject:
            literal = Literal(item, datatype=XSD.string)
            self.graph.add((n[self.name], n.has_bioproject, literal))

    def biosample(self, biosample):
        """
        Convert all BioSample ids into RDF and adds them to the graph

        Args:
            biosample: a collection of BioSample ids associated with the Genome
        """

        for item in biosample:
            literal = Literal(item, datatype=XSD.string)
            self.graph.add((n[self.name], n.has_biosample, literal))

    def strain(self, strain):
        """
        Convert all strain names into RDF and adds them to the graph

        Args:
            strain: a collection of strain names associated with the Genome
        """

        for item in strain:
            literal = Literal(item, datatype=XSD.string)
            self.graph.add((n[self.name], n.has_strain, literal))

    def organism(self, organism):
        """
        Convert organism into RDF and adds it to the graph

        Args:
            organism (str): name of the organism of the Genome
        """

        self.graph.add((n[self.name], n.is_genome_of, n[organism]))
        self.graph.add((n[organism], n.has_genome, n[self.name]))

    def isolation_host(self, from_host):
        """
        Convert all host data into RDF and adds them to the graph

        Args:
            from_host: a collection of hosts that the Genome has been sampled from
        """

        for item in from_host:
            node = _sparql.find_from_host(item)
            self.graph.add((n[self.name], n.has_isolation_attribute, n[node]))
            self.graph.add((n[node], n.is_isolation_attribute_of, n[self.name]))

    def isolation_source(self, from_source):
        """
        Convert all source data into RDF and adds them to the graph

        Args:
            from_source: a collection of biological sources that the Genome has been sampled from
        """

        for item in from_source:
            node = _sparql.find_source(item)
            self.graph.add((n[self.name], n.has_isolation_attribute, n[node]))
            self.graph.add((n[node], n.is_isolation_attribute_of, n[self.name]))

    def syndrome(self, syndrome):
        """
        Convert all syndrome data into RDF and adds them to the graph

        Args:
            syndrome: a collection of syndromes associated with the Genome
        """

        for item in syndrome:
            node = _sparql.find_syndrome(item)
            self.graph.add((n[self.name], n.has_isolation_attribute, n[node]))
            self.graph.add((n[node], n.is_isolation_attribute_of, n[self.name]))

    def Htype(self, Htype=None):
        """
        Convert H serotype into RDF and adds it to the graph

        Args:
            Htype: the id of the H-antigen associated with the Genome, or None if not provided
        """

        if Htype:
            self.graph.add((n[self.name], n.has_Htype, n["H" + str(Htype)]))
            self.graph.add((n["H" + str(Htype)], n.is_Htype_of, n[self.name]))

    def Otype(self, Otype=None):
        """
        Convert O serotype into RDF and adds it to the graph

        Args:
            Otype: the id of the O-antigen associated with the Genome, or None if not provided
        """

        if Otype:
            self.graph.add((n[self.name], n.has_Otype, n["O" + str(Otype)]))
            self.graph.add((n["O" + str(Otype)], n.is_Otype_of, n[self.name]))

    def User(self, User):
        """
        Converts User id into RDF and adds it to the graph

        Args:
            User: the id of the user who uploaded the Genome, restricting permissions to that user
        """

        self.graph.add((n[self.name], n.is_owned_by, n[User]))
        self.graph.add((n[User], n.owns, n[self.name]))


class PendingGenome(Genome):
    """
    A genome that has not completed sequence analysis.

    Attributes:
        Same as Genome
    """

    def rdf(self):
        """
        Convert all PendingGenome variables to RDF, tags the Genome as a PendingGenome, and adds them to the graph
        """

        super(PendingGenome, self).rdf()
        self.graph.add((n[self.name], rdf.type, n.pending_genome))


class CompletedGenome(Genome):
    """
    A genome with completed sequence analysis

    Attributes:
        Same as Genome
    """

    def rdf(self):
        """
        Convert all CompletedGenome variables to RDF, tags the Genome as a CompletedGenome, and adds them to the graph

        """

        super(CompletedGenome, self).rdf()
        self.graph.add((n[self.name], rdf.type, n.completed_genome))