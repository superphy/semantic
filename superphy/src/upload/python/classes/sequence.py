#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from named_individual import *

"""
Sequence: a DNA sequence belonging to a genome

"""

class Sequence(NamedIndividual):
    """
    A sequence and its associated informatics metadata.
    """

    def __init__(self, graph, name, genome, sequences, bp, contigs, checksum, is_from):
        """
        Create a Sequence with its associated metadata.
        Args:
            graph(rdflib.Graph): container object to store RDF triples
            name (str): name of the sequence (Genome name with "_seq" appended to it
            genome (str): the genome that the sequence belongs to
            sequence (list[str]): a list of sequences (contigs) that the sequence record has
            bp (int): number of base pairs in the sequence
            contigs (int): number of contiguous sequences in the sequence
            checksum (str): the hash generated from a hashlib used to confirm sequences are unique
            is_from (str): identifier explaining origin of sequence (e.g. WGS, PLASMID, CORE)
        """
        super(Sequence, self).__init__(graph, name)
        self.genome = genome
        self.sequences = sequences
        self.bp = bp
        self.contigs = contigs
        self.checksum = checksum
        self.is_from = is_from


    def rdf(self):
        """
        Convert all Sequence variables to RDF and adds it to the graph
        """
        super(Sequence, self).rdf()
        self.graph.add((n[self.name], rdf.type, n.Sequence))

        for sequence in self.sequences:
            self.graph.add((n[self.name], n.has_value, Literal(str(sequence), datatype=XSD.string)))

        self.graph.add((n[self.name], n.has_base_pair, Literal(str(self.bp), datatype=XSD.int)))
        self.graph.add((n[self.name], n.has_contigs, Literal(str(self.contigs), datatype=XSD.int)))
        self.graph.add((n[self.name], n.has_checksum, Literal(str(self.checksum), datatype=XSD.string)))
        self.graph.add((n[self.name], n.is_from, Literal(str(self.is_from), datatype=XSD.string)))
        self.graph.add((n[self.genome], n.has_sequence, n[self.name]))
        self.graph.add((n[self.name], n.is_sequence_of, n[self.genome]))


    def add_seq_validation(self, boolean):
        """
        Converts result of sequence validation to RDF and adds it to the graph
        Args:
            boolean (bool): indicates whether or not the sequence has passed validation
        """
        self.graph.add((n[self.genome], n.has_valid_sequence, Literal(str(boolean), datatype=XSD.string)))

    def add_hits(self, hits):
        """
        Converts a list of hits from sequence validation to RDF and adds it to the graph
        Args:
            hits(list[str]): list of validating regions that the sequence had 90%+ identity with
        """
        for hit in hits:
            self.graph.add((n[self.name], n.has_hit, Literal(str(hit), datatype=XSD.string)))