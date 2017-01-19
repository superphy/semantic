#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
A module containing some generic utility functions for the project

"""

import inspect
import os
import string

def generate_path(filename):
    """
    Generates the absolute filepath based on the location of the caller of this
    function

    Args:
        filename (str): relative location of the caller

    Returns: absolute filepath for the given filename based on the location of
    the caller

    """
    frame = inspect.stack()
    (frame, filepath, line_number, function_name, lines, index) = frame[1]
    del frame
    return os.path.join(os.path.dirname(filepath), filename)

def strip_non_alphabetic(str_):
    """
    Strips all non-alphabetic characters (not present in string.ascii_letters)
    from a given string

    Args:
        str: any string

    Returns: a string with only characters from string.ascii_letters

    """
    all_ = string.maketrans('', '')
    nochars = all_.translate(all_, string.ascii_letters)
    return str_.translate(all_, nochars)

def strip_non_numeric(str_):
    """
    Strips all non-numeric characters (not present in string.digits) from a
    given string

    Args:
        str: any string

    Returns: a string with only characters from string.digits

    """
    all_ = string.maketrans('', '')
    nodigs = all_.translate(all_, string.digits)
    return str_.translate(all_, nodigs)

def generate_output(graph):
    """
    Returns RDF Graph data in the turtle format and clears the Graph

    Args:
        graph (rdflib.Graph): container object to store RDF triples
    """

    output = graph.serialize(format="turtle")
    return output

def generate_file_output(graph, destination):
    """
    Export RDF Graph data to a turtle file at the given destination

    Args:
        graph (rdflib.Graph): container object to store RDF triples
        destination (str): an internal filepath relative to the  __init__.py
        file this module belongs to
    """

    graph.serialize(destination=destination, format="turtle")

def parse_nih_name(description):
    """
    Parses a String of a nih name (eg. record.description after Bio.SeqIO.parse)
    and returns a dictionary of the substrings we're interesting FOR creating
    uriSpecies

    Args:
        description (str): a record.description
        ex. gi|427200135|gb|ANLJ01000508.1| Escherichia coli 89.0511 gec890511.contig.603_1, whole genome shotgun sequence
    Returns:
        (dict) with keys: accession_id, species, assembly, contig
        ex. {'accession_id': 'ANLJ01000508', 'contig': '000508', 'assembly': 'ANLJ01', 'species': '89.0511'}

    TODO:
        -add code to parse other nih naming conventions
        -what happens when no species name??
    """

    identifiers = {'accession_id' : description.split("|")[3].split(".")[0]}
    identifiers['species'] = description.split("|")[4].split(" ")[3]
    identifiers['assembly'] = identifiers['accession_id'][0:6]
    identifiers['contig'] = identifiers['accession_id'][6:12]
    return identifiers

def generate_uri(uri, s=''):
    """
    Takes a string as one would define for .ttl files and returns a URI for rdflib.

    Args:
        uri (str): a string following .ttl convention for a URI
        ex. g:Identifier as shorthand for http://www.biointerchange.org/gfvo#Identifier
    Returns:
        (rdflib.URIRef) with URI needed to add to rdflib.Graph
    """

    from rdflib import Namespace, URIRef, Literal
    from ConfigParser import SafeConfigParser

    #if you call with a uri already
    if isinstance(uri, URIRef):
        return URIRef(str(uri) + s)

    parser = SafeConfigParser()
    parser.read('defaults.cfg')
    prefix = uri.split(':')[0]
    postfix = uri.split(':')[1]
    if prefix == '': #this is our : case
        return URIRef(parser.get('Namespaces', 'root') + postfix)
    else:
        return URIRef(parser.get('Namespaces', prefix) + postfix)

def from_nuccore(accession):
    """Obtains the FASTA sequence via the NCBI Genbank Nucleotide database
    using Entrez EUtils. If there
    is nothing found for the sequence, raise a ValueError.

    Args:
        accession (str): genbank accession id
    """

    from Bio import Entrez, SeqIO

    Entrez.email = "superphy.info@gmail.com"
    handle = None
    i = 0

    while i < 3:
        try:
            handle = Entrez.efetch(
                db="nuccore",
                id=accession,
                rettype="fasta",
                retmode="text"
            )
            SeqIO.write(handle, 'tmp/' + accession + '.fasta', 'fasta')
            return 'tmp/' + accession + '.fasta'
        except:
            i += 1
            continue
    try:
        handle is None
    except NameError:
        raise TypeError("Could not retrieve file for analysis")
