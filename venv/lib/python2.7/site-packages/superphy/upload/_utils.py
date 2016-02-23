#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
A module containing some generic utility functions for the project

"""

import inspect
import os
import string

__author__ = "Stephen Kan"
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Stephen Kan"
__email__ = "stebokan@gmail.com"

def generate_path(filename):
    """
    Generates the absolute filepath based on the location of the caller of this function

    Args:
        filename (str): relative location of the caller

    Returns: absolute filepath for the given filename based on the location of the caller

    """
    frame = inspect.stack()
    (frame, filepath, line_number, function_name, lines, index) = frame[1]
    del frame
    return os.path.join(os.path.dirname(filepath), filename)

def strip_non_alphabetic(str):
    """
    Strips all non-alphabetic characters (not present in string.ascii_letters) from a given string

    Args:
        str: any string

    Returns: a string with only characters from string.ascii_letters

    """
    all = string.maketrans('','')
    nochars = all.translate(all, string.ascii_letters)
    return str.translate(all, nochars)

def strip_non_numeric(str):
    """
    Strips all non-numeric characters (not present in string.digits) from a given string

    Args:
        str: any string

    Returns: a string with only characters from string.digits

    """
    all = string.maketrans('','')
    nodigs = all.translate(all, string.digits)
    return str.translate(all, nodigs)

def generate_output(graph):
    """
    Returns RDF Graph data in the turtle format and clears the Graph

    Args:
        graph (rdflib.Graph): container object to store RDF triples
    """

    output = graph.serialize(format="turtle")
    graph.remove( (None, None, None) )
    return output

def generate_file_output(graph, destination):
    """
    Export RDF Graph data to a turtle file at the given destination

    Args:
        graph (rdflib.Graph): container object to store RDF triples
        destination (str): an internal filepath relative to the  __init__.py file this module belongs to
    """

    graph.serialize(destination=destination, format="turtle")
    graph.remove( (None, None, None) )