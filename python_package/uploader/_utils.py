__author__ = 'Stephen Kan'

import os
import inspect
import string

def generate_path(filename):
    frame = inspect.stack()
    (frame, filepath, line_number, function_name, lines, index) = frame[1]
    del frame
    return os.path.join(os.path.dirname(filepath), filename)

def strip_non_alphabetic(str):
    all = string.maketrans('','')
    nochars = all.translate(all, string.ascii_letters)
    return str.translate(all, nochars)

def strip_non_numeric(str):
    all = string.maketrans('','')
    nodigs = all.translate(all, string.digits)
    return str.translate(all, nodigs)

def generate_output(graph):
    """
    Returns RDF Graph data in the turtle format and clears the Graph
    """

    output = graph.serialize(format="turtle")
    graph.remove( (None, None, None) )
    return output

def generate_file_output(graph, destination):
    """
    Export RDF Graph data to a turtle file at the given destination

    Args:
        destination: an internal filepath relative to the  __init__.py file this module belongs to
    """

    graph.serialize(destination=destination, format="turtle")
    graph.remove( (None, None, None) )