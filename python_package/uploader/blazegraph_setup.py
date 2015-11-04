__author__ = 'Stephen Kan'

from classes import Host, HostCategory, FromSource, IsolationSyndrome, Microbe, Htype, Otype, generate_file_output
from rdflib import Graph
import json
import os
import inspect

"""
This module converts predefined metadata for Hosts, HostCategory, FromSource, IsolationSyndrome, and Microbe (all
from classes.py) into RDF triples and exports it as a turtle file ready for converting into Blazegraph.

TODO: add a way to modify source files by reading the json, check for duplicates, and write a new JSON
      Would be nice to do this serially
"""
g = Graph()

def convert_host_categories():
    host_categories = generate_json("data/host_categories.txt")

    for host_category in host_categories:
        name, label = host_category
        HostCategory(g, name, label).rdf()


def convert_hosts():
    hosts = generate_json("data/hosts.txt")

    for host in hosts:
        name, label, sci_name, com_name, host_category = host
        Host(g, name, label, sci_name, com_name, host_category).rdf()


def convert_sources():
    sources = generate_json("data/sources.txt")

    for source in sources:
        name, label, host_category = source
        FromSource(g, name, label, host_category).rdf()


def convert_syndromes():
    syndromes = generate_json("data/syndromes.txt")

    for syndrome in syndromes:
        name, label, host_category = syndrome
        IsolationSyndrome(g, name, label, host_category).rdf()


def convert_microbes():
    microbes = generate_json("data/microbes.txt")

    for microbe in microbes:
        name, label, sci_name, com_name = microbe
        Microbe(g, name, label, sci_name, com_name).rdf()


def generate_serotypes():
    for num in range(1,56):
        Htype(g, num).rdf()
    Htype(g, "Unknown").rdf()
    Htype(g, "-").rdf()

    for num in range(1,187):
        Otype(g, num).rdf()
    Otype(g, "Unknown").rdf()


def generate_json(filename):
    path = os.path.join(os.path.dirname(inspect.getfile(inspect.currentframe())), filename)
    file = open(path, "r+")
    return json.load(file)


def generate_all():
    convert_host_categories()
    convert_hosts()
    convert_microbes()
    convert_sources()
    convert_syndromes()
    generate_serotypes()
    generate_file_output(g, os.path.join(os.path.dirname(inspect.getfile(inspect.currentframe())), 'ontologies/setup.ttl'))
