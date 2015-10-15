__author__ = 'Stephen Kan'

from superphy_classes import Host, HostCategory, FromSource, IsolationSyndrome, Microbe, Htype, Otype, generate_output
import json

"""
This module converts predefined metadata for Hosts, HostCategory, FromSource, IsolationSyndrome, and Microbe (all
from superphy_classes.py) into RDF triples and exports it as a turtle file ready for converting into Blazegraph.

TODO: add a way to modify source files by reading the json, check for duplicates, and write a new JSON
      Would be nice to do this serially
"""


def convert_host_categories():
    file = open("data/host_categories.txt", "r+")
    host_categories = json.load(file)

    for host_category in host_categories:
        name, label = host_category
        HostCategory(name, label).rdf()


def convert_hosts():
    file = open("data/hosts.txt", "r+")
    hosts = json.load(file)

    for host in hosts:
        name, label, sci_name, com_name, host_category = host
        Host(name, host_category, label, sci_name, com_name).rdf()


def convert_sources():
    file = open("data/sources.txt", "r+")
    sources = json.load(file)

    for source in sources:
        name, label, host_category = source
        FromSource(name, label, host_category).rdf()


def convert_syndromes():
    file = open("data/syndromes.txt", "r+")
    syndromes = json.load(file)

    for syndrome in syndromes:
        name, label, host_category = syndrome
        IsolationSyndrome(name, label, host_category).rdf()


def convert_microbes():
    file = open("data/microbes.txt", "r+")
    microbes = json.load(file)

    for microbe in microbes:
        name, label, sci_name, com_name = microbe
        Microbe(name, label, sci_name, com_name).rdf()


def generate_serotypes():
    for num in range(1,56):
        Htype(num).rdf()
    Htype("Unknown").rdf()
    Htype("-").rdf()

    for num in range(1,187):
        Otype(num).rdf()
    Otype("Unknown").rdf()


def generate_all():
    convert_host_categories()
    convert_hosts()
    convert_microbes()
    convert_sources()
    convert_syndromes()
    generate_serotypes()
    generate_output("outputs/setup.ttl")
