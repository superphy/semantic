#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
This module converts predefined metadata for Hosts, HostCategory, FromSource,
IsolationSyndrome, and Microbe (all from classes.py) into RDF triples and
exports it as a turtle file ready for converting into Blazegraph.

TODO: add a way to modify source files by reading the json,
    check for duplicates, and write a new JSON.
    Would be nice to do this serially
"""
import json
import os

from rdflib import Graph

from superphy.upload._utils import generate_file_output, generate_path
from superphy.upload.classes import FromSource, Host, HostCategory, Htype, \
    IsolationSyndrome, Microbe, Otype


__author__ = "Stephen Kan"
__copyright__ = """
    © Copyright Government of Canada 2012-2015. Funded by the Government of
    Canada Genomics Research and Development Initiative
    """
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Stephen Kan"
__email__ = "stebokan@gmail.com"


class BlazegraphSetup(object):
    """
    A class that sets up JSON files containing curated data and exports them
    into an turtle file.
    """

    def __init__(self):
        """
        Sets up an instance of the class with a rdflib.Graph() instance

        """
        self.graph = Graph()

    def convert_host_categories(self):
        """
        Loads host category data from JSON into the rdflib.Graph()

        """
        host_categories = self.import_json("data/host_categories.txt")

        for host_category in host_categories:
            name, label = host_category
            HostCategory(self.graph, name, label).rdf()

        del host_categories


    def convert_hosts(self):
        """
        Loads host data from JSON into the rdflib.Graph()

        """
        hosts = self.import_json("data/hosts.txt")

        for host in hosts:
            name, label, sci_name, com_name, host_category = host
            Host(
                self.graph, name, label, sci_name, com_name, host_category
            ).rdf()

        del hosts


    def convert_sources(self):
        """
        Loads isolation source data from JSON into the rdflib.Graph()

        """
        sources = self.import_json("data/sources.txt")

        for source in sources:
            name, label, host_category = source
            FromSource(self.graph, name, label, host_category).rdf()

        del sources


    def convert_syndromes(self):
        """
        Loads isolation syndrome data from JSON into the rdflib.Graph()

        """
        syndromes = self.import_json("data/syndromes.txt")

        for syndrome in syndromes:
            name, label, host_category = syndrome
            IsolationSyndrome(self.graph, name, label, host_category).rdf()

        del syndromes


    def convert_microbes(self):
        """
        Loads microbe data from JSON into the rdflib.Graph()

        """
        microbes = self.import_json("data/microbes.txt")

        for microbe in microbes:
            name, label, sci_name, com_name = microbe
            Microbe(self.graph, name, label, sci_name, com_name).rdf()

        del microbes

    def generate_serotypes(self):
        """
        Generates a range of H and O serotypes to add to the rdflib.Graph

        """
        for num in range(1, 56):
            Htype(self.graph, str(num)).rdf()
        Htype(self.graph, "NM").rdf()

        for num in range(1, 188):
            Otype(self.graph, str(num)).rdf()
        Otype(self.graph, "NT").rdf()
    @classmethod
    def import_json(cls, filename):
        """
        Imports JSON data from the specified file into Python

        Args:
            filename (str): the relative filepath to this python function

        Returns: a Python object composed of the data from the JSON data

        """
        path = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            filename
        )
        with open(generate_path(path), "r+") as file_:
            return json.load(file_)

    def setup_curated_data(self):
        """
        Converts all curated data stored in JSON format into a turtle file
        ready for uploading into Blazegraph
        """
        self.convert_host_categories()
        self.convert_hosts()
        self.convert_microbes()
        self.convert_sources()
        self.convert_syndromes()
        self.generate_serotypes()
        generate_file_output(self.graph, generate_path('ontologies/setup.ttl'))
