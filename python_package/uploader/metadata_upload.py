#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
This module uploads genome metadata that has been prepared by another script and converted into a JSON format

TODO: join the two workflows someday as to avoid file IO issues and to allow for parallel operation and batching
"""

from ijson.backends import YAJLImportError

try:
    import ijson.backends.yajl2 as ijson
except YAJLImportError:
    import ijson.backends.yajl as ijson

import sys
import traceback
from collections import defaultdict
from rdflib import Graph

from _utils import generate_path, generate_output
from _eutils import return_elink_uid, return_esearch_uid
from _sparql import check_NamedIndividual

from classes import PendingGenome
from blazegraph_upload import BlazegraphUploader

reload(sys)
sys.setdefaultencoding("utf-8")

__author__ = "Stephen Kan"
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Stephen Kan"
__email__ = "stebokan@gmail.com"


class MetadataUploader(object):
    """

    """
    def __init__(self, filename, organism):
        """
        Initializes the class with zeroed counters and containers.

        Args:
            filename (str): the relative filepath from this script
            organism (str): the SPARQL URI for the organism of interest being uploaded

        """
        self.progress = 0
        self.error = 0
        self.filename = filename
        self.organism = organism

    def upload(self):
        with open(generate_path(self.filename), "r") as fd:
            data = ijson.parse(fd)
            metadata = None

            for prefix, event, value in data:
                if ("." not in prefix) and (prefix is not "") and (event == "start_map"):
                    self.add_to_graph(metadata)
                    metadata = GenomeMetadata(prefix, self.organism)

                if prefix.endswith(".displayname"):
                    metadata.add_genome_parameter(prefix.split(".",3)[1],value)

            self.add_to_graph(metadata)

        print "%d genomes parsed, %d errors occurred." % (self.progress, self.error)


    def add_to_graph(self, metadata):
        if metadata:
            try:
                self.progress += 1
                print "%d: downloading files" % self.progress
                self.create_pending_genome(metadata)
            except Exception as e:
                self.error_logging(metadata.name)

    def error_logging(self, name):
        self.error += 1
        with open(generate_path("outputs/errors.txt"), "a") as f:
            f.write("%s \n\n %s \n ================================ \n\n" % (name, traceback.format_exc()))
        print "Error %d occurred." % self.error

    def create_pending_genome(self, metadata):
        g = Graph()
        name = metadata.name

        if check_NamedIndividual(name):
            print "%s already in Blazegraph." % name
        else:
            self.get_ncbi_ids(metadata)
            kwargs = metadata.build_genome_kwargs()
            PendingGenome(g, **kwargs).rdf()
            BlazegraphUploader().upload_data(generate_output(g))


    def get_ncbi_ids(self, metadata):
        nuccore_id = return_esearch_uid("nuccore", metadata.name)

        for item in return_elink_uid("nuccore", "bioproject", nuccore_id):
            metadata.add_genome_parameter("bioproject", item)

        for item in return_elink_uid("nuccore", "biosample", nuccore_id):
            metadata.add_genome_parameter("biosample", item)


class GenomeMetadata(object):

    genome_params = {"isolation_date": "date", "isolation_location": "location", "isolation_host": "host",
                     "isolation_source": "source"}

    def __init__(self, prefix, organism):
        self.name = prefix
        self.organism = organism
        self.dict = defaultdict(set)


    def add_genome_parameter(self, key, value):
        self.dict[key].add(value)


    def build_genome_kwargs(self):
        kwargs = {'name':self.name, 'organism':self.organism}

        for key, value in self.dict.iteritems():
            if key in self.genome_params:
                param = self.genome_params.get(key)
                kwargs.update({param: value})

            elif key == "serotype":
                kwargs.update(self.get_serotypes(value))

            else:
                kwargs.update({key: value})

        return kwargs

    def get_serotypes(self, serotypes):
        Otype = None
        Htype = None

        for serotype in serotypes:
            if "ONT" in serotype:
                pass
            else:
                Otype = serotype.split(":")[0][1:]

            if "NA" in serotype:
                pass
            elif "NM" in serotype:
                Htype = "-"
            else:
                Htype = serotype.split(":")[1][1:]

        return {"Otype": Otype, "Htype": Htype}



MetadataUploader("samples/meta_pipe_result.json", "ecoli").upload()