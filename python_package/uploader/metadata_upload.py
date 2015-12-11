#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""This module uploads genome metadata that has been prepared by another script and converted into a JSON format.
ijson is used to parse the JSON file because it reads from file line by line, instead of storing to memory first.

TODO: join the two workflows someday as to avoid file IO issues and to allow for parallel operation and batching

Classes:
    MetadataUploader: uploads metadata from a stored JSON object
    GenomeMetadata: stores genome metadata for uploading
"""

from collections import defaultdict
import sys
import traceback

from ijson.backends import YAJLImportError
try:
    import ijson.backends.yajl2 as ijson
except YAJLImportError:
    import ijson.backends.yajl as ijson
from rdflib import Graph


from _eutils import return_elink_uid, return_esearch_uid
from _sparql import check_NamedIndividual
from _utils import generate_output, generate_path
from classes import PendingGenome
from blazegraph_upload import BlazegraphUploader

# UTF-8 is necessary to handle the encoding for non-ASCII characters in user-inputed strings (user-inputted fields are
# a big offender of this (e.g. location names)
reload(sys)
sys.setdefaultencoding("utf-8")

__author__ = "Stephen Kan"
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Stephen Kan"
__email__ = "stebokan@gmail.com"


class MetadataUploader(object):
    """A class for parsing and uploading genomic metadata from a specially prepared JSON file.

    TODO: link Genbank data extraction directly to MetadataUploader, so that we can avoid File IO issues.
    We may still need files if memory is insufficient; batch up metadata into files and upload files sequentially?
    """
    def __init__(self, filename, organism):
        """Initializes the class with zeroed counters and containers.

        Args:
            filename (str): the relative filepath from this script
            organism (str): the SPARQL URI for the organism of interest being uploaded
        """
        self.progress = 0
        self.error = 0
        self.filename = filename
        self.organism = organism

    def upload(self):
        """Uploads the contents of the given file by parsing it as an ijson stream.

        Prints out ending message regarding number of genomes processed and errors encountered
        """
        with open(generate_path(self.filename), "r") as fd:
            data = ijson.parse(fd)
            self.parse_metadata(data)

        print "%d genomes parsed, %d errors occurred." % (self.progress, self.error)

    def parse_metadata(self, data):
        """Takes in data from an ijson stream and parses for genomic metadata based on attached tags

        Args:
            data: an ijson stream created from an open file
        """
        metadata = None
        for prefix, event, value in data:
            if ("." not in prefix) and (prefix is not "") and (event == "start_map"):
                self.add_to_graph(metadata)
                metadata = GenomeMetadata(prefix, self.organism)

            if prefix.endswith(".displayname"):
                metadata.add_genome_parameter(prefix.split(".", 3)[1], value)
        self.add_to_graph(metadata)

    def add_to_graph(self, metadata):
        """Attempts to upload data to Blazegraph via conversion to RDF and turtle. Tracks errors made during this process
        as it likely has to do with either missing curated data, a blocked or inaccessible NCBI site, and issues with
        formatting and information availability.

        Args:
            metadata (GenomeMetadata): a GenomeMetadata object used to store relevant key-value pairs from the parse
        """
        if metadata:
            try:
                self.progress += 1
                print "%d: downloading files" % self.progress
                self.create_pending_genome(metadata)
            except Exception as e:
                self.error_logging(metadata.name)

    def error_logging(self, name):
        """Records the trackback of any error messages to an log file so that if any are encountered, the log file will
        retain pertinent information for debugging

        Args:
            name(str): The genome that is currently being uploaded
        """
        self.error += 1
        with open(generate_path("outputs/errors.txt"), "a") as f:
            f.write("%s \n\n %s \n ================================ \n\n" % (name, traceback.format_exc()))
        print "Error %d occurred." % self.error


    def create_pending_genome(self, metadata):
        """Creates a PendingGenome object to export the data out in the turtle format with the appropriate RDF tags and
        uploads it into Blazegraph.

        Args:
            metadata(GenomeMetadata): An instance that contains metadata pertaining to an individual genome
        """
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
        """Gets the bioproject and biosample ids linked to a genome's accession ID (pertinent to NCBI sequences only)

        Args:
            metadata(GenomeMetadata): An instance that contains metadata pertaining to an individual genome
        """
        nuccore_id = return_esearch_uid("nuccore", metadata.dict["accession"])

        for item in return_elink_uid("nuccore", "bioproject", nuccore_id):
            metadata.add_genome_parameter("bioproject", item)

        for item in return_elink_uid("nuccore", "biosample", nuccore_id):
            metadata.add_genome_parameter("biosample", item)


class GenomeMetadata(object):
    """
    A container class used to store information about each individual genome before uploading to Blazegraph.
    Refactored out because the data was acting as a data clump.
    """
    def __init__(self, name, organism):
        """Initializes the class with appropriate fields and variables being created

        Args:
            name(str): The name of the genome this class represents
            organism(str): The organism from which the genome came from.
        """
        self.name = name
        self.organism = organism
        self.dict = defaultdict(set)
        self.dict['accession'].add(name)


    def add_genome_parameter(self, key, value):
        """Adds a new genome data entry to the class.

        Args:
            key(str): The type of data being added
            value(str): The value of the data
        """
        self.dict[key].add(value)


    def build_genome_kwargs(self):
        """Converts all data stored by the class into a dict used for creating a superphy.uploader.classes.Genome instance

        Returns: the kwargs dict to be used as an argument in the constructor of classes.Genome
        """
        kwargs = {'name':self.name, 'organism':self.organism}

        for key, value in self.dict.iteritems():
            if key == "serotype":
                kwargs.update(self.get_serotypes(value))

            else:
                kwargs.update({key: value})

        return kwargs

    def get_serotypes(self, serotypes):
        """Extracts information about O and H type serotype from a serotype data field

        Args:
            serotypes(str): a O???:H?? serotype string where the ??? indicates numbers of interest

        Returns: a dict containing the extracted information
        """
        Otype = None
        Htype = None

        for serotype in serotypes:
            if "ONT" in serotype or "OR" in serotype:
                Otype = "NT"
            else:
                Otype = serotype.split(":")[0][1:]

            if "HNM" in serotype:
                Htype = "NM"
            else:
                Htype = serotype.split(":")[1][1:]

        return {"Otype": Otype, "Htype": Htype}

#MetadataUploader("samples/meta_pipe_result.json", "ecoli").upload()