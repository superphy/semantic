#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""This module uploads genome metadata that has been prepared by another script
and converted into a JSON format.
ijson is used to parse the JSON file because it reads from file line by line,
instead of storing to memory first.

TODO: join the two workflows someday as to avoid file IO issues and to allow
for parallel operation and batching

Classes:
    MetadataUploader: uploads metadata from a stored JSON object
    GenomeMetadata: stores genome metadata for uploading
"""

from collections import defaultdict
import sys
import traceback
import json
import re

from ijson.backends import YAJLImportError
try:
    import ijson.backends.yajl2 as ijson
except YAJLImportError:
    import ijson.backends.yajl as ijson
from rdflib import Graph

#from Bio.Blast import NCBIXML

from superphy.upload._eutils import return_elink_uid, return_esearch_uid
from superphy.upload._sparql import check_named_individual, has_ref_gene, \
    _sparql_query
from superphy.upload._utils import generate_output, generate_path
from superphy.upload.classes import PendingGenome, Gene, GeneLocation
from superphy.upload.contig_upload import ContigUploader, ContigsWrapper
from superphy.upload.blazegraph_upload import BlazegraphUploader
'''
    # UTF-8 is necessary to handle the encoding for non-ASCII characters
    in user-inputed strings (user-inputted fields are
    # a big offender of this (e.g. location names)
'''

#reload(sys)
#sys.setdefaultencoding("utf-8")

__author__ = "Stephen Kan"
__copyright__ = """
    © Copyright Government of Canada 2012-2015. Funded by the Government of
    Canada Genomics Research and Development Initiative
    """
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Stephen Kan"
__email__ = "stebokan@gmail.com"

class MetadataUploader(object):
    """
    A class for parsing and uploading metadata from a specially prepared JSON
    file.

    TODO: link Genbank data extraction directly to MetadataUploader, so that we
    can avoid File IO issues.
    We may still need files if memory is insufficient; batch up metadata into
    files and upload files sequentially?
    """
    def __init__(self, filename):
        """Initializes the class with zeroed counters and containers.

        Args:
            filename (str): the relative filepath from this script
            organism (str): the SPARQL URI for the organism of interest being
            uploaded
        """
        self.progress = 0
        self.error = 0
        self.filename = filename

    def error_logging(self, name):
        """Records the trackback of any error messages to an log file so that
        if any are encountered, the log file will
        retain pertinent information for debugging

        Args:
            name(str): The genome that is currently being uploaded
        """
        self.error += 1
        with open(generate_path("outputs/errors.txt"), "a") as file_:
            file_.write(
                '%s \n\n %s \n '
                '================================ \n\n'
                % (name, traceback.format_exc())
            )

        print "Error %d occurred." % self.error


class Metadata(object):
    """
    A container class used to store information about an entity (such as a
    genome or gene) before uploading to Blazegraph.
    Refactored out because the data was acting as a data clump.
    """
    def __init__(self, name):
        """Initializes the class with appropriate fields and variables being
        created

        Args:
            name(str): The name of the entity
        """
        self.name = name
        self.dict = defaultdict(set)

    def add_parameter(self, key, value):
        """Adds a new genome data entry to the class.

        Args:
            key(str): The type of data being added
            value(str): The value of the data
        """
        self.dict[key].add(value)

    def build_kwargs(self):
        """Converts all data stored by the class into a dict used for creating
        a superphy.uploader.classes.(some entity) instance

        Returns: the kwargs dict to be used as an argument in the constructor
        of classes.(some entity)
        """
        kwargs = {'name':self.name}

        for key, value in self.dict.iteritems():
            kwargs.update({key: value})

        return kwargs


class GenomeMetadataUploader(MetadataUploader):
    """A class for parsing and uploading genome metadata from a specially
    prepared JSON file.

    TODO: link Genbank data extraction directly to MetadataUploader, so that we
    can avoid File IO issues.
    We may still need files if memory is insufficient; batch up metadata into
    files and upload files sequentially?
    """

    def __init__(self, filename, organism):
        super(GenomeMetadataUploader, self).__init__(filename)
        self.organism = organism

    def upload(self):
        """Uploads the contents of the given file by parsing it as an ijson
        stream.

        Prints out ending message regarding number of genomes processed and
        errors encountered
        """
        with open(generate_path(self.filename), "r") as fd:
            data = ijson.parse(fd)
            self.parse_metadata(data)

        print "%d genomes parsed, %d errors occurred." % (self.progress, self.error)

    def parse_metadata(self, data):
        """Takes in data from an ijson stream and parses for genomic metadata
        based on attached tags

        Args:
            data: an ijson stream created from an open file
        """
        metadata = None
        for prefix, event, value in data:
            if ("." not in prefix) and (prefix is not "") and \
            (event == "start_map"):
                self.add_to_graph(metadata)
                metadata = GenomeMetadata(prefix, self.organism)

            if prefix.endswith(".displayname"):
                metadata.add_parameter(prefix.split(".", 3)[1], value)
        self.add_to_graph(metadata)

    def add_to_graph(self, metadata):
        """Attempts to upload data to Blazegraph via conversion to RDF and
        turtle. Tracks errors made during this process
        as it likely has to do with either missing curated data, a blocked or
        inaccessible NCBI site, and issues with
        formatting and information availability.

        Args:
            metadata (GenomeMetadata): a GenomeMetadata object used to store
            relevant key-value pairs from the parse
        """
        if metadata:
            try:
                self.progress += 1
                print "%d: downloading files" % self.progress
                self.create_pending_genome(metadata)
            except Exception:
                self.error_logging(metadata.name)

    def create_pending_genome(self, metadata):
        """Creates a PendingGenome object to export the data out in the turtle
        format with the appropriate RDF tags and
        uploads it into Blazegraph.

        Args:
            metadata(GenomeMetadata): An instance that contains metadata
            pertaining to an individual genome
        """
        graph = Graph()
        name = metadata.name

        if check_named_individual(name):
            print "%s already in Blazegraph." % name
        else:
            self.get_ncbi_ids(metadata)
            kwargs = metadata.build_kwargs()
            PendingGenome(graph, **kwargs).rdf()
            BlazegraphUploader().upload_data(generate_output(graph))

    @classmethod
    def get_ncbi_ids(cls, metadata):
        """Gets the bioproject and biosample ids linked to a genome's accession
        ID (pertinent to NCBI sequences only)

        Args:
            metadata(GenomeMetadata): An instance that contains metadata
            pertaining to an individual genome
        """
        nuccore_id = return_esearch_uid("nuccore", metadata.dict["accession"])

        for item in return_elink_uid("nuccore", "bioproject", nuccore_id):
            metadata.add_parameter("bioproject", item)

        for item in return_elink_uid("nuccore", "biosample", nuccore_id):
            metadata.add_parameter("biosample", item)


class GenomeMetadata(Metadata):
    """
    A container class used to store information about each individual genome
    before uploading to Blazegraph.
    Refactored out because the data was acting as a data clump.
    """
    def __init__(self, name, organism):
        """Initializes the class with appropriate fields and variables being
        created

        Args:
            name(str): The name of the genome this class represents
            organism(str): The organism from which the genome came from.
        """
        super(GenomeMetadata, self).__init__(name)
        self.organism = organism
        self.dict = defaultdict(set)
        self.dict['accession'].add(name)

    def build_kwargs(self):
        """Converts all data stored by the class into a dict used for creating
        a superphy.uploader.classes.Genome instance

        Returns: the kwargs dict to be used as an argument in the constructor
        of classes.Genome
        """
        kwargs = {'name':self.name, 'organism':self.organism}

        for key, value in self.dict.iteritems():
            if key == "serotype":
                kwargs.update(self.get_serotypes(value))

            else:
                kwargs.update({key: value})

        return kwargs

    @classmethod
    def get_serotypes(cls, serotypes):
        """Extracts information about O and H type serotype from a serotype
        data field

        Args:
            serotypes(str): a O???:H?? serotype string where the ??? indicates
            numbers of interest

        Returns: a dict containing the extracted information
        """
        o_type = None
        h_type = None

        for serotype in serotypes:
            if "NT" in serotype or "R" in serotype:
                o_type = "NT"
            else:
                o_type = serotype.split(":")[0][1:]

            if "NM" in serotype:
                h_type = "NM"
            elif "NA" in serotype:
                h_type = None
            else:
                h_type = serotype.split(":")[1][1:]

        return {"Otype": o_type, "Htype": h_type}


class GeneMetadataUploader(MetadataUploader):
    """
    A class for parsing and uploading gene metadata from its XML file and
    tab-delimited data
    """

    def __init__(self, filename, kind):
        """
        Initializes GeneMetadataUploader class

        Args:
            kind (str): Describes the type of genes in the given file
            (e.g. virulence factors, antimicrobial resistance, etc)
        """
        super(GeneMetadataUploader, self).__init__(filename)
        self.kind = kind
        self.dict = {}


    def upload_genes(self):
        with open(generate_path(self.filename), "r") as file_:
            data = json.load(file_)
            if self.kind == "virulence_factor":
                self.parse_vf(data)
            elif self.kind == "antimicrobial_resistance":
                self.parse_amr(data)

    def parse_vf(self, data):
        """
        Args:
            data: an object from json.load (a dictionary) of gene metadata
        """
        for key in data:
            gene_name = key
            if "/" in gene_name: #if gene has multiple names
                gene_name = gene_name.split("/")[0]
            metadata = GeneMetadata(gene_name)
            metadata.add_parameter("gene", gene_name)
            metadata.add_parameter("gene_type", self.kind)
            for heading in data[key]:
                if data[key][heading] != "":
                    # Not a literal
                    # value = self.remove_bad_chars(str(data[key][heading]))
                    value = str(data[key][heading])
                    metadata.add_parameter(heading, value)

            self.add_to_graph(metadata)


    def parse_amr(self, data):
        """
        Parses a specially prepared amr json file for uploading amr gene Metadata
        Args:
            data: an object from json.load (a dictionary) of gene metadata
        """
        for key in data:
            if "ARO_name" in data[key]:
                gene_name = self.remove_bad_chars(str(data[key]["ARO_name"]))
                metadata = GeneMetadata(gene_name)
                metadata.add_parameter("gene", str(data[key]["ARO_name"]))
                metadata.add_parameter("gene_type", self.kind)
                for key2 in data[key]:
                    if key2 == "ARO_category":
                        for category_key in data[key][key2]:
                            metadata.add_parameter(
                                "category_uri",
                                self.remove_bad_chars(
                                    data[key][key2][category_key]["category_aro_name"]
                                )
                            )
                            metadata.add_parameter(
                                "category", data[key][key2][category_key]["category_aro_name"])
                    elif ("ARO" in key2) and not "description" in key2:
                        value = self.remove_bad_chars(str(data[key][key2]))
                        metadata.add_parameter(key2.lower(), value)
                self.add_to_graph(metadata)


    @classmethod
    def remove_bad_chars(cls, str_):
        """
        Removes / or white spaces from a string and returns the new string.

        Args:
            s(str): a string
        """
        new_str = re.sub(r"\/| |\(", "_", str_)
        new_str = re.sub(r"\.|\,|\)", "", new_str)
        new_str = re.sub(r"\'\'|\'", "_prime", new_str)
        return new_str


    def add_to_graph(self, metadata):
        """Attempts to upload data to Blazegraph via conversion to RDF and
        turtle. Tracks errors made during this process
        as it likely has to do with either missing curated data, a blocked or
        inaccessible NCBI site, and issues with
        formatting and information availability.

        Args:
            metadata (GeneMetadata): a GeneMetadata object used to store
            relevant key-value pairs from the parse
        """
        if metadata:
            try:
                self.progress += 1
                print "%d: downloading files" % self.progress
                self.create_gene(metadata)
            except Exception:
                self.error_logging(metadata.name)


    @classmethod
    def create_gene(cls, metadata):
        """
        Creates a Gene object to export the data out in the turtle format
        with the appropriate RDF tags and
        uploads it into Blazegraph.

        Args:
            metadata(GeneMetadata): An instance that contains metadata
            pertaining to an individual gene
        """
        graph = Graph()
        name = metadata.name

        if check_named_individual(name):
            print "%s already in Blazegraph." % name
        else:
            kwargs = metadata.build_kwargs()
            Gene(graph, **kwargs).rdf()
            BlazegraphUploader().upload_data(generate_output(graph))

class GeneMetadata(Metadata):
    """
    A container class used to store information about an individual gene before
    upoading to Blazegraph.
    """
    def __init__(self, name):
        """
        Initializes the class with appropriate fields and variables being
        created.

        Args:
            name(str): Name of the gene that this class represents
        """
        super(GeneMetadata, self).__init__(name)
        self.dict = defaultdict(set)
        self.dict["name"] = name

    def build_kwargs(self):
        """
        Converts all data stored by the class into a dict used for creating
        a superphy.uploader.classes.Genome instance

        Returns: the kwargs dict to be used as an argument in the constructor
        of classes.Genome
        """
        kwargs = {'name': self.name}

        for key, value in self.dict.iteritems():
            kwargs.update({key: value})

        return kwargs


#MetadataUploader("samples/meta_pipe_result.json", "ecoli").upload()

###### For Testing purposes ######

if __name__ == "__main__":
 #  # For genome testing
    # MD = GenomeMetadataUploader("samples/meta_pipe_result.json", "Human")
    # MD.upload()

 #    # For gene testing
    # GMD1 = GeneMetadataUploader('data/superphy_vf.json', "virulence_factor")
    # GMD1.upload_genes()

    GMD2 = GeneMetadataUploader('data/card.json', "antimicrobial_resistance")
    GMD2.upload_genes()
