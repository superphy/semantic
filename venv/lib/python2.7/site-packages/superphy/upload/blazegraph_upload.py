#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Module containing some utility functions to upload data and files of various formats ont Blazegraph
"""

import os

import requests

from _utils import generate_path
from superphy.shared.endpoint import file_update

__author__ = "Stephen Kan"
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Stephen Kan"
__email__ = "stebokan@gmail.com"

class BlazegraphUploader(object):
    """
    A class that sets up data upload to Blazegraph via the initialized specified namespace
    """

    def __init__(self):
        pass

    def upload_all_ontologies(self):
        """
        Uploads all ontologies in the specified folder.

        The format of the ontology is automatically interpreted by Blazegraph based on the file extension. If any
        format fails, it is probably because of an extension mismatch (for example, Turtle files are not .owl as the
        WC3 standardized file format for RDF and OWL is RDF/XML.

        """
        folder = generate_path("ontologies")
        files = os.listdir(folder)
        for file in files:
            path = os.path.join(folder, file)
            print "importing %s" % file
            file_update(path)

    def upload_file(self, filepath):
        """
            Relic - replaced by shared module
        """
        return endpoint.file_update(filepath)

    def upload_data(self, data):
        """Uploads raw data onto Blazegraph. To ensure that Blazegraph interprets properly, it is necessary to specify
        the format in a Context-Header.

        Accepted formats are listed on this site: https://wiki.blazegraph.com/wiki/index.php/REST_API#MIME_Types

        Currently, the only data type needed is turtle, so this function is not usable for other formats.

        Args:
            data (turtle): a turtle data object

        Prints out the response object from Blazegraph
        """
        headers = {'Content-Type':'application/x-turtle'}
        r = requests.post(os.getenv('SUPERPHY_RDF_URL', "http://localhost:9000/blazegraph/namespace/superphy/sparql"), data=data, headers=headers)
        print r.content
