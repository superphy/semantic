#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Module containing some utility functions to upload data and files of various
formats ont Blazegraph
"""

import os

import requests

from superphy.upload._utils import generate_path
from superphy.shared.endpoint import file_update

__author__ = "Stephen Kan"
__copyright__ = """
    © Copyright Government of Canada 2012-2015. Funded by the Government of
    Canada Genomics Research and Development Initiative
    """
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Stephen Kan"
__email__ = "stebokan@gmail.com"

class BlazegraphUploader(object):
    """
    A class that sets up data upload to Blazegraph via the initialized
    specified namespace
    """

    def __init__(self):
        pass

    @classmethod
    def upload_all_ontologies(cls):
        """
        Uploads all ontologies in the specified folder.

        The format of the ontology is automatically interpreted by Blazegraph
        based on the file extension. If any format fails, it is probably
        because of an extension mismatch (for example, Turtle files are not
        .owl as the WC3 standardized file format for RDF and OWL is RDF/XML.

        """
        folder = generate_path("ontologies")
        files = os.listdir(folder)
        for file_ in files:
            path = os.path.join(folder, file_)
            print "importing %s" % file
            file_update(path)

    @classmethod
    def upload_file(cls, filepath):
        """
            Relic - replaced by shared module
        """
        return file_update(filepath)

    @classmethod
    def upload_data(cls, data):
        """
        Uploads raw data onto Blazegraph. To ensure that Blazegraph interprets
        properly, it is necessary to specify the format in a Context-Header.

        Accepted formats are listed on this site:
        https://wiki.blazegraph.com/wiki/index.php/REST_API#MIME_Types

        Currently, the only data type needed is turtle, so this function is not
        usable for other formats.

        Args:
            data (turtle): a turtle data object

        Prints out the response object from Blazegraph
        """
        headers = {'Content-Type':'application/x-turtle'}
        request = requests.post(
            os.getenv(
                'SUPERPHY_RDF_URL',
                "http://localhost:9000/blazegraph/namespace/superphy/sparql"
            ),
            data=data,
            headers=headers
        )
        print request.content
