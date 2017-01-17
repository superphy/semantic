"""
#!/usr/bin/python
#Filename: endpoint.py
#Author: Bryce Drew , Stephen Kan
#Date: Sept. 29, 2015
#Functionality:
    #Functions for interacting with SPARQL endpoint. You can query and update
    the Db.
    #These should be called with SPARQL queries / updates as the data
    parameters, and the URI of your sparql endpoint as your URL.
    #See example functions below.
#Responsibilities:
    #The caller is responsible for sanitizing the input::These functions do NOT sanitize your input.
    #The caller is responsible for formatting the data into sparql queries
    #The caller is responsible for running a db server at the url provided.
"""

import os
import requests
from SPARQLWrapper import SPARQLWrapper, JSON

class Endpoint(object):
    """
    Interface between blazegraph and python
    """
    DEFAULT = "http://10.139.14.172:9000/blazegraph/namespace/superphy/sparql"
    @classmethod
    def query(cls, data, url=os.getenv('SUPERPHY_RDF_URL', DEFAULT)):
        """
        Sends the data passed in to the endpoint indicated as a sparql query.
        #Takes a sparql query and a url of your sparql endpoint.
        #Returns an object containing the information you requested.
        -- http://rdflib.readthedocs.org/en/latest/gettingstarted.html
        """
        sparql = SPARQLWrapper(url)
        sparql.setQuery(data)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        return results

    @classmethod
    def ask(cls, data):
        """
        Handles ask requests to the endpoint. This is different than queries,
        and it's return type is always boolean. You have to give it a valid
        'ask' query. You can't just give this function a SELECT ... query.
        """
        results = cls.query(data)
        return results['boolean']

    @classmethod
    def update(cls, data, url=os.getenv('SUPERPHY_RDF_URL', DEFAULT)):
        """
        #This function uses the requests library, as it was already
        implemented when we switched to RDFLib (SPARQLWrapper). It uses
        'update' instead of query. Update changes the graph, and you should
        sanitize your input to avoid vandalism.
        """
        payload = {'update': data}
        r = requests.post(url, payload)
        return r.content
    @classmethod
    def file_update(cls, path):
        """
            Will upload the file at given abolute path
        """
        return cls.update("""LOAD <file:%s>;""" % (path))
