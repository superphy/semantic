#!/usr/bin/python
#Filename: endpoint.py
#Author: Bryce Drew , Stephen Kan
#Date: Sept. 29, 2015
#Functionality:
    #Functions for interacting with SPARQL endpoint. You can query and update the Db. 
    #These should be called with SPARQL queries / updates as the data parameters, and the URI of your sparql endpoint as your URL.
    #See example functions below.
#Responsibilities:
    #The caller is responsible for sanitizing the input::These functions do NOT sanitize your input.
    #The caller is responsible for formatting the data into sparql queries
    #The caller is responsible for running a db server at the url provided.

import requests
import subprocess
import logging
import re
from urllib import urlopen
from SPARQLWrapper import SPARQLWrapper, JSON
from rdflib.plugins.stores.sparqlstore import SPARQLUpdateStore
from rdflib import Graph
from rdflib.namespace import Namespace
from superphy.config import parser
from pprint import pprint



#Change this to change where the sparql endpoint is pointing
_url = "http://localhost:9999/bigdata/namespace/superphy/sparql"
_ontology_file = "file:////home/drewb/Desktop/User_Login_GraphDB/ontology/User_Ontology_RDF_XML.owl"

logger = logging.getLogger(__name__)

#Takes a sparql query and a url of your sparql endpoint. 
#Returns an object containing the information you requested. -- http://rdflib.readthedocs.org/en/latest/gettingstarted.html
def query(data):
    sparql = SPARQLWrapper(_url)
    sparql.setQuery(data)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    return results

def ask(data):
    results = query(data)
    return results['boolean']

#This function uses the requests library, as it was already implemented when we switched to RDFLib (SPARQLWrapper). It uses 'update' instead of query.
#Update changes the graph, and you should sanitize your input to avoid vandalism.

def update(data):
    payload = {'update': data}
    r = requests.post(_url, payload)
    return r.content

#Example data_location: "file:////home/ubiquitin/Documents/Ontologies/RDF Schemas for the Project/owl.ttl"
def file_update(data_location):
    payload = {'uri': data_location}
    r = requests.post(_url, payload)
    return r

#These are example uses of these functions.

#Note: Should printing really be in this file? Probably not
'''
#Quick debugging print function. Takes the object returned from bgquery, and the name of all variables you want to display from your query. Unfortunatly, you have to provide the names.
def print_query(results,*args):
    for result in results["results"]["bindings"]:
        string = ""
        for i, thing in enumerate(args):
            string = string + result[thing]['value'] + " "
        print string
'''
'''
def print_spo(results):
    print_query(results,"s","p","o")
'''
def start():
    """Starts the endpoint client"""
    subprocess.call('bash/start_blazegraph')


class SuperphySparqlStoreError(Exception):
    """RDF store error class

    """

    def __init__(Exception):
        pass


def superphy_namespace():
    """Creates RDFlib Namespace class for the superphy namespace

    """
    return Namespace('https://github.com/superphy#')


def strip_superphy_namespace(URIs):
    """Strips the 'https://github.com/superphy#' part from the uri

    Args:
        URIs(list): Can be list of strings or list of tuples/lists

    Return:
        list: Same format as input

    """

    modified = []
    for row in URIs:
        if isinstance(row, basestring):
            modified.append(re.sub(r'^https://github.com/superphy#', '', row))
        else:
            modified.append([ re.sub(r'^https://github.com/superphy#', '', u) for u in row ])

    return modified


class SuperphyStore(SPARQLUpdateStore):
    """Wrapper for RDFlib SPARQLUpdateStore 

    Initializes a RDFLib SPARQLUpdateStore connected to blazegraph API endpoint
    (Note: SPARQLUpdateStore is a subclass of SPARQLWrapper)

    Example:

    Args:



    """

    def __init__(self, *args, **kwargs):
        """constructor

        Obtains rdf_url from environment variables. Calls parent constructor
        with this enpoint url.

        Args:
            All args are passed to the SPARQLUpdateStore constructor

        """

        # Example on how to parse subclass arguments
  #       try:
     #      self._w = kwargs.pop('w')
        # except KeyError:
     #      pass


        # Retrieve superphy config
        self.config = parser.read()

        # Initialize store
        response = urlopen(self.config['rdf_url'])
        status = response.getcode()
        if status != 200:
            raise SuperphySparqlStoreError("SPARQL API Endpoint URL {} not found".format(self.config['rdf_url']))

        # Add super constructor arguments so that it is connected to endpoint
        kwargs.update({ 'queryEndpoint': self.config['rdf_url'], 'update_endpoint': self.config['rdf_url'], 
            'context_aware': False })
        
        # Call parent constructor
        super(SuperphyStore,self).__init__(*args, **kwargs)

        # Set default namespaces 
        self.bind('superphy', u'https://github.com/superphy#')

        # Set crudentials here
        # self.setCredentials(user, passwd)
        # self.setHTTPAuth('DIGEST')
        

class SuperphyGraph(object):
    """Container for RDFlib Graph 

    Initializes a RDFLib Graph with SPARQLUpdateStore connected to blazegraph API endpoint
    (Note: SPARQLUpdateStore is a subclass of SPARQLWrapper)

    Example:

    Args:



    """

    def __init__(self, logger=None):
        """


        Args:
            param1 (str): Description of `param1`.
            param2 (Optional[int]): Description of `param2`. Multiple
                lines are supported.
            param3 (List[str]): Description of `param3`.

        """

        self.logger = logger or logging.getLogger(__name__)

        # Initialize store connected to blazegraph (URL comes from ENV variable)
        self._superphyStore = SuperphyStore()

        # Initialize graph object connected to store
        self._graph = Graph(self._superphyStore)

    @property
    def graph(self):
        """
        RDF graph object

        """
        return self._graph

    @property
    def superphyStore(self):
        """
        Returns SuperphyStore object

        """
        return self._superphyStore
    
    
