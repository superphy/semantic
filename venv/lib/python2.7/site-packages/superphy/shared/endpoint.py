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
import os
from SPARQLWrapper import SPARQLWrapper, JSON

#Change this to change where the sparql endpoint is pointing

#Takes a sparql query and a url of your sparql endpoint. 
#Returns an object containing the information you requested. -- http://rdflib.readthedocs.org/en/latest/gettingstarted.html
def query(data, url = os.getenv('SUPERPHY_RDF_URL', "http://localhost:9000/blazegraph/namespace/superphy/sparql") ):
	sparql = SPARQLWrapper(url)
	sparql.setQuery(data)
	sparql.setReturnFormat(JSON)
	results = sparql.query().convert()
	return results

def ask(data, url = os.getenv('SUPERPHY_RDF_URL', "http://localhost:9000/blazegraph/namespace/superphy/sparql")):
	results = query(data)
	return results['boolean']

#This function uses the requests library, as it was already implemented when we switched to RDFLib (SPARQLWrapper). It uses 'update' instead of query.
#Update changes the graph, and you should sanitize your input to avoid vandalism.

def update(data, url = os.getenv('SUPERPHY_RDF_URL', "http://localhost:9000/blazegraph/namespace/superphy/sparql")):
	payload = {'update': data}
	r = requests.post(url, payload)
	return r.content

def file_update(path):
	"""
		Will upload the file at given abolute path
	"""
	return update("""LOAD <file:%s>;""" % (path))
