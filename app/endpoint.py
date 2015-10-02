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
from SPARQLWrapper import SPARQLWrapper, JSON

#Change this to change where the sparql endpoint is pointing
_url = "http://localhost:9999/bigdata/sparql" 
_ontology_file = "file:////home/drewb/Desktop/User_Login_GraphDB/ontology/User_Ontology_RDF_XML.owl"

#Takes a sparql query and a url of your sparql endpoint. 
#Returns an object containing the information you requested. -- http://rdflib.readthedocs.org/en/latest/gettingstarted.html
def query(data):
	sparql = SPARQLWrapper(_url)
	sparql.setQuery(data)
	sparql.setReturnFormat(JSON)
	results = sparql.query().convert()
	return results

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

#Quick debugging print function. Takes the object returned from bgquery, and the name of all variables you want to display from your query. Unfortunatly, you have to provide the names.
def print_query(results,*args):
	for result in results["results"]["bindings"]:
		string = ""
		for i, thing in enumerate(args):
			string = string + result[thing]['value'] + " "
		print string

def print_spo(results):
	print_query(results,"s","p","o")

#Inserts an example user into the DB

def insert_user_ontology():
	file_update(_ontology_file)
if __name__ == "__main__":
	insert_user_ontology()