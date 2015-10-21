#!/usr/bin/python
#Filename: user_query_strings.py
#Author: Bryce Drew
#Date: Sept. 30, 2015
#Functionality:
	#These are sparql queries. They are designed to be specific for the user ontology:

def last_user():
	return """
	PREFIX user: <https://github.com/superphy#User>
	PREFIX RDF_type: <http://www.w3.org/1999/02/22-rdf-syntax-ns#type>
	SELECT ?s
	WHERE {
	  ?s RDF_type: user:
	}
	ORDER BY DESC(?s)
	LIMIT 1
	"""

def email_exists(email): 
	return """ASK {?x <https://github.com/superphy#hasEmail> '%s'}""" % (email)

def insert_user(sparql_id,email):
	return """
	PREFIX user: <https://github.com/superphy#User>
	PREFIX owl_NamedIndividual: <http://www.w3.org/2002/07/owl#>
	PREFIX RDF_type: <http://www.w3.org/1999/02/22-rdf-syntax-ns#type>
	PREFIX indv: <https://github.com/superphy#%s.User>
	PREFIX email: <https://github.com/superphy#hasEmail>

	INSERT DATA{
		indv: RDF_type: owl_NamedIndividual:.
		indv: RDF_type: user:.
		indv: email: '%s'
}""" % (sparql_id, email)

def insert_example_user():
	sparql_id = 0
	email = "brycedrew.test@test.com"
	return insert_user(sparql_id,email)