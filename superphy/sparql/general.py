#!/usr/bin/python
#Filename: general_query_strings.py
#Author: Bryce Drew
#Date: Sept. 30, 2015
#Functionality:
	#These are sparql queries. They are designed to be general non-specific queries:

from superphy import endpoint

def delete_all_triples(): endpoint.update("""DELETE {?s?p?o} WHERE {?s?p?o}""")

def add_literal(indv): endpoint.update("""
	PREFIX user: <https://github.com/superphy#User>
	PREFIX owl_NamedIndividual: <http://www.w3.org/2002/07/owl#>
	PREFIX RDF_type: <http://www.w3.org/1999/02/22-rdf-syntax-ns#type>
	PREFIX indv: <%s>
	PREFIX email: <https://github.com/superphy#hasEmail>

	INSERT DATA{
		indv: email: '%s'
}""" % (sparql_id, email.lower()))

#The sparql statement will get all triples
def get_all_triples(): #Verified
	return endpoint.query("""SELECT * {?s ?p ?o}""")

#The sparql statement will get all triples with object literals
def get_object_literals(): #Verified
	return endpoint.query("""SELECT ?s ?p ?o WHERE {?s ?p ?o FILTER  ISLITERAL(?o)}""")

#The sparql statement will get all tripples with at least 1 literal
def get_all_literals():
	return endpoint.query("""
	SELECT ?s ?p ?o
	WHERE {d
	  {?s ?p ?o
	   FILTER  ISLITERAL(?o)
	  }
	  UNION
	  {?s ?p ?o
	   FILTER  ISLITERAL(?p)
	  }
	  UNION
	  {?s ?p ?o
	   FILTER  ISLITERAL(?s)
	  }
	}""")

#Tested
#The sparql statement will get all triples with 3 URIs.
def get_all_uri_triples():
	return endpoint.query("""
	SELECT ?s ?p ?o
	WHERE {
	  {?s ?p ?o}
	  MINUS
	  {?s ?p ?o
	    FILTER ISLITERAL(?o)
	  }
	  MINUS
	  {
	    FILTER ISLITERAL(?p)
	  }
	  MINUS
	  {
	    FILTER ISLITERAL(?s)
	  }
	}""")