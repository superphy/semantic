#!/usr/bin/python
#Filename: user_query_strings.py
#Author: Bryce Drew
#Date: Sept. 30, 2015
#Functionality:
	#These are sparql queries. They are designed to be specific for the user ontology:

def email_exists(email): 
	return """ASK {?x <https://github.com/superphy#hasEmail> '%s'}""" % (email)

def insert_example_user():
	return """PREFIX user: <https://github.com/superphy#n000>
		PREFIX hasEmail: <https://github.com/superphy#hasEmail>
		PREFIX hasPassword: <https://github.com/superphy#hasPassword>
		PREFIX hasFirstName: <https://github.com/superphy#hasFirstName>
		PREFIX hasLastName: <https://github.com/superphy#hasLastName>
		PREFIX hasOrg: <https://github.com/superphy#hasOrganization>
		PREFIX authentic: <https://github.com/superphy#isAuthenticated>

		INSERT DATA{
		user: hasEmail: "BryceDrew.1" .
		user: hasPassword: "secretpassword" .
		user: hasFirstName: "Bryce" .
		user: hasLastName: "Drew" .
		user: hasOrg: "Laboratory for Foodborne Zoonoses Canada" .
		user: authentic: 'false' .
		user: <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2002/07/owl#NamedIndividual> .
		user: <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <https://github.com/superphy#User> .
	}"""
#TESTED!
def insert_user(node, email, password, first_name, last_name, org):
	return """PREFIX user: <https://github.com/superphy#%s>
		PREFIX hasEmail: <https://github.com/superphy#hasEmail>
		PREFIX hasPassword: <https://github.com/superphy#hasPassword>
		PREFIX hasFirstName: <https://github.com/superphy#hasFirstName>
		PREFIX hasLastName: <https://github.com/superphy#hasLastName>
		PREFIX hasOrg: <https://github.com/superphy#hasOrganization>
		PREFIX authentic: <https://github.com/superphy#isAuthenticated>

		INSERT DATA{
		user: hasEmail: "%s" .
		user: hasPassword: "%s" .
		user: hasFirstName: "%s" .
		user: hasLastName: "%s" .
		user: hasOrg: "%s" .
		user: authentic: 'false' .
		user: <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2002/07/owl#NamedIndividual> .
		user: <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <https://github.com/superphy#User> .
	}""" % (node, email, password, first_name, last_name, org)