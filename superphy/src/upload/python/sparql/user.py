#!/usr/bin/python
#Filename: user_query_strings.py
#Author: Bryce Drew
#Date: Sept. 30, 2015
#Functionality:
    #These are sparql queries. They are designed to be specific for the user ontology:

from sparql.endpoint import Endpoint

def add_literal_by_email(email, predicate, literal):
    update = """
    PREFIX email: <https://github.com/superphy#hasEmail>
    PREFIX p: <%s>

    INSERT{
        ?indv p: '%s'
    }
    WHERE {
        ?indv email: '%s'
    }""" % (predicate, literal, email)
    Endpoint.update(update)
    return update

def last_user():
    return Endpoint.query("""
    PREFIX user: <https://github.com/superphy#User>
    PREFIX RDF_type: <http://www.w3.org/1999/02/22-rdf-syntax-ns#type>
    SELECT ?s
    WHERE {
      ?s RDF_type: user:
    }
    ORDER BY DESC(?s)
    LIMIT 1
    """)

def email_exists(email):
    return Endpoint.ask("""ASK {?x <https://github.com/superphy#hasEmail> '%s'}""" % (email))



def insert_user(sparql_id, email):
    return Endpoint.update("""
    PREFIX user: <https://github.com/superphy#User>
    PREFIX owl_NamedIndividual: <http://www.w3.org/2002/07/owl#>
    PREFIX RDF_type: <http://www.w3.org/1999/02/22-rdf-syntax-ns#type>
    PREFIX indv: <https://github.com/superphy#%s.User>
    PREFIX email: <https://github.com/superphy#hasEmail>

    INSERT DATA{
        indv: RDF_type: owl_NamedIndividual:.
        indv: RDF_type: user:.
        indv: email: '%s'
}""" % (sparql_id, email.lower()))

#Called by sign_up app view
def insert_next_user(email="no"):
    #See if email exists
    if email_exists(email.lower()):
        #This doesn't acutally stop anything with the sql though, so don't
        #think this is a good place to put a validator. -Bryce
        print "Already Exists!"
        return False
    #Get next user id
    results = last_user()
    user_id = 0
    for result in results["results"]["bindings"]:
        for i, thing in enumerate('s'):
            user_id = int(result[thing]['value'][28:-5]) + 1
            break
        break
    #Fill debug email for testing purposes.
    if email == "no":
        email = "test.%s@test.com" % user_id #This is here to provide an easy sample email when debugging.
    return insert_user(user_id, email.lower())
