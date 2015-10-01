__author__ = 'ubiquitin'

import requests
import os
from SPARQLWrapper import SPARQLWrapper, JSON

def upload_all_ontologies():
    faldo = "file:" + os.path.join(os.path.dirname(__file__), 'faldo.ttl')
    gfvo = "file:" + os.path.join(os.path.dirname(__file__), 'gfvo.xml')
    Superphy = "file:" + os.path.join(os.path.dirname(__file__), 'Superphy.ttl')

    ontologies = {faldo, gfvo, Superphy}

    for ontology in ontologies:
        bg_url = "http://localhost:9999/bigdata/namespace/superphy/sparql"
        data = {'uri': ontology}
        r = requests.post(url=bg_url, data=data)
        print r.content

def example_query():
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")

    sparql.setQuery("""
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX : <https://github.com/superphy/>

        SELECT ?Genome ?PropertyType ?propertyValue
        WHERE {
            ?Genome rdf:type :completed_genome .
            ?Genome ?PropertyType ?propertyValue .
        }
    """)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    for result in results["results"]["bindings"]:
        print(result["Genome"]["value"], result["PropertyType"]["value"], result["propertyValue"]["value"])

def example_modify_pending():
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")

    sparql.setQuery("""
        prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX : <https://github.com/superphy/>

        DELETE { ?Genome rdf:type :pending_genome .}
        INSERT { ?Genome rdf:type :completed_genome .}
        WHERE { ?Genome rdf:type :pending_genome . }
    """)

    sparql.setReturnFormat(JSON)
    sparql.method = 'POST' #Need this to use SPARQL UPDATE Queries
    results = sparql.query().convert()

    print results

def example_modify_completed():
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")

    sparql.setQuery("""
        prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX : <https://github.com/superphy/>

        DELETE { ?Genome rdf:type :completed_genome .}
        INSERT { ?Genome rdf:type :pending_genome .}
        WHERE { ?Genome rdf:type :completed_genome . }
    """)

    sparql.setReturnFormat(JSON)
    sparql.method = 'POST' #Need this to use SPARQL UPDATE Queries
    results = sparql.query().convert()

    print results
