__author__ = 'ubiquitin'


import requests
import os


def upload_all_ontologies():
    faldo = "file:" + os.path.join(os.path.dirname(__file__), 'ontologies/faldo.ttl')
    gfvo = "file:" + os.path.join(os.path.dirname(__file__), 'ontologies/gfvo.xml')
    Superphy = "file:" + os.path.join(os.path.dirname(__file__), 'ontologies/Superphy.ttl')
    results = "file:" + os.path.join(os.path.dirname(__file__), 'results.ttl')


    ontologies = {faldo, gfvo, Superphy, results}

    for ontology in ontologies:
        bg_url = "http://localhost:9999/bigdata/namespace/superphy/sparql"
        data = {'uri': ontology}
        r = requests.post(url=bg_url, data=data)
        print r.content

def upload_ontology(filepath):
    ontology = "file:" + os.path.join(os.path.dirname(__file__), filepath)
    bg_url = "http://localhost:9999/bigdata/namespace/superphy/sparql"
    data = {'uri': ontology}
    r = requests.post(url=bg_url, data=data)
    print r.content

#upload_all_ontologies()
upload_ontology("results.ttl")