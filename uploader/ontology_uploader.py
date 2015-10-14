__author__ = 'Stephen Kan'


import requests
import os


def upload_all_ontologies():
    faldo = "file:" + os.path.join(os.path.dirname(__file__), 'ontologies/faldo.ttl')
    gfvo = "file:" + os.path.join(os.path.dirname(__file__), 'ontologies/gfvo.xml')
    Superphy = "file:" + os.path.join(os.path.dirname(__file__), 'ontologies/Superphy.ttl')
    setup = "file:" + os.path.join(os.path.dirname(__file__), 'outputs/setup.ttl')


    ontologies = {faldo, gfvo, Superphy, setup}

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
