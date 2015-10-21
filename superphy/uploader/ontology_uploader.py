__author__ = 'Stephen Kan'

import requests
import os
import inspect


bg_url = "http://localhost:9999/bigdata/namespace/superphy/sparql"
ontologies = {'ontologies/faldo.ttl','ontologies/gfvo.xml','ontologies/Superphy.ttl','outputs/setup.ttl'}


def upload_all_ontologies():
    for ontology in ontologies:
        upload_file(ontology)

def upload_file(filepath):
    file = "file:" + os.path.join(os.path.dirname(inspect.getfile(inspect.currentframe())), filepath)
    data = {'uri': file}
    r = requests.post(bg_url, data)
    print r.content

def upload_data(data):
    headers = {'Content-Type':'application/x-turtle'}
    r = requests.post(bg_url, data=data, headers=headers)
    print r.content
