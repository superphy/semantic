__author__ = 'Stephen Kan'

import requests
import os
import inspect

bg_url = "http://localhost:9999/bigdata/namespace/superphy/sparql"

def upload_all_ontologies():
    folder = os.path.join(os.path.dirname(inspect.getfile(inspect.currentframe())), "ontologies")
    files = os.listdir(folder)

    for file in files:
        ontology = os.path.join(folder, file)
        upload_file(ontology)

def upload_file(filepath):
    file = "file:" + filepath
    print filepath
    data = {'uri': file}
    r = requests.post(bg_url, data)
    print r.content

def upload_data(data):
    headers = {'Content-Type':'application/x-turtle'}
    r = requests.post(bg_url, data=data, headers=headers)
    print r.content
