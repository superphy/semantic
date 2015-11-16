__author__ = 'Stephen Kan'

import requests
import os
from caller_path_gen import path

bg_url = "http://localhost:9999/bigdata/namespace/superphy/sparql"

def upload_all_ontologies():
    folder = path("ontologies")
    files = os.listdir(folder)

    for file in files:
        ontology = os.path.join(folder, file)
        upload_file(ontology)

def upload_file(filepath):
    file = "file:" + filepath
    data = {'uri': file}
    r = requests.post(bg_url, data)
    print r.content

def upload_data(data):
    headers = {'Content-Type':'application/x-turtle'}
    r = requests.post(bg_url, data=data, headers=headers)
    print r.content

def create_namespace():
    headers = {'Content-Type':'application/xml'}
    data = "".join(line for line in open(path("data/namespace.xml")))
    r = requests.post('http://localhost:9999/bigdata/namespace', data=data, headers=headers)
    return r.content