__author__ = 'ubiquitin'

import json
import string
from Bio import Entrez

def count_entries():

    i = 0

    with open("samples/meta_pipe_result.json") as json_file:
        json_data = json.load(json_file)

        for accession in json_data:
            i += 1

    print i

def get_nuccore_id(accession):
    Entrez.email = "stebokan@gmail.com"
    handle = Entrez.esearch(db="nuccore", retmax=5, term=accession)
    record = Entrez.read(handle)

    print ("Number of Nuccore Files: " + record["Count"])

    for id in record["IdList"]:
        print (id)

def get_DBfile(nuccore_ids):
    Entrez.email = "stebokan@gmail.com"
    handle = Entrez.efetch(db="nuccore", id=nuccore_ids, rettype="gb", retmode="xml")
    records = Entrez.parse(handle)

    for record in records:
        print "BioProject Id: " + only_digits(record["GBSeq_xrefs"][0]["GBXref_id"])
        print "BioSample Id: " + only_digits(record["GBSeq_xrefs"][1]["GBXref_id"])

def get_DBlink(nuccore_id):
    BPid = None
    BSid = None
    Entrez.email = "stebokan@gmail.com"
    handle = Entrez.elink(dbfrom="nuccore", db="bioproject", id=nuccore_id)
    records = Entrez.parse(handle)

    for record in records:
        BPid = record["LinkSetDb"][0]["Link"][0]["Id"]
        print "BioProject Id: " + record["LinkSetDb"][0]["Link"][0]["Id"]


    try:
        handle = Entrez.elink(dbfrom="nuccore", db="biosample", id=nuccore_id)
        records = Entrez.parse(handle)

        for record in records:
            print "BioSample Id: " + record["LinkSetDb"][0]["Link"][0]["Id"]
    except IndexError:
        handle = Entrez.elink(dbfrom="bioproject", db="biosample", id=BPid)
        records = Entrez.parse(handle)

        for record in records:
            print "BioSample Id: " + record["LinkSetDb"][0]["Link"][0]["Id"]

def only_digits(str):
    all = string.maketrans('','')
    nodigs = all.translate(all, string.digits)
    return str.translate(all, nodigs)


get_nuccore_id("CAFL00000000")
get_nuccore_id("JHIX00000000")
get_nuccore_id("CP001855")

get_DBfile("354810037")

get_DBlink("607502257")

get_nuccore_id("BA000007")
get_DBlink("47118301")

""" ================================================================================================================ """

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

