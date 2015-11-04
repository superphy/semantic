__author__ = 'Stephen Kan'


from SPARQLWrapper import SPARQLWrapper, JSON

def find_from_host(host):
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")
    queryString = "PREFIX : <https://github.com/superphy#>" \
                  "SELECT ?p WHERE {?s ?o " + '"' + host + '"' + "^^xsd:string . ?s :is_object_of ?p}"
    sparql.setQuery(queryString)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    for result in results["results"]["bindings"]:
        return result["p"]["value"]

def find_syndrome(syndrome):
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")
    queryString = "PREFIX : <https://github.com/superphy#>" \
                  "SELECT ?s WHERE {?s ?o " + '"' + syndrome + '"' + "^^xsd:string .}"
    sparql.setQuery(queryString)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    for result in results["results"]["bindings"]:
        return result["s"]["value"]

def find_source(source):
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")
    queryString = "PREFIX : <https://github.com/superphy#>" \
                  "SELECT ?s WHERE {?s ?o " + '"' + source + '"' + "^^xsd:string .}"
    sparql.setQuery(queryString)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    for result in results["results"]["bindings"]:
        return result["s"]["value"]

def check_NamedIndividual(name):
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")
    queryString = "PREFIX : <https://github.com/superphy#>" \
                  "ASK { :"+name+" rdf:type owl:NamedIndividual .}"
    sparql.setQuery(queryString)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    return results["boolean"]

def find_missing_sequences():
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")
    queryString = "PREFIX : <https://github.com/superphy#>" \
                  "PREFIX gfvo: <http://www.biointerchange.org/gfvo#>" \
                  "PREFIX owl: <http://www.w3.org/2002/07/owl#>" \
                  "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>" \
                  "SELECT ?s ?p WHERE { ?s rdf:type gfvo:Genome . MINUS { ?s :has_sequence ?o }}"
    sparql.setQuery(queryString)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    return (result["s"]["value"].rsplit("#", 1)[1] for result  in results["results"]["bindings"])
