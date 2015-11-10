__author__ = 'Stephen Kan'

from SPARQLWrapper import SPARQLWrapper, JSON


def find_from_host(host):
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")
    queryString = (
        'PREFIX : <https://github.com/superphy#>'
        'SELECT ?p WHERE {?s ?o "%s"^^xsd:string . ?s :is_object_of ?p}' % host
    )
    sparql.setQuery(queryString)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    for result in results["results"]["bindings"]:
        return result["p"]["value"]


def find_syndrome(syndrome):
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")
    queryString = (
        'PREFIX : <https://github.com/superphy#>'
        'SELECT ?s WHERE {?s ?o "%s"^^xsd:string .}' % syndrome
    )
    sparql.setQuery(queryString)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    for result in results["results"]["bindings"]:
        return result["s"]["value"]


def find_source(source):
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")
    queryString = (
        'PREFIX : <https://github.com/superphy#>'
        'SELECT ?s WHERE {?s ?o "%s"^^xsd:string .}' % source
    )
    sparql.setQuery(queryString)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    for result in results["results"]["bindings"]:
        return result["s"]["value"]


def check_NamedIndividual(name):
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")
    queryString = (
        'PREFIX : <https://github.com/superphy#>'
        'ASK { :%s rdf:type owl:NamedIndividual .}' % name
    )
    sparql.setQuery(queryString)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    return results["boolean"]


def find_missing_sequences():
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")
    queryString = (
        "PREFIX : <https://github.com/superphy#>"
        "PREFIX gfvo: <http://www.biointerchange.org/gfvo#>"
        "PREFIX owl: <http://www.w3.org/2002/07/owl#>"
        "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>"
        "SELECT ?s ?acc WHERE { ?s rdf:type gfvo:Genome . ?s :has_accession ?acc MINUS { ?s :has_sequence ?o }}"
    )
    sparql.setQuery(queryString)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    return ((result["s"]["value"].rsplit("#", 1)[1], result["acc"]["value"]) for result in results["results"]["bindings"])


def find_unlabelled_sequences():
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")
    queryString = (
        "PREFIX : <https://github.com/superphy#>"
        "PREFIX gfvo: <http://www.biointerchange.org/gfvo#>"
        "PREFIX owl: <http://www.w3.org/2002/07/owl#>"
        "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>"
        "SELECT ?s WHERE { ?s rdf:type :Sequence . MINUS { ?s :is_from ?o }}"
    )
    sparql.setQuery(queryString)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    return (result["s"]["value"].rsplit("#", 1)[1] for result in results["results"]["bindings"])


def find_duplicate_biosamples():
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")
    queryString = (
        'PREFIX owl: <http://www.w3.org/2002/07/owl#>\n'
        'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n'
        'PREFIX : <https://github.com/superphy#>\n'
        'PREFIX gfvo: <http://www.biointerchange.org/gfvo#>\n'
        'SELECT ?BioSample (GROUP_CONCAT( ?Genome ; SEPARATOR = "#") AS ?Genomes) (COUNT (?Genome) AS ?Elements)\n'
        'WHERE { ?Genome rdf:type gfvo:Genome . ?Genome :has_biosample ?BioSample . '
        'MINUS { ?Genome :has_sequence ?Sequence . ?Sequence :is_from "WGS"^^xsd:string .}}\n'
        'GROUP BY ?BioSample HAVING ( ?Elements > 1)'
    )

    sparql.setQuery(queryString)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    return ((result["BioSample"]["value"], result["Genomes"]["value"].split("#", )[1::2])
            for result in results["results"]["bindings"])


def find_core_genomes():
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")
    queryString = (
        'PREFIX owl: <http://www.w3.org/2002/07/owl#>\n'
        'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n'
        'PREFIX : <https://github.com/superphy#>\n'
        'PREFIX gfvo: <http://www.biointerchange.org/gfvo#>\n'
        'SELECT ?BioSample ?Sequence\n'
        'WHERE { ?Genome rdf:type gfvo:Genome . ?Genome :has_biosample ?BioSample .'
        '?Genome :has_sequence ?Sequence . ?Sequence :is_from "CORE"^^xsd:string .}\n'
    )

    sparql.setQuery(queryString)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    dict = {}

    for result in results["results"]["bindings"]:
        dict[result["BioSample"]["value"]] = result["Sequence"]["value"].split("#", )[1].split("_seq")[0]

    return dict


def find_core_genome(biosample):
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")
    queryString = (
        'PREFIX owl: <http://www.w3.org/2002/07/owl#>\n'
        'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n'
        'PREFIX : <https://github.com/superphy#>\n'
        'PREFIX gfvo: <http://www.biointerchange.org/gfvo#>\n'
        'SELECT ?Genome \n'
        'WHERE { ?Genome rdf:type gfvo:Genome . ?Genome :has_biosample "%s"^^xsd:string. '
        '?Genome :has_sequence ?Sequence . ?Sequence :is_from "CORE"^^xsd:string .}\n' % biosample
    )

    sparql.setQuery(queryString)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    return [result["Genome"]["value"].split("#", 1)[1] for result in results["results"]["bindings"]]


def testing():
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")
    queryString = (
        'PREFIX owl: <http://www.w3.org/2002/07/owl#>\n'
        'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n'
        'PREFIX : <https://github.com/superphy#>\n'
        'PREFIX gfvo: <http://www.biointerchange.org/gfvo#>\n'
        # 'INSERT DATA{ :ASDFASDFTESTING rdf:type gfvo:Genome}'
        'DELETE { :ASDFASDFTESTING ?property ?value }\n'
        'WHERE { :ASDFASDFTESTING ?property ?value }'
    )

    sparql.method = 'POST'
    sparql.setQuery(queryString)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    print results

def delete_instance(name):
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")
    queryString = (
        'PREFIX owl: <http://www.w3.org/2002/07/owl#>\n'
        'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n'
        'PREFIX : <https://github.com/superphy#>\n'
        'PREFIX gfvo: <http://www.biointerchange.org/gfvo#>\n'
        'DELETE { :%s ?property ?object . ?subject ?property :%s . }\n'
        'WHERE { { :%s ?property ?object } UNION { ?subject ?property :%s } }' % (name, name, name, name)
    )

    sparql.method = 'POST'
    sparql.setQuery(queryString)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    print results

def insert_accession_sequence(core, plasmid, plasmid_seq):
    sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")
    queryString = (
        'PREFIX owl: <http://www.w3.org/2002/07/owl#>\n'
        'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n'
        'PREFIX : <https://github.com/superphy#>\n'
        'PREFIX gfvo: <http://www.biointerchange.org/gfvo#>\n'
        'INSERT DATA { :%s :has_accession "%s"^^xsd:string. '
        ':%s :has_sequence :%s . :%s :is_sequence_of :%s . }\n'
        % (core, plasmid, core, plasmid_seq, plasmid_seq, core)
    )

    sparql.method = 'POST'
    sparql.setQuery(queryString)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    print results