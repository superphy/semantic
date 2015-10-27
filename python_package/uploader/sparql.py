__author__ = 'Stephen Kan'


from SPARQLWrapper import SPARQLWrapper, JSON

sparql = SPARQLWrapper("http://localhost:9999/bigdata/namespace/superphy/sparql")

'''
sparql.setQuery("""
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX : <https://github.com/superphy#>

    SELECT ?Genome ?PropertyType ?propertyValue
    WHERE {
        ?Genome rdf:type :pending_genome .
        ?Genome ?PropertyType ?propertyValue .
    }
""")

sparql.setReturnFormat(JSON)
results = sparql.query().convert()

for result in results["results"]["bindings"]:
    print(result["Genome"]["value"], result["PropertyType"]["value"], result["propertyValue"]["value"])
'''




def find_from_host(host):
    queryString = "PREFIX : <https://github.com/superphy#>" \
                  "SELECT ?p WHERE {?s ?o " + '"' + host + '"' + "^^xsd:string . ?s :is_object_of ?p}"
    sparql.setQuery(queryString)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    for result in results["results"]["bindings"]:
        return result["p"]["value"]

def find_syndrome(syndrome):
    queryString = "PREFIX : <https://github.com/superphy#>" \
                  "SELECT ?s WHERE {?s ?o " + '"' + syndrome + '"' + "^^xsd:string .}"
    sparql.setQuery(queryString)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    for result in results["results"]["bindings"]:
        return result["s"]["value"]

def find_source(source):
    queryString = "PREFIX : <https://github.com/superphy#>" \
                  "SELECT ?s WHERE {?s ?o " + '"' + source + '"' + "^^xsd:string .}"
    sparql.setQuery(queryString)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    for result in results["results"]["bindings"]:
        return result["s"]["value"]

def check_genome(name):
    queryString = "PREFIX : <https://github.com/superphy#>" \
                  "ASK { :"+name+" rdf:type owl:NamedIndividual .}"
    sparql.setQuery(queryString)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    return results["boolean"]

"""
print find_from_host("cat")
print find_syndrome("Meningitis")
print find_source("Blood")


if check_genome("JFJL00000000"):
    print "hello"
if check_genome("asdfasdfasdf"):
    print "world"
"""