__author__ = 'ubiquitin'


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

sparql.setQuery("""
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX : <https://github.com/superphy#>

    SELECT ?p
    WHERE {
        ?s ?o "cat"^^xsd:string .
        ?s :is_object_of ?p
    }
""")

sparql.setReturnFormat(JSON)
results = sparql.query().convert()

for result in results["results"]["bindings"]:
    print result["p"]["value"]

