#rdf namespaces
namespaces = {
    'root' : 'https://www.github.com/superphy#',
    'ge' : 'http://purl.obolibrary.org/obo/GENEPIO_',
    'g' : 'http://www.biointerchange.org/gfvo#',
    'obi' : 'http://purl.obolibrary.org/obo/OBI_',
    'envo' : 'http://purl.obolibrary.org/obo/ENVO_',
    'doid' : 'http://purl.obolibrary.org/obo/DOID_',
    'faldo' : 'http://biohackathon.org/resource/faldo#',
    'ncbi' : 'http://purl.obolibrary.org/obo/NCBI_Taxon_',
    'so' : 'http://purl.obolibrary.org/obo/SO_',
    'dc' : 'http://purl.org/dc/elements/1.1/',
    'rdf' : 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'
}

#database defaults
database = {}
database['blazegraph_url'] = 'http://10.139.14.172:9000/blazegraph/namespace/superphy/sparql'
#database['blazegraph_url'] = 'http://localhost:9000/blazegraph/namespace/superphy/sparql'
