#!/usr/bin/python

import superphy.shared.endpoint as endpoint

## Returns all virulence factors
def virulence_factors():
	return endpoint.query("""
	PREFIX  :     <https://github.com/superphy#>
	PREFIX  rdfs: <http://www.w3.org/2000/01/rdf-schema#>
	PREFIX  owl:  <http://www.w3.org/2002/07/owl#>
	PREFIX  rdf:  <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
	PREFIX  gfvo: <http://www.biointerchange.org/gfvo#>

	SELECT ?vf 
	WHERE {
		?vf rdf:type gfvo:gene .
		?vf rdf:type :virulence_factor .
	}""")