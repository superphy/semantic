#!/usr/bin/python
import os
import superphy.shared.endpoint as endpoint

## Returns all genes
def get_all_genes(type="all"):
    string = ""
    if type == "vf":
        string = "?Gene rdf:type :virulence_factor ."
    elif type == "amr":
        string = "?Gene rdf:type :antimicrobial_resistance ."

    return endpoint.query("""
    PREFIX  :      <https://github.com/superphy#>
    PREFIX  rdfs:  <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX  owl:   <http://www.w3.org/2002/07/owl#>
    PREFIX  rdf:   <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX  gfvo:  <http://www.biointerchange.org/gfvo#>
    PREFIX  faldo: <http://biohackathon.org/resource/faldo#>

    SELECT  ?Gene
    (GROUP_CONCAT (DISTINCT ?_Gene_Name ; separator=',\\n') AS ?Gene_Name)
    (GROUP_CONCAT (DISTINCT ?_Category ; separator=',\\n') AS ?Category)
    (GROUP_CONCAT (DISTINCT ?_Sub_Category ; separator=',\\n') AS ?Sub_Category)
    WHERE
      { 
        { ?Gene rdf:type gfvo:gene .
          ?Gene :has_name ?_Gene_Name .
          %s}
        OPTIONAL
          { ?Gene :has_category ?_Category}
        OPTIONAL
          { ?Gene :has_sub_category ?_Sub_Category}
      }
    GROUP BY ?Gene
    ORDER BY (?Gene)
    """ % (string) , url = os.getenv('SUPERPHY_RDF_URL'))


## Returns a certain gene
def get_gene(name):
    return endpoint.query("""
    PREFIX  :      <https://github.com/superphy#>
    PREFIX  rdfs:  <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX  owl:   <http://www.w3.org/2002/07/owl#>
    PREFIX  rdf:   <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX  gfvo:  <http://www.biointerchange.org/gfvo#>
    PREFIX  faldo: <http://biohackathon.org/resource/faldo#>

    SELECT  ?Gene
    (GROUP_CONCAT (DISTINCT ?_Gene_Name ; separator=',\\n') AS ?Gene_Name)(GROUP_CONCAT (DISTINCT ?_Accession ; separator=',\\n') AS ?Accession) (GROUP_CONCAT (DISTINCT ?_Sub_Category ; separator=',\\n') AS ?Sub_Category)(GROUP_CONCAT (DISTINCT ?_Category_Id ; separator=',\\n') AS ?Category_Id)(GROUP_CONCAT (DISTINCT ?_ARO_Accession ; separator=',\\n') AS ?ARO_Accession)(GROUP_CONCAT (DISTINCT ?_ARO_Id ; separator=',\\n') AS ?ARO_Id)(GROUP_CONCAT (DISTINCT ?_VFO_Id ; separator=',\\n') AS ?VFO_Id)
    WHERE
      { 
        { ?Gene rdf:type gfvo:gene .
          ?Gene :has_name "%s"^^xsd:string .
          ?Gene :has_name ?_Gene_Name}
        OPTIONAL
          { ?Gene :has_category ?_Category_Id}
        OPTIONAL
          { ?Gene :has_has_sub_category ?_Sub_Category}
        OPTIONAL
          { ?Gene :has_vfo_id ?_VFO_Id}
        OPTIONAL
          { ?Gene :has_aro_accession ?_ARO_Accession}
        OPTIONAL
          { ?Gene :has_aro_id ?_ARO_Id}
      }
    GROUP BY ?Gene
    ORDER BY (?Gene)
    """ % (name) , url = os.getenv('SUPERPHY_RDF_URL'))

