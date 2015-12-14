#!/usr/bin/python

from superphy import endpoint

def get_genome_meta_data():
    string = """
PREFIX  :     <https://github.com/superphy#>
PREFIX  rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX  owl:  <http://www.w3.org/2002/07/owl#>
PREFIX  rdf:  <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX  gfvo: <http://www.biointerchange.org/gfvo#>

SELECT  ?Genome_Uri 
(GROUP_CONCAT (DISTINCT ?_Syndrome ; separator=',\\n') AS ?Syndromes)(GROUP_CONCAT (DISTINCT ?_Accession ; separator=',\\n') AS ?Accession) (GROUP_CONCAT (DISTINCT ?_Biosample_Id ; separator=',\\n') AS ?Biosample_Id)(GROUP_CONCAT (DISTINCT ?_Bioproject_Id ; separator=',\\n') AS ?Bioproject_Id)(GROUP_CONCAT (DISTINCT ?_Strain ; separator=',\\n') AS ?Strain)(GROUP_CONCAT (DISTINCT ?_Serotype_O ; separator=',\\n') AS ?Serotype_O)(GROUP_CONCAT (DISTINCT ?_Serotype_H ; separator=',\\n') AS ?Serotype_H)(GROUP_CONCAT (DISTINCT ?_Scientific_Name ; separator=',\\n') AS ?Scientific_Name)(GROUP_CONCAT (DISTINCT ?_Common_Name ; separator=',\\n') AS ?Common_Name)(GROUP_CONCAT (DISTINCT ?_Isolation_Date ; separator=',\\n') AS ?Isolation_Date)(GROUP_CONCAT (DISTINCT ?_Geographic_Location ; separator=',\\n') AS ?Geographic_Location)
WHERE
  { { ?Genome_Uri rdf:type gfvo:Genome}
    OPTIONAL
      { ?Genome_Uri :has_bioproject ?_Bioproject_Id}
    OPTIONAL
      { ?Genome_Uri :has_biosample ?_Biosample_Id}
    OPTIONAL
      { ?Genome_Uri :has_Htype ?_Serotype_H_Uri .
        ?_Serotype_H_Uri rdfs:label ?_Serotype_H
      }
    OPTIONAL
      { ?Genome_Uri :has_Otype ?_Serotype_O_Uri .
        ?_Serotype_O_Uri rdfs:label ?_Serotype_O
      }
    OPTIONAL
      { ?Genome_Uri :has_geographic_location ?_Geographic_Location}
    OPTIONAL
      { ?Genome_Uri :has_accession ?_Accession}
    OPTIONAL
      { ?Genome_Uri :has_strain ?_Strain}
    OPTIONAL
      { ?Genome_Uri :has_attribute ?_From_Host_Uri .
        ?_From_Host_Uri rdf:type :isolation_from_host .
        ?_From_Host_Uri :has_attribute ?_Host_Uri .
        ?_Host_Uri :scientific_name ?_Scientific_Name . 
        ?_Host_Uri :common_name ?_Common_Name
      }
    OPTIONAL
      { ?Genome_Uri :has_isolation_date ?_Isolation_Date}
    OPTIONAL
      { ?Genome_Uri :has_isolation_attribute ?_Syndrome_Uri .
        ?_Syndrome_Uri rdf:type :isolation_syndrome .
        ?_Syndrome_Uri rdfs:label ?_Syndrome
      }
  }
GROUP BY ?Genome_Uri
ORDER BY (?Genome_Uri)
LIMIT 3
"""
    #string = string.replace("\r\n", "\n")
    return endpoint.query(string)