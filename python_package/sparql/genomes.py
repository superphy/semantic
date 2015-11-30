#!/usr/bin/python

from superphy import endpoint

def get_genome_meta_data(order, offset, limit):
    print "order, offset, limit"
    print order, offset, limit
    return endpoint.query("""
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX : <https://github.com/superphy#>
PREFIX gfvo: <http://www.biointerchange.org/gfvo#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

#Author: Bryce Drew
#Date: November 10

SELECT 
?Genome_Uri
?Syndrome
?Accession
?Biosample_Id
?Bioproject_Id
?Strain
?Serotype_O
?Serotype_H
?Scientific_Name
?Common_Name
?Isolation_Date
?Geographic_Location

#?Stx1_Subtype
#?Stx2_Subtype

WHERE{
  #Genome
  {
    ?Genome_Uri rdf:type gfvo:Genome .
  }
  #Bioproject_Id
  OPTIONAL{
    ?Genome_Uri :has_bioproject ?Bioproject_Id .
  }
  #Biosample_Id
  OPTIONAL{
    ?Genome_Uri :has_biosample ?Biosample_Id .
  }
  #Serotype_H
  OPTIONAL{
    ?Genome_Uri :has_Htype ?Serotype_H_Uri .
    ?Serotype_H_Uri rdfs:label ?Serotype_H .
  }
  #Serotype_O
  OPTIONAL{
    ?Genome_Uri :has_Otype ?Serotype_O_Uri .
    ?Serotype_O_Uri rdfs:label ?Serotype_O .
  }
  #Geographic_location
  OPTIONAL{ 
    ?Genome_Uri :has_geographic_location ?Geographic_Location .
  }
  #Accession number
  OPTIONAL{ 
    ?Genome_Uri :has_accession ?Accession .
  }
  #Strain
  OPTIONAL{ 
    ?Genome_Uri :has_strain ?Strain .
  }
  #Scientific & Common Name
  OPTIONAL{
    ?Genome_Uri :has_attribute ?From_Host_Uri .
    ?From_Host_Uri rdf:type :isolation_from_host .
    ?From_Host_Uri :has_attribute ?Host_Uri .
    ?Host_Uri :scientific_name ?Scientific_Name .
    ?Host_Uri :common_name ?Common_Name .
  }
  #Isolation_Date
  OPTIONAL{
    ?Genome_Uri :has_isolation_date ?Isolation_Date .
  }
  #Syndrome
  OPTIONAL{ 
    ?Genome_Uri :has_isolation_attribute ?Syndrome_Uri .
    ?Syndrome_Uri rdf:type :isolation_syndrome .  
    ?Syndrome_Uri rdfs:label ?Syndrome
  }
}
ORDER BY %s 
OFFSET %s
LIMIT %s 

""" % (order, offset, limit) )