#!/usr/bin/python
import os
import superphy.shared.endpoint as endpoint

from prefixes import prefixes

## Returns all genes
def get_all_genes(type_="all"):
    gene_type = ""
    if type_ == "vf":
        gene_type = "?Gene rdf:type :virulence_factor ."
    elif type_ == "amr":
        gene_type = "?Gene rdf:type :antimicrobial_resistance ."

    query = prefixes + """
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
    """ % (gene_type)

    return endpoint.query(query)


## Returns a certain gene
def get_gene(name):
    query = prefixes + """
    SELECT  ?Gene
    (GROUP_CONCAT (DISTINCT ?_Gene_Name ; separator=',\\n') AS ?Gene_Name)
    (GROUP_CONCAT (DISTINCT ?_Accession ; separator=',\\n') AS ?Accession)
    (GROUP_CONCAT (DISTINCT ?_Sub_Category ; separator=',\\n') AS ?Sub_Category)
    (GROUP_CONCAT (DISTINCT ?_Category_Id ; separator=',\\n') AS ?Category_Id)
    (GROUP_CONCAT (DISTINCT ?_ARO_Accession ; separator=',\\n') AS ?ARO_Accession)
    (GROUP_CONCAT (DISTINCT ?_ARO_Id ; separator=',\\n') AS ?ARO_Id)
    (GROUP_CONCAT (DISTINCT ?_VFO_Id ; separator=',\\n') AS ?VFO_Id)
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
    """ % (name)

    return endpoint.query(query)


## Returns the instances of a particular gene in a genome
def find_regions(gene, genome):
    query = prefixes + """
    SELECT  ?Region 
    WHERE
      { 
        { ?Region rdf:type faldo:Region .
          ?Gene :has_copy ?Region .
          ?Contig :has_gene ?Region .
          ?Contig :is_contig_of ?Genome .
          ?Gene :has_name "%s"^^xsd:string . 
          #?Genome :has_accession "%s"^^xsd:string . 
          }
      }
    """ % (gene, genome)
    return endpoint.query(query)

def get_regions(genome_list, gene_list):
    genes = ""
    genomes = ""
    for gene in gene_list:
        genes = genes + ' "%s" ' % gene
    for genome in genome_list:
        genomes = genomes + ' :%s ' % genome
    query = prefixes + """
    SELECT ?Region ?Gene_Name ?Genome
    WHERE
      {
        VALUES ?Genome { %s }
        VALUES ?Gene_Names { %s }
        ?Region a faldo:Region .
        ?Gene :has_copy ?Region ;
              :has_name ?Gene_Name .
        ?Contig :has_gene ?Region ;
                :is_contig_of ?Genome .
        FILTER (?Gene_Name = ?Gene_Names)
    }
    """ % (genomes, genes)

    print "query is:", query

    return endpoint.query(query)

def get_categories(type):
    query = prefixes
    if type == "vf":
        query += """
        SELECT ?Category ?Subcategory
        WHERE {
            ?cat rdfs:subClassOf vfo:category .
            ?subcat rdfs:subClassOf ?cat .
            ?cat rdfs:label ?Category .
            ?subcat rdfs:label ?Subcategory .
            FILTER (?cat != ?subcat && ?cat != vfo:category) .
            MINUS { ?cat owl:equivalentClass ?equivclass . }
        }
        """
    elif type == "amr":
        ## Temporarily loads vf data until amr categories are figured out
        query += """
        SELECT ?Category ?Subcategory
        WHERE {
            ?cat rdfs:subClassOf vfo:category .
            ?subcat rdfs:subClassOf ?cat .
            ?cat rdfs:label ?Category .
            ?subcat rdfs:label ?Subcategory .
            FILTER (?cat != ?subcat && ?cat != vfo:category) .
            MINUS { ?cat owl:equivalentClass ?equivclass . }
        }
        """
    else:
        raise ValueError("Non-valid gene type inserted.")

    return endpoint.query(query)
