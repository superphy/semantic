"""
makecsv.py
Created csv from json
"""

#!/usr/bin/env python
# -*- coding: UTF-8 -*-


import json
from SuperPhy.models.sparql.genomes import *
from SuperPhy.models.upload._sparql import *

class Makecsv():

    @classmethod
    def default(cls, results, extra=None):
        # need to use .loads instead of .load as there is no read function defined for
        #`results`
        response = {}
        response = json.loads(results)

        #the JSON stores all the metadata in the heads,vars list
        metadata = response["head"]["vars"]

        #we are using accession as the key to our table, so we don't want to output
        #it as a column header, we want it as a row header, and to ensure that the
        #data are always printed in the correct order
        metadataWithoutAccession = []
        csvLine = ""
        for header in metadata:
            if header == "Accession":
                next
            else:
                metadataWithoutAccession.append(header)
                #generate the headers of our table
                csvLine += "\t" + header

        #the JSON stores every genome as an entry in results,bindings
        genomeObj = response["results"]["bindings"]

        #we want to print every genome as a new row, and the metadata in the same order
        #using the accession number as the genome name

        
        csvLine += "\n" # appending a new line at the end of the headers
        
        for genome in genomeObj:
            #get the accession number of each genome, making sure no newlines exist
            acc = genome["Accession"]["value"].rstrip()
            #start the next line of output
            csvLine += acc

            #get all of the other metadata for printing
            for m in metadataWithoutAccession:
                #some values have include a newline
                #we would prefer not to have that for our table output
                nextM = genome[m]["value"].replace("\n"," ")
                csvLine += "\t" + nextM

            #end the current line, prepare for the next one
            csvLine +="\n"

        return csvLine