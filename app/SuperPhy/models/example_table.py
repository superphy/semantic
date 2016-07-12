#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import json
from subprocess import check_output
from SuperPhy.models.sparql.genomes import *
from SuperPhy.models.upload._sparql import *



"""
Example using the Accession = A = ["ADWR00000000", "AQFH00000000", \
 "AICG00000000", "AJLU00000000", "AJLU00000000", "AVZM00000000"] dataset

Using the test data at:
http://10.139.14.125:5000/data/testmetadata
"""


#The following will return the data from the server (in JSON)
results = check_output(["curl", "http://10.139.14.125:5000/data/testmetadata"])

#need to use .loads instead of .load as there is no read function defined for
#`results`
jsResults = json.loads(results)


#the JSON stores all the metadata in the heads,vars list
metadata = jsResults["head"]["vars"]


#we are using accession as the key to our table, so we don't want to output
#it as a column header, we want it as a row header, and to ensure that the
#data are always printed in the correct order
metadataWithoutAccession = []
headerLine = ""
for m in metadata:
    if m == "Accession":
        next
    else:
        metadataWithoutAccession.append(m)
        #generate the header of our table
        headerLine += "\t" + m


#print the header to our table
print headerLine



#the JSON stores every genome as an entry in results,bindings
jsGenomes = jsResults["results"]["bindings"]


#we want to print every genome as a new row, and the metadata in the same order
#using the accession number as the genome name
dataLine = ""
for genome in jsGenomes:
    #get the accession number of each genome, making sure no newlines exist
    acc = genome["Accession"]["value"].rstrip()
    #start the next line of output
    dataLine += acc

    #get all of the other metadata for printing
    for m in metadataWithoutAccession:
        #some values have include a newline
        #we would prefer not to have that for our table output
        nextM = genome[m]["value"].replace("\n"," ")
        dataLine += "\t" + nextM

    #end the current line, prepare for the next one
    dataLine +="\n"

print dataLine








