#!/usr/bin/env python
# -*- coding: UTF-8 -*-

    #
    # Example of uploading a new genome into Blazegraph.
    # We need to do three things.
    # 1. Create a new graph
    # 2. Add triples representing the genome data using the ontologies of our project
    # 3. Insert the data into the database
    #

import rdflib
from blazegraph_upload import *
from classes.sequence import *


g = rdflib.Graph()
seq = Sequence(g, "newSequence_seq", "ATCCnewGenome",(">contig1", "ATGC", ">contig2", "GGGG"), 42, 1, "fakeCheckSum", "WGS")
seq.rdf()
output = g.serialize(format='turtle')

uploader = BlazegraphUploader
uploader.upload_data(output)










