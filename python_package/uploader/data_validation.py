__author__ = 'Stephen Kan'

import sparql

class Data_Validator(object):
    def __init__(self):
        pass

    def merge_biosample_duplicates(self):
        results = sparql.find_duplicate_biosamples()
        core_genomes = sparql.find_core_genomes()

        for result in results:
            (biosample, accessions) = result

