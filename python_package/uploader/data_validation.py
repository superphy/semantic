__author__ = 'Stephen Kan'

import sparql

class Data_Validator(object):

    def __init__(self):
        self.unmerged = set()

    def merge_biosample_duplicates(self):
            results = sparql.find_duplicate_biosamples()

            for result in results:
                (biosample, accessions) = result
                core_genome = sparql.find_core_genome(biosample)

                if len(core_genome) is 1:
                    accessions.remove(core_genome[0])

                    for accession in accessions:
                        sparql.delete_instance(accession)
                        seq = "".join([accession, "_seq"])
                        print seq
                else:
                    self.unmerged.add(biosample)

            print [entry for entry in self.unmerged]






Data_Validator().merge_biosample_duplicates()