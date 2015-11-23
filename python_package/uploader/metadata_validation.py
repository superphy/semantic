__author__ = 'Stephen Kan'

import _sparql

class MetadataValidator(object):

    def __init__(self):
        self.unmerged = set()
        self.unmergedlist = []

    def merge_biosample_duplicates(self):
            results = _sparql.find_duplicate_biosamples()

            for result in results:
                (biosample, accessions) = result
                core_genome = _sparql.find_core_genome(biosample)

                if len(core_genome) is 1:
                    accessions.remove(core_genome[0])

                    for accession in accessions:
                        _sparql.delete_instance(accession)
                        seq = "".join([accession, "_seq"])
                        _sparql.insert_accession_sequence(core_genome[0], accession, seq)

                else:
                    self.unmerged.add(biosample)
                    self.unmergedlist.append(accessions)

            print [entry for entry in self.unmerged]

if __name__ == "__main__":
    MetadataValidator().merge_biosample_duplicates()