#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
A module used to validate uploaded genome metadata on Blazegraph.

TODO: add more functionality; at the moment there is only a method to merge Plasmid and Core genome metadata from
the same BioSample ID --> may have some issues with resequencing (check dates and whatnot)
"""

import _sparql

class MetadataValidator(object):
    """
    A class for finding, resolving, and logging incidents of invalid sequences.
    """

    def __init__(self):
        """
        Sets up the class.

        """
        self.unmerged = set()
        self.unmergedlist = []

    def merge_biosample_duplicates(self):
        """
        Checks Blazegraph for any BioSample ids linked to multiple genome metadata entries and resolves them by
        joining all entries under a core genome

        ASSUMPTION: there is only one core genome per BioSample

        """
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