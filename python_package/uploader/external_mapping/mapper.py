#!/usr/bin/env python

"""mapper

Converts input metadata JSON input with varisous metadata type-value pairs
to URIs used in the Superphy DB

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.


"""

import logging
import json
import argparse
import itertools
from pprint import pprint
from superphy.endpoint import SuperphyGraph
from ontology import HostOntology, SourceOntology, GenomeRecord
import cleanup_routines
import validation_routines

__author__ = "Matt Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "CC-SY"
__version__ = "4.0"
__maintainer__ = "Matt Whiteside"
__email__ = "matthew.whiteside@canada.ca"


class SuperphyMapperError(Exception):
    """Exceptions are documented in the same way as classes.

    """

    def __init__(Exception, m):
        pass


class Mapper(object):
    """Maps various metadata JSON inputs to SuperphyDB URIs


    """

    def __init__(self, decisiontree_json_file, overrides_json_file, logger=None):
        """


        Args:
            logger (Optional[object]): Pointer to logger op
            param3 (List[str]): Description of `param3`.

        """

        self.logger = logger or logging.getLogger(__name__)

        # Initialize object to track unknown attributes encountered in input
        self.unknowns = UnknownRecord()
            
        # Initialize attribute ontology objects
        self.graph = SuperphyGraph().graph
        # host
        h = HostOntology(self.graph)
        # source
        s = SourceOntology(self.graph)

        self._ontologies = {
            'host': h,
            'source': s,
        }

        # Initialize standard parsing routines
        self._default_cleanup_routines = [
            'basic_formatting'
        ]

        self._default_validation_routines = [
            'non_value'
        ]

        # Load & validate DecisionTree JSON input
        with open(decisiontree_json_file) as data_file:    
            self.decision_tree = json.load(data_file)
        
        if not self._valid_dt():
            raise SuperphyMapperError("Invalid decision tree. Check log for details.")
            
        # Load & validate Overrides JSON input
    

    def ontology(self, name):
        """Get method for ontology object

        Args:
            name(str): ontology dict key name

        Returns:
            AttributeOntology subclass

        """

        return self._ontologies[name]

    
    def map(self, meta_json_file):
        """

        Format:
            {
                genome_accession1: { 
                    term1: ['string1', 'string2'],
                    term2: 'string'
                },
                genome_accession2: { 
                    term1: ['string1', 'string2'],
                    term2: 'string'
                }
            }


        """

        genomes = []

        # Load JSON metadata input
        meta = {}
        with open(meta_json_file) as data_file:
            meta = json.load(data_file)

        if not isinstance(meta, dict):
            raise SuperphyMapperError("Invalid input meta file. Expecting JSON object.")

    
        # Iterate through genome entries
        for g_key in meta:

            # Start genome record
            g = GenomeRecord(g_key)

            # Parse genome attributes
            for att in meta[g_key]:
                val = meta[g_key][att]
                if isinstance(val, list):
                    for v in val:
                        self.assign(g, att, v)
                else:
                    self.assign(g, att, val)

            # Save genome
            genomes.append(g)


        # Output issues

        # Output reports

        pass

    def assign(self, genome, att, val):
        """Using decision tree rules, assign term-value pair to superphy ontology

        Args:
            genome(GenomeRecord): object for recording genome assignments
            att(str): foreign metadata term
            val(str): foreign metadata value

        Returns:
            bool: True if assign was successful
        
        """

        # Skip empty values
        if not val:
            return True

        # Retrieve decision routines
        dt = self.decision_tree[att]

        if not dt:
            # Encountered unknown term which is not handled in decision tree file
            self.unknowns.att(att, val)
            return False

        else:
            if dt['keep']:
                # This foreign metadata term is a keeper
                
                # Run value through all standard cleanup routines
                clean_val = val 
                for m in self._default_cleanup_routines:
                    tmp, clean_val = getattr(cleanup_routines, m)(clean_val)

                # Run value through custom cleanup routines, stopping at first one
                # that is successful
                for m in dt['cleanup_routines']:
                    success, clean_val = getattr(cleanup_routines, m)(clean_val)
                    if success:
                        break

                self.logger.debug("Cleanup routines turned %s into %s"%(val, clean_val))

                # Find matching Superphy term & value for this attribute-value pair
                # Default validation routines do things like 'skip' over missing values
                assigned = False
                vr = itertools.chain(self._default_validation_routines, dt['validation_routines'])
                for m in vr:
                    pprint( m)
                    superphy_tuples = getattr(validation_routines, m)(clean_val, self)

                    if superphy_tuples:
                        if superphy_tuples != 'skip':
                            # Store one or more superphy key-value pairs

                            for superphy_term,superphy_value in superphy_tuples:
                                # Add each assigned key-value pair to genome record

                                ontology_obj = self._ontologies[superphy_term]
                                if ontology_obj:
                                    try:
                                        att_obj = ontology_obj.individual(superphy_value)
                                        getattr(genome, superphy_term)(att_obj)
                                        assigned = True
                                    except:
                                        # Unknown ontology value, likely needs to be added to DB

                                        # Record for posterity
                                        self.unknowns.val(att, clean_val)

                                else:
                                    # Unrecognized superphy attribute, probably typo in decision tree json file
                                    raise SuperphyMapperError("Unrecognized superphy metadata attribute %s"%superphy_term)

                        else:
                            self.unknowns.skipped(att, val)
                        break

                # Was foreign attribute assigned to at least one superphy key-value pair?
                if assigned:
                    return True
                else:
                    return False
                
            else:
                # Record all discarded terms
                self.unknowns.discarded(att, val)
                return True





    def _valid_dt(self):
        """Checks decision tree dict has correct format and
        uses defined methods.

        Format:
            term1: { 
                keep: 1,
                cleanup_routines: []
                validation_routines: []
            },
            term2: {
                keep: 0
            }

        Args:
            decision_tree(dict): Dictionary loaded from JSON file

        Returns:
            bool: True if format is valid

        """

        for term, tree  in self.decision_tree.iteritems():

            recognized_fields = [(cleanup_routines, 'cleanup_routines'), (validation_routines, 'validation_routines')]
            
            if tree['keep']:
                # This term will be parsed

                # Check that block has required fields and that methods in the cleanup and validation lists are defined
                for module,field in recognized_fields:
                    if field in tree:
                        if isinstance(tree[field], list):
                            for m in tree[field]:
                                if not hasattr(module, m):
                                    self.logger.warn("Invalid decision tree format: attribute %s field %s contains an undefined method %s."%(term,field,m))
                                    return False

                        else:
                            self.logger.warn("Invalid decision tree format: attribute %s field %s is not an array."%(term,field))
                            return False

                    else:
                        self.logger.warn("Invalid decision tree format: attribute %s missing required field %s."%(term,field))
                        return False

            else:
                # This term will be discarded
                # The block should be empty, if its not warn user

                for f in recognized_fields:
                    if f in tree:
                        self.logger.warn("Invalid decision tree format: discarded attribute %s should be empty."%term)
                        return False


        return True

    # def _valid_m(self, meta):
    #     """Checks meta dict has correct format

    #     Format:
    #         {
    #             genome_accession1: { 
    #                 term1: [],
    #                 term2: 'string'
    #             },
    #             genome_accession2: { 
    #                 term1: [],
    #                 term2: 'string'
    #             }
    #         }

    #     Args:
    #         meta(dict): Dictionary loaded from JSON file

    #     Returns:
    #         bool: True if format is valid
    
    #     """

    #     if not isinstance(meta, 'dict'):
    #         self.logger.warn("Invalid input meta file. Expecting JSON object.")
    #         return False

    #     for g in meta:



    #     return True 



class UnknownRecord(object):
    """Counts unknown metadata attributes and values, as well
    as discarded values


    Attributes:
        vals(dict): Counts of unknown values assigned to known attributes
        atts(dict): Counts of unknown attributes
        discarded(dict): Counts of which attributes have been discarded

    """

    def __init__(self):
        """Constructor

       
        Args:
            None

        """
        self.vals = {}
        self.atts = {}
        self.discarded = {}


    def val(self, att, val):
        """Increment count for unknown value

        Args:
            att(str): attribute name
            val(str): value name

        """
        if self.vals.get(att, {}).get(val, False):
            self.vals[att][val] += 1
        else:
            self.vals[att][val] = 1


    def att(self, att, val):
        """Increment count for unknown attribute

        Args:
            att(str): attribute name
            val(str): value name

        """
        if self.atts.get(att, {}).get(val, False):
            self.atts[att][val] += 1
        else:
            self.atts[att][val] = 1


    def discarded(self, att, val):
        """Increment count for discarded attribute & value

        Args:
            att(str): attribute name
            val(str): value name

        """
        if self.discarded.get(att, {}).get(val, False):
            self.discarded[att][val] += 1
        else:
            self.discarded[att][val] = 1


    def write():
        """Write summary of count values


        """

        pass



if __name__ == "__main__":
    """Run mapper with given inputs

    """
    logging.basicConfig(level=logging.WARNING)

    # Parse command-line args
    parser = argparse.ArgumentParser(description='Foreign metadata mapping to Superphy ontology terms')
    parser.add_argument('decision_tree_json', action="store")
    parser.add_argument('overrides_json', action="store")
    parser.add_argument('meta_json', action="store")

    files = parser.parse_args()

    mapper = Mapper(files.decision_tree_json, files.overrides_json)

    mapper.map(files.meta_json)