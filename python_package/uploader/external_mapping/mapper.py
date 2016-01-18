#!/usr/bin/env python

"""mapper

Converts input metadata JSON input with various metadata type-value pairs
to URIs used in the Superphy DB.

Example:
    python mapper.py decision

Rules for assigning foreign metadata key/value pairs to superphy URIs are defined in decision_tree JSON file.
This decision_tree structure list specific functions in cleanup_routines.py that santizes the input (maps synonyms to a single
term, formats, etc) and specific functions in validation_routines.py that can recognize the targeted foreign data terms values
 and can return superphy URIs.

Decision Tree JSON format:

    {
        'foreign_term1': {
            'keep': 0
        },
        {
            'keep': 1,
            'cleanup_routines': ['routine1','routine2'],
            'validation_routines': ['routine3']
        }
    }

Input metadata JSON format:

    {
        'genome_accession1': { 
            'term1': ['string1', 'string2'],
            'term2': 'string'
        },
        'genome_accession2': { 
            'term1': ['string1', 'string2'],
            'term2': 'string'
        }
    }

Occasionally there will be conflicts in foreign meta-data (e.g. suggestions of multiple hosts, incompatible source and host).
Overrides JSON file allows the user to define correct annotations for specific genomes in order to resolive the conflicts.

Overrides JSON format:




"""

import logging
import json
import argparse
import itertools
import sys
from superphy.endpoint import SuperphyStore
from ontology import HostOntology, SourceOntology, LocationOntology, GenomeRecord
import cleanup_routines
import validation_routines

# REMOVE not needed in production
from pprint import pprint

__author__ = "Matt Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Matt Whiteside"
__email__ = "matthew.whiteside@canada.ca"


class SuperphyMapperError(Exception):
    """Errors encountered in workflow e.g. malformed input

    """
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
        self.store = SuperphyStore()
        # host
        h = HostOntology(self.store, self.logger)
        # source
        s = SourceOntology(self.store, self.logger)
        # location
        l = LocationOntology(self.store, self.logger)

        self._ontologies = {
            'host': h,
            'source': s,
            'location': l
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
            try:     
                self.decision_tree = json.load(data_file)
                #self.decision_tree = self._byteify(self.decision_tree)
            except Exception, m:
                self.logger.error("Error encountered during JSON parsing: %s"%m)
                raise SuperphyMapperError("Invalid decision tree. Check log for details.")

        
        if not self._valid_dt():
            raise SuperphyMapperError("Invalid decision tree. Check log for details.")
            
        # Load & validate Overrides JSON input
    
    def _byteify(self, input):
        """Convert utf-8 input to strings

        """
        if isinstance(input, dict):
            return { self._byteify(key): self._byteify(value) for key,value in input.iteritems() }
        elif isinstance(input, list):
            return [ self._byteify(element) for element in input ]
        elif isinstance(input, unicode):
            return input.encode('utf-8')
        else:
            return input


    def ontology(self, name):
        """Get method for ontology object

        Args:
            name(str): ontology dict key name

        Returns:
            AttributeOntology subclass

        """

        return self._ontologies[name]

    
    def run(self, meta_json_file, outfile):
        """Top-level method that runs foreign-to-superphy mapping
        algorithm

        Args:
            meta_json_file(str): filepath to input JSON file
            outfile(str): filepath to output JSON file

        Returns:
            bool: True if successful


        Input format:
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
            #meta = self._byteify(meta)

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


        if self.unknowns.unresolved_terms():
            self.logger.info(self.unknowns.summary())
            m = "Unknown values found in input. These must be handled"
            self.logger.warn(m)
            raise SuperphyMapperError(m)
            return False


        # Apply overrides that resolve issues like multiple hosts
        # for specific genomes
        #self.overrides()

        # Output issues
        all_ok = True
        for g in genomes:
            valid, message = g.is_valid()

            if not valid:
                m = "Genome %s failed checks:\n\t%s"%(g.uniquename(), message)
                self.logger.warn(m)
                raise SuperphyMapperError(m)


        if not all_ok:
            return False


        # Output foreign-to-superphy mapping
        mapping = {}
        for g in genomes:
            mapping[g.uniquename()] = g.output()

        # Convert to JSON
        with open(outfile, 'w') as outfile:
            json.dump(mapping, outfile)

        return True


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
                    self.logger.debug("Trying validation routine: %s on value: %s"%(m,clean_val))
                    superphy_tuples = getattr(validation_routines, m)(clean_val, self)

                    if superphy_tuples:
                        self.logger.debug("Validation routine result: %s"%superphy_tuples)
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
                                    except Exception:
                                        # Unknown ontology value, likely needs to be added to DB
                                        self.logger.exception("Creation of annotation record for genome %s and term %s/%s failed."%(genome.uniquename(), superphy_term, superphy_value))

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
                    # Value is not matched by any validation method
                    self.unknowns.val(att, clean_val)
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

        for term, tree in self.decision_tree.iteritems():

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
        discarded(dict): Counts of which attributes have been entirely discarded
        skipped(dict): Counts of discarded values in attributes that were at least partially parsed

    """

    def __init__(self):
        """Constructor

       
        Args:
            None

        """
        self.vals = {}
        self.atts = {}
        self.discarded = {}
        self.skipped = {}


    def val(self, att, val):
        """Increment count for unknown value

        Args:
            att(str): attribute name
            val(str): value name

        """

        if att in self.vals:
            if val in self.vals[att]:
                self.vals[att][val] += 1
            else:
                self.vals[att][val] = 1
        else:
            self.vals[att] = { val: 1 }


    def att(self, att, val):
        """Increment count for unknown attribute

        Args:
            att(str): attribute name
            val(str): value name

        """
        
        if att in self.atts:
            if val in self.atts[att]:
                self.atts[att][val] += 1
            else:
                self.atts[att][val] = 1
        else:
            self.atts[att] = { val: 1 }



    def discarded(self, att, val):
        """Increment count for discarded attribute & value

        Args:
            att(str): attribute name
            val(str): value name

        """

        if att in self.discarded:
            if val in self.discarded[att]:
                self.discarded[att][val] += 1
            else:
                self.discarded[att][val] = 1
        else:
            self.discarded[att] = { val: 1 }


    def skipped(self, att, val):
        """Increment count for skipped value in an attribute

        Args:
            att(str): attribute name
            val(str): value name

        """

        if att in self.skipped:
            if val in self.skipped[att]:
                self.skipped[att][val] += 1
            else:
                self.skipped[att][val] = 1
        else:
            self.skipped[att] = { val: 1 }


    def unresolved_terms(self):
        """Unresolved terms/values are foreign metadata terms/values
        that have been encountered in the input but have no defined
        rules in the decision_tree object

        Returns:
            bool: True if unknown terms/values have been recorded

        """

        if self.atts or self.vals:
            return True

        return False


    def summary(self):
        """Generate summary of unresolved terms

        Returns:
            str: Summary text

        """

        tot = 0
        summary_string = "\n\nList of unresolved foreign data.  These must be dealt with before Mapper.run() can successfully complete.\n"
        summary_string += "\nUnresolved foreign attributes:\n-------------------------------\n"

        n = 0
        if self.atts:
            for a in self.atts:
                for v in self.atts[a]:
                    n += 1
                    summary_string += "%i. %s - %s (%i occurences)\n"%(n, a, v, self.atts[a][v])
                    

        else:
            summary_string += "None\n"

        tot += n

        summary_string += "\nUnresolved foreign values:\n-------------------------------\n(known attribute, unknown value)\n"
        n = 0
        if self.vals:
            for a in self.vals:
                for v in self.vals[a]:
                    n += 1
                    summary_string += "%i. %s - %s (%i occurences)\n"%(n, a, v, self.vals[a][v])
        else:
            summary_string += "None\n"

        tot += n
        summary_string += "\nTOTAL UNRESOLVED FOREIGN TERMS & VALUES: %i\n"%tot

        summary_string += "\nList of intentially discarded foreign terms and values.  Only for reference.\n"
        summary_string += "\nDiscarded attributes:\n-------------------------------\n"

        n = 0
        if self.discarded:
            for a in self.discarded:
                for v in self.discarded[a]:
                    n += 1
                    summary_string += "%i. %s - %s (%i occurences)\n"%(n, a, v, self.discarded[a][v])
                    

        else:
            summary_string += "None\n"

        summary_string += "\nDiscarded values:\n-------------------------------\n(entire attribute is not discarded, only certain values)\n"
        n = 0
        if self.skipped:
            for a in self.skipped:
                for v in self.skipped[a]:
                    n += 1
                    summary_string += "%i. %s - %s (%i occurences)\n"%(n, a, v, self.skipped[a][v])
        else:
            summary_string += "None\n"


        return summary_string



if __name__ == "__main__":
    """Run mapper with given inputs

    """
    #logging.basicConfig(level=logging.DEBUG, stream=sys.stderr)
    logging.basicConfig(level=logging.DEBUG)

    # Parse command-line args
    parser = argparse.ArgumentParser(description='Foreign metadata mapping to Superphy ontology terms')
    parser.add_argument('decision_tree_json', action="store")
    parser.add_argument('overrides_json', action="store")
    parser.add_argument('meta_json', action="store")
    parser.add_argument('output', action="store")

    files = parser.parse_args()

    mapper = Mapper(files.decision_tree_json, files.overrides_json)

    mapper.run(files.meta_json, files.output)