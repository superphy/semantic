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
from record import HostOntology, SourceOntology, GenomeRecord

__author__ = "Matt Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "CC-SY"
__version__ = "4.0"
__maintainer__ = "Matt Whiteside"
__email__ = "matthew.whiteside@canada.ca"


class SuperphyMapperError(Exception):
    """Exceptions are documented in the same way as classes.

    """

    def __init__(Exception):
        pass


class Mapper(object):
    """Maps various metadata JSON inputs to SuperphyDB URIs


    """

    def __init__(self, logger=None):
        """


        Args:
            logger (Optional[object]): Pointer to logger op
            param3 (List[str]): Description of `param3`.

        """

        self.logger = logger or logging.getLogger(__name__)

        # Initialize object to unknown attributes encountered in input
        self.unknowns = UnknownRecord()
            
        # Initialize attribute ontology objects
        # host
        h = HostOntology()
        # source
        s = SourceOntology()

        self._ontologies = {
            'host': h,
            'source': s,
        }

        # Initialize standard parsing routines
        self._default_cleanup_routines = [
            'basic_formatting'
        ]

    
    def parse(meta_json_file, decisiontree_json_file, overrides_json_file):
        """

        """

        # Load JSON input

        # Load & validate DecisionTree JSON input

        # Load & validate Overrides JSON input

        # Iterate through inputs
        for g in inputs:
            for att in g.meta:
                pass

    


        # Output issues

        # Output reports

        pass

    def _parse_attribute():
        """
        
        """

        # Default cleanup routines, does things like strip whitespace, changes case etc

        # Custom cleanup routines are targeted for specific values
        # Stop at first successful routine

        # Find mapped superphy term and value

        # Add to genome record




        pass


    def _add_attribute(self, record, att, val):
        """

        """

        ontology_obj = self._ontologies[att]
        if ontology_obj:
            try:
                att_obj = ontology_obj.individual(val)
                setattr(record, att, val)
            except:
                # Unknown ontology value, likely needs to be added to DB

                # Record for posterity
                self.unknowns.val(att, val)

        else:
            # Unrecognized superphy attribute, probably typo in decision tree json file
            raise SuperphyMapperError("Unrecognized superphy metadata attribute %s"%att)



        pass


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
        if val in self.vals[att]:
            self.vals[att]++
        else:
            self.vals[att] = 1


    def att(self, att):
        """Increment count for unknown attribute

        Args:
            att(str): attribute name

        """
        if att in self.atts:
            self.atts[att]++
        else:
            self.atts[att] = 1


    def discarded(self, att, val):
        """Increment count for discarded attribute & value

        Args:
            att(str): attribute name
            val(str): value name

        """
        if self.discarded.get(att, {}).get(val, False):
            self.discarded[att][val]++
        else:
            self.discarded[att][val] = 1


    def write():
        """Write summary of count values


        """

        pass


