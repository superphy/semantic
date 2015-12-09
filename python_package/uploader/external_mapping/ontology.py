#!/usr/bin/env python

"""Container objects for temporarily storing metadata.

Classes store metadata and verify correctness and possible conflicts in metadata

Example:
    Usage:

        $ python example_google.py


"""

import abc
import itertools
from pprint import pformat
from superphy.endpoint import strip_superphy_namespace

# REMOVE not needed in production
from pprint import pprint


__author__ = "Matt Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Matt Whiteside"
__email__ = "matthew.whiteside@canada.ca"


def make_unicode(input):
    """Encode byte string in utf-8

    Strings coming from json will be in utf-8 format
    encode uri's in unicode type to match

    """
    if type(input) != unicode:
        input =  input.decode('utf-8')
        return input
    else:
        return input


class SuperphyMetaError(Exception):
    """Conflicts or missing metadata

    """

    def __init__(Exception):
        pass


class AttributeOntology(object):
    """Abstract class that defines interface for storing sections of the Superphy metadata ontology


    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, graph):
        """Constructor

        Args:
            graph (object): SupephyGraph object

        """
        self._ontology = self._initialize_ontology(graph)
        self._attribute_class = 'Attribute'

        pass


    @abc.abstractproperty
    def attribute_class(self):
        """Attribute class represented by ontology

        """
        return 'Should never get here'


    @abc.abstractmethod
    def _initialize_ontology(self, graph):
        """Retrieves ontology terms from Superphy blazegraph
        and initializes internal dictionary object

        Args:
            graph (object): SupephyGraph object

        Returns:
            dictionary: Ontology uri -> list of categories

        """
        pass


    @abc.abstractmethod
    def individual(self, uri):
        """Checks if uri is part of ontology and
        returns instance of attribute class if found.

        Args:
            uri (string): uri in ontology

        Returns:
            object: Instance of attribute class matching uri

        """
        pass


class HostOntology(AttributeOntology):
    """Interface for Host ontology


    """
    
    def __init__(self, graph):
        """Constructor

        Retrieves ontology terms from DB

        Args:
            graph (object): SupephyGraph object

        """
        self._ontology = self._initialize_ontology(graph)

        self._attribute_class = 'HostIndividual'


    def attribute_class(self):
        """Attribute class represented by ontology

        """
        return self._attribute_class()


    def _initialize_ontology(self, graph):
        """Retrieves ontology terms from Superphy blazegraph
        and initializes internal dictionary object

        Args:
            graph (object): SupephyGraph object

        Returns:
            dictionary: Ontology uri -> list of categories

        """
        
        hosts = graph.query(
            """SELECT ?h ?c ?l
               WHERE {
                ?h a superphy:Host ;
                superphy:has_host_category ?c ;
                rdfs:label ?l .
               }""")
       
        hosts = strip_superphy_namespace(hosts)
        
        ontology = {}
        for h in hosts:
            uri = str(h[0])
            cat = str(h[1])
            lab = str(h[2])
            ontology[uri] = {
                '_uri': uri,
                '_category': cat,
                '_label': lab
            }

        return ontology


    def individual(self, uri):
        """If uri is part of ontology and
        returns instance of attribute class if found.

        Throws exception if none found

        Args:
            uri (string): uri in ontology

        Returns:
            object: Instance of attribute class matching uri

        """
        if uri in self._ontology:
            host = self._ontology[uri]
            constuctor = globals()[self._attribute_class]
            return constuctor(host)
        else:
            raise SuperphyMetaError("Unrecognized host uri: %s"%uri)


    def uri_list(self):
        """Get list of value uri's for superphy term

        Args:
            None

        Returns:
            List: list of strings

        """

        return list(self._ontology.keys())


class HostIndividual(object):
    """Represents instance of host attribute

        Tracks category and label associated with host

    """

    def __init__(self, *args, **kwargs):
        """Constructor

        Args:
            args (dictionary): Expects keys '_uri', '_label', '_category'
            kwargs (keywords): See above
        """
        for dictionary in args:
            for key in dictionary:
                setattr(self, key, dictionary[key])
        for key in kwargs:
            setattr(self, key, kwargs[key])


    def uri(self):
        return self._uri

    def category(self):
        return self._category

    def label(self):
        return self._label

    

class SourceOntology(AttributeOntology):
    """Interface for Source ontology


    """
    
    def __init__(self, graph):
        """Constructor

        Retrieves ontology terms from DB

        Args:
            graph (object): SupephyGraph object

        """
        self._ontology = self._initialize_ontology(graph)

        self._attribute_class = 'SourceIndividual'


    def attribute_class(self):
        """Attribute class represented by ontology

        """
        return self._attribute_class()


    def _initialize_ontology(self, graph):
        """Retrieves ontology terms from Superphy blazegraph
        and initializes internal dictionary object

        Args:
            graph (object): SupephyGraph object

        Returns:
            dictionary: Ontology uri -> list of categories

        """
        
        sources = graph.query(
            """SELECT ?h ?c ?l
               WHERE {
                ?h a superphy:isolation_from_source ;
                superphy:has_host_category ?c ;
                rdfs:label ?l .
               }""")
       
        sources = strip_superphy_namespace(sources)
        
        ontology = {}
        for s in sources:
            uri = str(s[0])
            cat = str(s[1])
            lab = str(s[2])
            if uri in ontology:
                # Add new category
                ontology[uri]['_category'].append(cat)
            else:
                # New source
                ontology[uri] = {
                    '_uri': uri,
                    '_category': [ cat ],
                    '_label': lab
                }

        return ontology


    def individual(self, uri):
        """If uri is part of ontology and
        returns instance of attribute class if found.

        Throws exception if none found

        Args:
            uri (string): uri in ontology

        Returns:
            object: Instance of attribute class matching uri

        """
        if uri in self._ontology:
            source = self._ontology[uri]
            constuctor = globals()[self._attribute_class]
            return constuctor(source)
        else:
            raise SuperphyMetaError("Unrecognized source uri: %s"%uri)


class SourceIndividual(object):
    """Represents instance of host attribute

        Tracks category and label associated with host

    """

    def __init__(self, *args, **kwargs):
        """Constructor

        Args:
            args (dictionary): Expects keys '_uri', '_label', '_category'
            kwargs (keywords): See above
        """
        for dictionary in args:
            for key in dictionary:
                setattr(self, key, dictionary[key])
        for key in kwargs:
            setattr(self, key, kwargs[key])


    def uri(self):
        return self._uri()

    def category(self):
        return self._category

    def label(self):
        return self._label

    def belongs(self, c):
        """Checks if category is in sources list of allowable categories

        Args:
            c(str): category name 

        Returns:
            bool

        """

        if c in self._category:
            return True
        else:
            return False
    

class GenomeRecord(object):
    """Stores genome metadata

    Metadata attributes are 

    Attributes:
        _host_list(List): List of HostIndividual objects assigned to genome
        _source_list(List): List of SourceIndividual objects assigned to genome
        _uniquename(str): genome accession string

    """

    def __init__(self, uniquename):
        """Constructor

        Initializes all superphy attributes to None.  As attributes
        are parsed they will be recorded using the setattr function
       
        Args:
            uniquename(string): URI for genome used in superphy

        """
        self._host_list = []
        self._source_list = []
        self._uniquename = uniquename
        

    def is_valid(self):
        """Checks for conflicts in metadata attributes

        Args:
            None

        Returns:
            tuple containing:
                bool: True if no conflicts found
                str: Further information if bool == False, else None

       
        """

        # There can only be one host
        # Host list is a set of unique HostIndividual objects
        if len(self._host_list) > 1:
            m = [ h.uri() for h in self._host_list ]
            return (False, 'Multiple hosts detected: ' + pformat(m))

        # There can only be one source
        if len(self._source_list) > 1:
            m = [ h.uri() for h in self._source_list ]
            return (False, 'Multiple sources detected: ' + pformat(m))

        # Check for category conflicts
        current_category = None
        if self._host_list:
            current_category = self._host_list[0].category()

        category_based_data = itertools.chain(self._source_list)
        for d in category_based_data:
            if not current_category:
                current_category = d.category()

            else:
                if not d.belongs(current_category):
                    # No overlap between categories
                    return (False, 'Category conflict! Data individual %s is lacking category %s'%(d.uri(), pformat(current_category)))


        return (True, None)


    def output(self):
        """Return genome metadata as dictionary

        Returns:
            dict: Superphy metadata keys and values


        """

        genome_dict = {}

        if self._host_list:
            genome_dict['isolation_host'] = self._host_list[0].uri()

        if self._source_list: 
             genome_dict['isolation_source'] = self._host_source[0].uri()


        return genome_dict


    def uniquename(self):
        """Return uniquename/accession for genome object

        """

        return self._uniquename
        

    def host(self, host_individual):
        """Add host to host list if its a new unique host

        """
        
        # Check if host already in list
        exists = False
        for i in self._host_list:
            if host_individual.uri == i.uri:
                exists = True
                break


        if not exists:
            self._host_list.append(host_individual)




