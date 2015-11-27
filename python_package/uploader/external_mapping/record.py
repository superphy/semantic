#!/usr/bin/env python

"""Container objects for temporarily storing metadata.

Classes store metadata and verify correctness and possible conflicts in metadata

Example:
    Usage:

        $ python example_google.py


"""

import abc
from rdflib.namespace import RDF
from superphy.endpoint import superphy_namespace, strip_superphy_namespace
from pprint import pprint


__author__ = "Matt Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Matt Whiteside"
__email__ = "matthew.whiteside@canada.ca"



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


    @property
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
            ontology[h[0]] = {
                '_uri': h[0],
                '_category': h[1],
                '_label': h[2]
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
            return constuctor(uri, host)
        else:
            raise SuperphyMetaError("Unrecognized host uri: %s"%uri)


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
            setattr(self, key, kwargs[key])):


    @property
    def uri(self):
        return self._uri()

    @property
    def category(self):
        return self._category

    @property
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


    @property
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
            if s[0] in ontology:
                # Add new category
                ontology[s[0]]['_category'].append(s[1])
            else:
                # New source
                ontology[s[0]] = {
                    '_uri': s[0],
                    '_category': [ s[1] ],
                    '_label': s[2]
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
            return constuctor(uri, source)
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
            setattr(self, key, kwargs[key])):


    @property
    def uri(self):
        return self._uri()

    @property
    def category(self):
        return self._category

    @property
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
        host(HostIndividual): Host assigned to genome
        source(SourceIndividual): Source assigned to genome

    """

    def __init__(self):
        """Constructor

        Initializes all superphy attributes to None.  As attributes
        are parsed they will be recorded using the setattr function
       
        Args:
            None

        """
        self.host = None
        self.source = None

    def is_valid():
        """Checks for conflicts in metadata attributes


        """

        # Check host

        # Check source
        pass


