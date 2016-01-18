#!/usr/bin/env python

"""Container objects for temporarily storing metadata.

Classes store metadata and verify correctness and possible conflicts in metadata

Example:
    Usage:

        $ python example_google.py


"""

import abc
import itertools
import json
from pprint import pformat
from rdflib import Literal
from geopy.geocoders import GoogleV3
from superphy.endpoint import strip_superphy_namespace
from superphy.config import parser
from superphy.uploader.namespaces import SUPERPHY, RDF

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
    pass


class AttributeOntology(object):
    """Abstract class that defines interface for storing sections of the Superphy metadata ontology


    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, store):
        """Constructor

        Args:
            store (object): Supephystore object

        """
        self._ontology = self._initialize_ontology(store)
        self._attribute_class = 'Attribute'

        pass


    @abc.abstractproperty
    def attribute_class(self):
        """Attribute class represented by ontology

        """
        return 'Should never get here'


    @abc.abstractmethod
    def _initialize_ontology(self, store):
        """Retrieves ontology terms from Superphy blazegraph
        and initializes internal dictionary object

        Args:
            store (object): Supephystore object

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
    
    def __init__(self, store, logger):
        """Constructor

        Retrieves ontology terms from DB

        Args:
            store (object): SupephyStore object

        """

        self.logger = logger

        self._ontology = self._initialize_ontology(store)

        self._attribute_class = 'HostIndividual'


    def attribute_class(self):
        """Attribute class represented by ontology

        """
        return self._attribute_class()


    def _initialize_ontology(self, store):
        """Retrieves ontology terms from Superphy blazegraph
        and initializes internal dictionary object

        Args:
            store (object): Supephystore object

        Returns:
            dictionary: Ontology uri -> list of categories

        """
        
        hosts = store.query(
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
    
    def __init__(self, store, logger):
        """Constructor

        Retrieves ontology terms from DB

        Args:
            store (object): Supephystore object

        """
        self.logger = logger

        self._ontology = self._initialize_ontology(store)

        self._attribute_class = 'SourceIndividual'


    def attribute_class(self):
        """Attribute class represented by ontology

        """
        return self._attribute_class()


    def _initialize_ontology(self, store):
        """Retrieves ontology terms from Superphy blazegraph
        and initializes internal dictionary object

        Args:
            store (object): Supephystore object

        Returns:
            dictionary: Ontology uri -> list of categories

        """
        
        sources = store.query(
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


    def uri_list(self):
        """Get list of value uri's for superphy term

        Args:
            None

        Returns:
            List: list of strings

        """

        return list(self._ontology.keys())


class SourceIndividual(object):
    """Represents instance of source attribute

        Tracks category and label associated with source

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
    

class SyndromeOntology(AttributeOntology):
    """Interface for Syndrome ontology


    """
    
    def __init__(self, store, logger):
        """Constructor

        Retrieves ontology terms from DB

        Args:
            store (object): Supephystore object

        """
        self.logger = logger

        self._ontology = self._initialize_ontology(store)

        self._attribute_class = 'SyndromeIndividual'


    def attribute_class(self):
        """Attribute class represented by ontology

        """
        return self._attribute_class()


    def _initialize_ontology(self, store):
        """Retrieves ontology terms from Superphy blazegraph
        and initializes internal dictionary object

        Args:
            store (object): Supephystore object

        Returns:
            dictionary: Ontology uri -> list of categories

        """
        
        syndrome = store.query(
            """SELECT ?h ?c ?l
               WHERE {
                ?h a superphy:isolation_syndrome ;
                superphy:has_host_category ?c ;
                rdfs:label ?l .
               }""")
       
        syndrome = strip_superphy_namespace(syndrome)
        
        ontology = {}
        for s in syndrome:
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
            syndrome = self._ontology[uri]
            constuctor = globals()[self._attribute_class]
            return constuctor(syndrome)
        else:
            raise SuperphyMetaError("Unrecognized syndrome uri: %s"%uri)


    def uri_list(self):
        """Get list of value uri's for superphy term

        Args:
            None

        Returns:
            List: list of strings

        """

        return list(self._ontology.keys())


class SyndromeIndividual(object):
    """Represents instance of syndrome attribute

        Tracks category and label associated with syndrome

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

    def belongs(self, c):
        """Checks if category is in syndrome list of allowable categories

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
        self._location_list = []
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

    def source(self, source_individual):
        """Add isolation_source to source list if its a new unique source

        """
        
        # Check if host already in list
        exists = False
        for i in self._source_list:
            if source_individual.uri == i.uri:
                exists = True
                break

        if not exists:
            self._source_list.append(source_individual)

    def location(self, location_individual):
        """Add isolation_location to location list if its a new unique source

        """
        
        # Check if host already in list
        exists = False
        for i in self._location_list:
            if location_individual.uri == i.uri:
                exists = True
                break

        if not exists:
            self._location_list.append(location_individual)




class LocationOntology(AttributeOntology):
    """Interface for Locations


    """
    
    def __init__(self, store, logger):
        """Constructor

        Interface to locations cache in DB

        Args:
            store (object): Supephystore object

        """
        self.logger = logger

        self._ontology = self._initialize_ontology(store)

        self._attribute_class = 'LocationIndividual'


    def attribute_class(self):
        """Attribute class represented by ontology

        """
        return self._attribute_class()


    def _initialize_ontology(self, store):
        """Prepares queries and initialize geocoder

        Args:
            store (object): Supephystore object

        Returns:
            None

        """

        self.queries = {
            'uri_lookup':
                """ASK { 
                    ?uri a superphy:geographic_location
                   }
                """,

            'search_term_lookup':
                """SELECT ?l
                   WHERE { 
                    ?l a superphy:geographic_location .
                    ?l superphy:matches_location_query ?search_term
                   }
                """
        }

        self.store = store

        # Initialize geocoder
        apikey = parser.read()['geocoding_apikey']
        if not(apikey):
            raise Exception("Undefined config variable: geocoding_apikey")
        self.geolocater = GoogleV3(api_key=apikey)


        return None


    def individual(self, search_term):
        """If search_term is part of ontology,
        returns instance of attribute class.

        Throws exception if none found

        Args:
            search_term (string): search_term string in ontology

        Returns:
            object: Instance of attribute class matching uri

        """

        g = self.store
        result = g.query(self.queries['search_term_lookup'], initBindings={'search_term': Literal(search_term)})
       
        if bool(result):
            
            constuctor = globals()[self._attribute_class]
            for uri in result:
                self.logger.debug("Adding LocationIndividual %s"%uri)
                return constuctor({'_uri': uri})
            
        else:
            raise SuperphyMetaError("Unrecognized location search_term: %s"%search_term)


    def _location_query(self, search_term): 
        """Returns location with matching search_term property

        Args:
            search_term[string]: location string

        Returns:
            rdflib.query.Result

        """

        g = self.store

        # Create variables
        st = Literal(search_term)
       
        result = g.query(self.queries['search_term_lookup'], initBindings={'search_term': st})

        return result



    def has_location(self, search_term):
        """Checks if location exists in db with matching search term

        Args:
            search_term[string]: location string

        Returns:
            bool

        """

        result = self._location_query(search_term)
        
        return bool(result)


    def get_location(self, search_term):
        """Adds new location to cache (if it does not exist already)

        Args:
            search_term[string]: location string

        Returns:
            string: new/existing location URI

        Raises:
            SuperphyMetaError: if location is invalid/geocoding failed

        """

        result = self._location_query(search_term)

        self.logger.debug("Found location: %s"%bool(result))

        if result:
            # Found term in DB matching search_term
            locations = strip_superphy_namespace(result)

            return locations[0]

        else:
            # No matching entry in cache with search_term

            # Use geocoder to translate/validate the search_term
            location = self.geolocater.geocode(search_term, exactly_one=True)

            self.logger.debug("Geocoded location: %s"%location)

            if not (location):
                raise SuperphyMetaError("Invalid location/geocoding failed for %s"%search_term)

            googleResponse = location.raw
            
            # Retreive place ID (only one result requested)
            place_id = googleResponse[u'place_id'].decode('utf-8','ignore')
            assert place_id, "place_id not found in google Maps API response: %r"%place_id

            uristr = "googlemaps_place_id_"+place_id
            uri = SUPERPHY[uristr] # use place_id as uri
            self.logger.debug("New geographic_location URI being added: %s"%uri.n3())

            result = self.store.query(self.queries['uri_lookup'], initBindings={'uri': uri})
            if result:
                # Found location object in DB with that uri
                self.logger.debug("Linking new location search term %s"%search_term)

                # Link this search term to formatted address
                self._add_search_term(uri, search_term)

                return uristr

            else:
                # Need to add new location object
                self.logger.debug("Adding new geographic_location URI %s"%uri.n3())

                self._add_db_location(location, uri, search_term)

                return uristr


    def _add_db_location(self, location, uri, search_term):
        """Insert new location object into DB

        Args:
            location[geopy.location.Location]: location object
            uri[string|rdflib.term.URIRef]: URI to use in DB
            search_term[string]: Search term string linked to geocoded address

        Returns:
            SparqlWrapper.Wrapper.QueryResult

        """

        # Create variables
        if isinstance(uri,basestring):
            uri = SUPERPHY[uri]

        googleResponse = location.raw
        jsonResponse = json.dumps(googleResponse)

        po = [
            (RDF.type, SUPERPHY.geographic_location),
            (SUPERPHY.matches_location_query, Literal(search_term)),
            (SUPERPHY.has_formatted_address, Literal(location.address)),
            (SUPERPHY.has_geocoding_result, Literal(jsonResponse))
        ]
      
       
        # Run query
        try:

            # Add insert stmts
            for p, o in po:
                self.store.add((uri, p, o))

            # Send request to server
            r = self.store.commit()
            self.logger.debug("Insert call returned response %s"%r.info())
            self.logger.debug("Insert url %s"%r.geturl())
            self.logger.debug("Insert results %s"%r.convert())
            self.logger.debug("QUERY:\n%s"%self.store.queryString)

            return r

        except Exception, e:
            self.logger.debug("Insert call raised exception %s"%e)
            raise e



    def _add_search_term(self, uri, search_term):
        """Link search term to location object into DB

        Args:
            location[geopy.location.Location]: location object
            uri[string|rdflib.term.URIRef]: URI to use in DB
            search_term[string]: Search term string linked to geocoded address

        Returns:
            None on success, raises exception on failure/error

        """

        """

            'insert_search_term':
                    INSERT DATA { 
                        ?location superphy:matches_location_query ?search_term .
                       }
                    
        
        """

        # Create variables
        if isinstance(uri,basestring):
            uri = URIRef(uri)

        po = [
            (URIRef('superphy:matches_location_query'), Literal(search_term)),
        ]
      
       
        # Run query
        try:

            # Add insert stmts
            for p, o in po:
                self.store.add((uri, p, o))

            # Send request to server
            r = self.store.store.commit()
            self.logger.debug("Insert call returned response %s"%r.info())

            return r

        except Exception, e:
            self.logger.debug("Insert call raised exception %s"%e)
            raise e



class LocationIndividual(object):
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


    

