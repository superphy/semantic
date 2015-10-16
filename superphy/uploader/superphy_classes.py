__author__ = 'Stephen Kan'

from rdflib import Graph, Namespace, Literal, XSD
import superphy_sparql

"""
This module converts inputted data into RDF triples in accordance to the Superphy ontology

This module uses the rdflib to handle RDF triples in a Graph instance in preparation for conversion
into a turtle or other RDF format.

The format to add triples would be:

    g.add( (subject, predicate, object) )

where each of the three terms ( can be constructed by using:

    namespace.objectname
    namespace[vars_containing_object_name]
    Literal(vars_containing_literal, datatype=XSD.typeofdata)

Examples:
    superphy = Namespace("https://github.com/superphy#")

    organism = "ecoli"
    tid = "562"

    g.add( (superphy[organism], superphy.has_taxonomy_id, Literal(str(tid), datatype=XSD.string) )

If you want to have an RDF object with many optional arguments, Genome is a good example.

Classes:
    NamedIndividual: superclass of all RDF Objects
    User: a Superphy user
    Organism: an organism
    Host: a host for another organism
    Microbe: a microbe
    Attribute: attributes of objects in the ontology (abstract class)
    HostCategory: host category for an object
    IsolationAttribute: details associated with isolation of a genome (abstract class)
    FromHost: the host from which the organism was extracted from
    FromSource: the biological source of the sample
    IsolationSyndrome: the syndromes associated with the sample
    Serotype: the serotype of the genome (abstract class)
    Otype: the id of the O-antigen associated with the sample
    Htype: the id of the H-antigen associated with the sample
    Genome: a genome sequenced from a sample
    PendingGenome: a genome that has not finished processing
    CompletedGenome: a genome that has finished processing

Methods:
    generate_output: exports uploaded RDFs to a turtle file
"""

# initialize a Graph
g = Graph()

# setting up namespaces for use
n = Namespace("https://github.com/superphy#")
owl = Namespace("http://www.w3.org/2002/07/owl#")
rdf = Namespace("http://www.w3.org/1999/02/22-rdf-syntax-ns#")
xml = Namespace("http://www.w3.org/XML/1998/namespace")
xsd = Namespace("http://www.w3.org/2001/XMLSchema#")
rdfs = Namespace("http://www.w3.org/2000/01/rdf-schema#")
gfvo = Namespace("http://www.biointerchange.org/gfvo#")

g.bind("", "https://github.com/superphy#")


class NamedIndividual(object):
    """
    The superclass for all RDF Objects (excluding Blank Nodes as they are not named)
    """

    def __init__(self, name):
        """
        Create a NamedIndividual with associated metadata
        """

        self.name = name

    def rdf(self):
        """
        Convert NamedIndividual metadata into RDF
        """

        g.add( (n[self.name], rdf.type, owl.NamedIndividual) )


class User(NamedIndividual):
    """
    This class is created when a new user is registered.
    """

    def __init__(self, email):
        """
        Create a new User with associated metadata

        Args:
            email (str): is both the individual name, and the only field literal
        """

        super(User,self).__init__(email)
        self.email = email

    def rdf(self):
        """
        Convert User metadata into RDF
        """

        super(User, self).rdf()
        g.add( (n[self.name], rdf.type, n.User) )
        g.add( (n[self.name], n.email, Literal(str(self.email), datatype = XSD.string)))


class Organism(NamedIndividual):
    """
    A organism
    """

    def __init__(self, name, label, scientific_name, common_name, taxonomy_id=None):
        """
        Create an Organism with associated metadata

        Args:
            name (str): the name of the Organism
            label (str): a label used by Meta::Miner to refer to the Organism
            scientific_name (str): the scientific name (genus species) of the Organism
            common_name (str): the common name of the Organism
            taxonomy_id: the taxonomy id for the Organism as given by NCBI
        """

        super(Organism, self).__init__(name)
        self.label = label
        self.scientific_name = scientific_name
        self.common_name = common_name
        self.taxonomy_id = taxonomy_id

    def rdf(self):
        """
        Convert Organism metadata into RDF
        """

        super(Organism, self).rdf()
        g.add( (n[self.name], rdf.type, n.Organism) )
        g.add( (n[self.name], rdfs.label, Literal(str(self.label), datatype=XSD.string)) )
        g.add( (n[self.name], n.scientific_name, Literal(str(self.scientific_name), datatype=XSD.string)) )
        g.add( (n[self.name], n.common_name, Literal(str(self.common_name), datatype=XSD.string)) )
        g.add( (n[self.name], n.has_taxonomy_id, Literal(str(self.taxonomy_id), datatype=XSD.string)) )


class Host(Organism):
    """
    A host for another organism
    """

    def __init__(self, name, label, scientific_name, common_name, host_category, taxonomy_id=None):
        """
        Create a Host instance with associated metadata

        Args:
            - refer to Organism
            host_category (str): the host category that Host belongs to
        """

        super(Host, self).__init__(name, label, scientific_name, common_name, taxonomy_id)
        self.host_category = host_category

    def rdf(self):
        """
        Convert Host metadta into RDF
        """

        super(Host, self).rdf()
        g.add( (n[self.name], rdf.type, n.Host) )
        g.add( (n[self.name], n.has_host_category, n[self.host_category]) )
        g.add( (n[self.host_category], n.is_host_category_of, n[self.name]) )
        FromHost(self.name, self.host_category).rdf()


class Microbe(Organism):
    """
    A microbe
    """

    def rdf(self):
        """
        Convert Microbe metadata into RDF
        """

        super(Microbe, self).rdf()
        g.add( (n[self.name], rdf.type, n.Microbe) )


class Attribute(NamedIndividual):
    """
    An attribute class used to describe objects in the ontology
    """

    def rdf(self):
        """
        Convert Attribute metadata into RDF
        """

        super(Attribute, self).rdf()
        g.add( (n[self.name], rdf.type, n.Attribute) )


class HostCategory(Attribute):
    """
    An host category attribute with associated metadata

    An attribute describing valid host categories for other attributes
    """

    def __init__(self, name, label):
        """
        Create a HostCategory with associated metadata

        Args:
            name (str): name of the HostCategory
            label (str): a label used by Meta::Miner to refer to the HostCategory
        """

        super(HostCategory, self).__init__(name)
        self.label = label

    def rdf(self):
        """
        Convert HostCategory metadata into RDF
        """

        super(Attribute, self).rdf()
        g.add( (n[self.name], rdf.type, n.host_category) )
        g.add( (n[self.name], rdfs.label, Literal(str(self.label), datatype=XSD.string)) )


class IsolationAttribute(Attribute):
    """
    An isolation attribute with associated metadata

    Attributes that are associated or describes how a genome sample was isolated
    """

    def rdf(self):
        """
        Convert IsolationAttribute metadata into RDF
        """

        super(IsolationAttribute, self).rdf()
        g.add( (n[self.name], rdf.type, n.isolation_attribute) )


class FromHost(IsolationAttribute):
    """
    A "from host" attribute with associated metadata
    """

    def __init__(self, host, host_category):
        """
        Create a FromHost with associated metadata

        Args:
            host (str): the host that FromHost references
            host_category (str): the host category that FromHost belongs to
        """

        self.name = "from_" + host
        super(FromHost, self).__init__(self.name)
        self.host = host
        self.host_category = host_category

    def rdf(self):
        """
        Convert FromHost metadata into RDF
        """

        super(FromHost, self).rdf()
        g.add( (n[self.name], rdf.type, n.isolation_from_host) )
        g.add( (n[self.name], n.has_object, n[self.host]) )
        g.add( (n[self.host], n.is_object_of, n[self.name]))
        g.add( (n[self.name], n.has_host_category, n[self.host_category]) )
        g.add( (n[self.host_category], n.is_host_category_of, n[self.name]) )


class FromSource(IsolationAttribute):
    """
    A biological source with associated metadata
    """

    def __init__(self, name, label, host_category):
        """
        Create a FromSource with associated metadata

        Args:
            name (str): the name of the FromSource
            label (str): a label used by Meta::Miner to refer to the FromSource
            host_category (str): the host category of the FromSource
        """

        super(FromSource, self).__init__(name)
        self.label = label
        self.host_category = host_category

    def rdf(self):
        """
        Convert FromSource metadata into RDF
        """

        super(FromSource, self).rdf()
        g.add( (n[self.name], rdf.type, n.isolation_from_source) )
        g.add( (n[self.name], rdfs.label, Literal(str(self.label), datatype=XSD.string)) )
        g.add( (n[self.name], n.has_host_category, n[self.host_category]) )
        g.add( (n[self.host_category], n.is_host_category_of, n[self.name]) )


class IsolationSyndrome(IsolationAttribute):
    """
    A syndrome with associated metadata
    """

    def __init__(self, name, label, host_category):
        """
        Create an IsolationSyndrome with associated metadata

        Args:
            name (str): the name of the IsolationSyndrome
            label (str): a label used by Meta::Miner to refer to the IsolationSyndrome
            host_category (str): the host category of the IsolationSyndrome
        """
        super(IsolationSyndrome, self).__init__(name)
        self.label = label
        self.host_category = host_category

    def rdf(self):
        """
        Convert IsolationSyndrome metadata into RDF
        """

        super(IsolationSyndrome, self).rdf()
        g.add( (n[self.name], rdf.type, n.isolation_syndrome) )
        g.add( (n[self.name], rdfs.label, Literal(str(self.label), datatype=XSD.string)) )
        g.add( (n[self.name], n.has_host_category, n[self.host_category]) )
        g.add( (n[self.host_category], n.is_host_category_of, n[self.name]) )


class Serotype(Attribute):
    """
    A serotype with associated metadata
    """

    def rdf(self):
        """
        Convert Serotype metadata into RDF
        """

        super(Serotype, self).rdf()
        g.add( (n[self.name], rdf.type, n.serotype) )


class Otype(Serotype):
    """
    A O-type antigen with associated metadata
    """

    def __init__(self, id):
        """
        Create a Otype with associated metadata

        Args:
            id (str): the id of the O antigen
        """

        self.name = "O" + str(id)
        super(Otype, self).__init__(self.name)

    def rdf(self):
        """
        Convert Otype metadata into RDF
        """
        super(Otype, self).rdf()
        g.add( (n[self.name], rdf.type, n.Otype) )


class Htype(Serotype):
    """
    A H-type antigen with associated metadata
    """

    def __init__(self, id):
        """
        Create a Htype with associated metadata

        Args:
            id (str): the id of the H antigen
        """

        self.name = "H" + str(id)
        super(Htype, self).__init__(self.name)

    def rdf(self):
        """
        Convert Htype metadata into RDF
        """

        super(Htype, self).rdf()
        g.add( (n[self.name], rdf.type, n.Htype) )


class Genome(NamedIndividual):
    """
    A genome with associated metadata

    """

    def __init__(self, name, **kwargs):
        """
        Create a Genome with associated metadata

        Since most of the arguments for genome are optional, **kwargs is used to pass them in, after
        filtering them via SearchParam's keywords to prevent spurious and potentially malicious keys
        from being kept.

        Attributes:
            searchparam: list of keywords to filter kwargs by

        Args:
            name (str): name of the Genome
            kwargs (dict): optional arguments for Genome, that will be filtered for
                           spurious entries

        """

        searchparam = ["date", "location", "accession", "bioproject", "biosample", "strain", "organism", "host",
                       "source", "syndrome", "Htype", "Otype", "User"]

        super(Genome, self).__init__(name)
        self.kwargs = {key:value for key, value in kwargs.items() if key in searchparam}

    def rdf(self):
        """
        Convert all Genome metadata to RDF

        To ensure that H and O types are assigned Unknowns, if kwargs does not include at least
        one of them, it would call the appropriate RDF method.

        As methods are called based on their keys, there is some coupling with minerJSON.py to
        ensure that the keys are properly named.
        """

        super(Genome, self).rdf()

        if "Htype" not in self.kwargs:
            self.Htype()
        if "Otype" not in self.kwargs:
            self.Otype()

        for key, value in self.kwargs.iteritems():
            getattr(self, key)(value)

    def date(self, date):
        """
        Convert all date entries into RDF

        Args:
            date: a collection of sampling dates for the Genome in the XSD date format (YYYY-MM-DD)
        """

        for item in date:
            literal = Literal(item, datatype=XSD.date)
            g.add( (n[self.name], n.has_isolation_date, literal) )

    def location(self, location):
        """
        Convert all location entries into RDF

        Args:
            location: a collection of sampling locations for the Genome
        """

        for item in location:
            literal = Literal(item, datatype=XSD.string)
            g.add( (n[self.name], n.has_geographic_location, literal) )

    def accession(self, accession):
        """
        Convert all NCBI Genbank/Nuccore accession numbers into RDF

        Args:
            accession: a collection of Nuccore accession ids associated with the Genome
        """

        for item in accession:
            literal = Literal(item, datatype=XSD.string)
            g.add( (n[self.name], n.has_accession, literal) )

    def bioproject(self, bioproject):
        """
        Convert all BioProject ids into RDF

        Args:
            bioproject: a collection of BioProject ids associated with the Genome
        """

        for item in bioproject:
            literal = Literal(item, datatype=XSD.string)
            g.add( (n[self.name], n.has_bioproject, literal) )

    def biosample(self, biosample):
        """
        Convert all BioSample ids into RDF

        Args:
            biosample: a collection of BioSample ids associated with the Genome
        """

        for item in biosample:
            literal = Literal(item, datatype=XSD.string)
            g.add( (n[self.name], n.has_biosample, literal) )

    def strain(self, strain):
        """
        Convert all strain names into RDF

        Args:
            strain: a collection of strain names associated with the Genome
        """

        for item in strain:
            literal = Literal(item, datatype=XSD.string)
            g.add( (n[self.name], n.has_strain, literal) )

    def organism(self, organism):
        """
        Convert organism into RDF

        Args:
            organism (str): name of the organism of the Genome
        """

        g.add( (n[self.name], n.is_genome_of, n[organism]))
        g.add( (n[organism], n.has_genome, n[self.name]))

    def host(self, from_host):
        """
        Convert all host data into RDF

        Args:
            from_host: a collection of hosts that the Genome has been sampled from
        """

        for item in from_host:
            node = superphy_sparql.find_from_host(item).split("#", 1)[1]
            g.add( (n[self.name], n.has_isolation_attribute, n[node]) )
            g.add( (n[node], n.is_isolation_attribute_of, n[self.name]) )

    def source(self, from_source):
        """
        Convert all source data into RDF

        Args:
            from_source: a collection of biological sources that the Genome has been sampled from
        """

        for item in from_source:
            node = superphy_sparql.find_source(item).split("#", 1)[1]
            g.add( (n[self.name], n.has_isolation_attribute, n[node]) )
            g.add( (n[node], n.is_isolation_attribute_of, n[self.name]) )

    def syndrome(self, syndrome):
        """
        Convert all syndrome data into RDF

        Args:
            syndrome: a collection of syndromes associated with the Genome
        """

        for item in syndrome:
            node = superphy_sparql.find_syndrome(item).split("#", 1)[1]
            g.add( (n[self.name], n.has_isolation_attribute, n[node]) )
            g.add( (n[node], n.is_isolation_attribute_of, n[self.name]) )

    def Htype(self, Htype=None):
        """
        Convert H serotype into RDF

        Args:
            Htype: the id of the H-antigen associated with the Genome, or None if not provided
        """

        if Htype is None:
            g.add((n[self.name], n.has_Htype, n.HUnknown))
            g.add((n.HUnknown, n.is_Htype_of, n[self.name]))
        else:
            g.add((n[self.name], n.has_Htype, n["H" + str(Htype)]))
            g.add((n["H" + str(Htype)], n.is_Htype_of, n[self.name]))

    def Otype(self, Otype=None):
        """
        Convert O serotype into RDF

        Args:
            Otype: the id of the O-antigen associated with the Genome, or None if not provided
        """

        if Otype is None:
            g.add((n[self.name], n.has_Otype, n.OUnknown))
            g.add((n.OUnknown, n.is_Otype_of, n[self.name]))
        else:
            g.add((n[self.name], n.has_Otype, n["O" + str(Otype)]))
            g.add((n["O" + str(Otype)], n.is_Otype_of, n[self.name]))

    def User(self, User):
        """
        Converts User id into RDF

        Args:
            User: the id of the user who uploaded the Genome, restricting permissions to him
        """

        g.add( (n[self.name], n.is_owned_by, n[User]) )
        g.add( (n[User], n.owns, n[self.name]) )


class PendingGenome(Genome):
    """
    A genome that has not completed sequence analysis.

    Attributes:
        Same as Genome
    """

    def rdf(self):
        """
        Convert all PendingGenome variables to RDF and tags the Genome as a PendingGenome
        """

        super(PendingGenome, self).rdf()
        g.add( (n[self.name], rdf.type, n.pending_genome) )


class CompletedGenome(Genome):
    """
    A genome with completed sequence analysis

    Attributes:
        Same as Genome
    """

    def rdf(self):
        """
        Convert all CompletedGenome variables to RDF and tags the Genome as a CompletedGenome

        """

        super(CompletedGenome, self).rdf()
        g.add( (n[self.name], rdf.type, n.completed_genome) )


def generate_output(destination):
    """
    Export RDF Graph data to a turtle file at the given destination

    Args:
        destination: an internal filepath relative to the  __init__.py file this module belongs to
    """

    g.serialize(destination=destination, format="turtle")
    g.remove( (None, None, None) )



""" =================================================== TESTING =================================================== """

"""
item = Host("hsapiens", "Homo sapiens (human)", "Homo sapiens", "human", "human")
item.rdf()
item = Microbe("ecoli", "Escherichia coli (E. coli)", "Escherichia coli", "E. coli")
item.rdf()

list(Htype(num).rdf()  for num in xrange(1,56))
list(Otype(num).rdf() for num in xrange(1,187))

genome = PendingGenome("NC_455763", **{"date":{"2010-10-09"}, "location":{"Canada"}, "Htype":"7", "from_host":{"Homo sapiens (human)"}, "organism": "ecoli", "from_source":{"Blood"}})
print genome.kwargs
genome.rdf()

kwargs = {'from_source': set([]), 'accession': set(['CP001855']), 'strain': set([u'NRG 857C']), 'Otype': u'83', 'bioproject': set(['41221']), 'date': set([]), 'syndrome': set([]), 'from_host': set([]), 'location': set([]), 'Htype': u'1', 'organism': 'ecoli', 'biosample': set(['2603727'])}
PendingGenome("CP001855", **kwargs).rdf()

g.serialize(destination="outputs/test.txt", format="turtle")

"""