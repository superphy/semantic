__author__ = 'Stephen Kan'

"""
NOTE: when creating objects that require the creation of reciprocal objects, always pair the g.add functions together
to maintain consistency. Even if there is a create_anotherfn, add both relationship triples at the same level of object creation

TODO: strip spaces and replace with underscores for URLnames (possibly via URLencode)
"""

from rdflib import Graph, Namespace, Literal, XSD
import superphy_SPARQL

# initialize a RDF Graph
g = Graph()

# setting up namespaces for use
n = Namespace("https://github.com/superphy#")
owl = Namespace("http://www.w3.org/2002/07/owl#")
rdf = Namespace("http://www.w3.org/1999/02/22-rdf-syntax-ns#")
xml = Namespace("http://www.w3.org/XML/1998/namespace")
xsd = Namespace("http://www.w3.org/2001/XMLSchema#")
rdfs = Namespace("http://www.w3.org/2000/01/rdf-schema#")
gfvo = Namespace("http://www.biointerchange.org/gfvo#")

g.bind("","https://github.com/superphy#")
g.bind("owl", "http://www.w3.org/2002/07/owl#")


def create_class(name):
    g.add( (n[name], rdf.type, owl.NamedIndividual) )


def create_organism(name, label, sci_name, com_name, taxonomy_id=None):
    create_class(name)
    g.add( (n[name], rdfs.label, Literal(str(label), datatype = XSD.string)) )
    g.add( (n[name], n.scientific_name, Literal(str(sci_name), datatype = XSD.string)) )
    g.add( (n[name], n.common_name, Literal(str(com_name), datatype = XSD.string)) )
    g.add( (n[name], n.has_taxonomy_id, Literal(str(taxonomy_id), datatype = XSD.string)) )


# assume: host_category is created
# TODO: add in a check for that
def create_host(name, host_category, label, sci_name, com_name, taxonomy_id=None):
    create_organism(name, label, sci_name, com_name, taxonomy_id)
    g.add( (n[name], rdf.type, n.Host) )
    g.add( (n[name], n.has_host_category, n[host_category]) )
    g.add( (n[host_category], n.is_host_category_of, n[name]) )
    create_from_host(name, host_category)

def create_from_host(host, host_category):
    name = "from_" + host
    create_class(name)
    g.add( (n[name], rdf.type, n.isolation_from_host) )
    g.add( (n[name], n.has_object, n[host]) )
    g.add( (n[host], n.is_object_of, n[name]))
    g.add( (n[name], n.has_host_category, n[host_category]) )
    g.add( (n[host_category], n.is_host_category_of, n[name]) )
    return name


def create_microbe(name, label, sci_name, com_name, taxonomy_id=None):
    create_organism(name, label, sci_name, com_name, taxonomy_id)
    g.add( (n[name], rdf.type, n.Microbe) )


# assumption: all other objects must have already been created; the genome metadata is always added last
# remember, date must be in dateTime format (YYYY-MM-DDTHH:MM:SS)
# TODO: refactor code to be less clunky (so much duplication!) and have more checks for created items
# TODO: include environment as a genome source (under organism or under source???)
def create_genome(name, date = None, location = None, accession = None, bioproject = None, biosample = None, strain = None,
                  organism = None, from_host = None, from_source = None, syndrome = None, Htype = None, Otype = None, User = None):
    create_class(name)
    g.add( (n[name], rdf.type, gfvo.Genome) )

    if date is not None:
        for d in date:
            g.add( (n[name], n.has_isolation_date, Literal(d, datatype = XSD.date)) )

    if location is not None:
        for l in location:
            g.add ( (n[name], n.has_geographic_location, Literal(l, datatype = XSD.string)) )

    if accession is not None:
        for a in accession:
            g.add( (n[name], n.has_accession, Literal(a, datatype = XSD.string)) )

    if bioproject is not None:
        for b in bioproject:
            g.add( (n[name], n.has_bioproject, Literal(b, datatype = XSD.string)) )

    if biosample is not None:
        for b in biosample:
            g.add( (n[name], n.has_biosample, Literal(b, datatype = XSD.string)) )

    if strain is not None:
        for s in strain:
            g.add( (n[name], n.has_strain, Literal(s, datatype = XSD.string)) )

    if organism is not None:
        g.add( (n[name], n.is_genome_of, n[organism]) )
        g.add( (n[organism], n.has_genome, n[name]))

    if from_host is not None:
        for h in from_host:
            node = superphy_SPARQL.find_from_host(h).split("#", 1)[1]
            g.add( (n[name], n.has_isolation_attribute, n[node]) )
            g.add( (n[node], n.is_isolation_attribute_of, n[name]) )

    if from_source is not None:
        for source in from_source:
            node = superphy_SPARQL.find_source(source).split("#", 1)[1]
            g.add( (n[name], n.has_isolation_attribute, n[node]) )
            g.add( (n[from_source], n.is_isolation_attribute_of, n[node]) )

    if syndrome is not None:
        for synd in syndrome:
            node = superphy_SPARQL.find_syndrome(synd).split("#", 1)[1]
            g.add( (n[name], n.has_isolation_attribute, n[node]) )
            g.add( (n[syndrome], n.is_isolation_attribute_of, n[node]))

    # remember, there can be unknown serovars.
    # Htype is just ID for now
    # TODO: standardize H/O types? Make sure they ALL keep using id #'s only
    add_Htype(name, Htype)

    add_Otype(name, Otype)

    if User is not None:
        g.add( (n[name], n.is_owned_by, n[User]) )
        g.add( (n[User], n.owns, n[name]) )


def add_Htype(name, Htype=None):
    if Htype is None:
        g.add((n[name], n.has_Htype, n.HUnknown))
        g.add((n.HUnknown, n.is_Htype_of, n[name]))
    else:
        g.add((n[name], n.has_Htype, n["H" + str(Htype)]))
        g.add((n["H" + str(Htype)], n.is_Htype_of, n[name]))


def add_Otype(name, Otype=None):
    if Otype is None:
        g.add((n[name], n.has_Otype, n.OUnknown))
        g.add((n.OUnknown, n.is_Otype_of, n[name]))
    else:
        g.add((n[name], n.has_Otype, n["O" + str(Otype)]))
        g.add((n["O" + str(Otype)], n.is_Otype_of, n[name]))


# TODO: definitely need to simplify and refactor genome uploader!!!
def create_pending_genome(name, date = None, location = None, accession = None, bioproject = None, biosample = None, strain = None,
                          organism = None, from_host = None, from_source = None, syndrome = None, Htype = None, Otype = None, User = None ):
    create_genome(name, date, location, accession, bioproject, biosample, strain, organism, from_host, from_source, syndrome, Htype, Otype, User)
    g.add( (n[name], rdf.type, n.pending_genome) )


def create_completed_genome(name, date = None, location = None, accession = None, bioproject = None, biosample = None, strain = None,
                            organism = None, from_host = None, from_source = None, syndrome = None, Htype = None, Otype = None, User = None):
    create_genome(name, date, location, accession, bioproject, biosample, strain, organism, from_host, from_source, syndrome, Htype, Otype, User)
    g.add( (n[name], rdf.type, n.completed_genome) )


def create_isolation_syndrome(name, label, host_category):
    create_class(name)
    g.add( (n[name], rdf.type, n.isolation_syndrome) )
    g.add( (n[name], rdfs.label, Literal(str(label), datatype = XSD.string)) )
    g.add( (n[name], n.has_host_category, n[host_category]) )
    g.add( (n[host_category], n.is_host_category_of, n[name]) )



def create_from_isolation_source(name, label, host_category):
    create_class(name)
    g.add( (n[name], rdf.type, n.isolation_from_source))
    g.add( (n[name], rdfs.label, Literal(str(label), datatype = XSD.string)) )
    g.add( (n[name], n.has_host_category, n[host_category]) )
    g.add( (n[host_category], n.is_host_category_of, n[name]) )


def create_host_category(name, label):
    create_class(name)
    g.add( (n[name], rdf.type, n.host_category) )
    g.add( (n[name], rdfs.label, Literal(str(label), datatype = XSD.string)) )


# can create any number of Otypes as required (1-187 and unknown)
def create_Otype(id):
    create_class("O" + str(id))
    g.add( (n["O" + str(id)], rdf.type, n["Otype"]) )

# can create any number of Htypes as required (1-56, -, and unknown)
def create_Htype(id):
    create_class("H" + str(id))
    g.add( (n["H" + str(id)], rdf.type, n["Htype"]) )


def generate_output(destination):
    g.serialize(destination=destination,format="turtle")

