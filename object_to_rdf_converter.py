__author__ = 'ubiquitin'

from rdflib import Graph, Namespace, Literal, XSD

# initialize a RDF Graph
g = Graph()

# setting up namespaces for use
n = Namespace("https://github.com/superphy/")
owl = Namespace("http://www.w3.org/2002/07/owl#")
rdf = Namespace("http://www.w3.org/1999/02/22-rdf-syntax-ns#")
xml = Namespace("http://www.w3.org/XML/1998/namespace")
xsd = Namespace("http://www.w3.org/2001/XMLSchema#")
rdfs = Namespace("http://www.w3.org/2000/01/rdf-schema#")
gfvo = Namespace("http://www.biointerchange.org/gfvo#")

g.bind("","https://github.com/superphy/")
g.bind("owl", "http://www.w3.org/2002/07/owl#")

def create_class(name):
    g.add( (n[name], rdf.type, owl.NamedIndividual) )


def create_organism(name, taxonomy_id=None):
    create_class(name)
    g.add( (n[name], n.has_taxonomy_id, Literal(str(taxonomy_id), datatype = XSD.string)) )


# assume: host_category is created
# TODO: add in a check for that
def create_host(name, host_category, taxonomy_id=None):
    create_organism(name, taxonomy_id)
    g.add( (n[name], rdf.type, n.Host) )
    g.add( (n[name], n.has_host_category, n[host_category]) )
    g.add( (n[host_category], n.is_host_category_of, n[name]) )
    from_host = create_from_host(name)
    g.add( (n[name], n.is_object_of, n[from_host]) )
    g.add( (n[from_host], n.has_object, n[name]) )


def create_from_host(host):
    name = "from_" + host
    create_class(name)
    g.add( (n[name], rdf.type, n.isolation_from_host) )
    return name


def create_microbe(name, taxonomy_id=None):
    create_organism(name, taxonomy_id)
    g.add( (n[name], rdf.type, n.Microbe) )


# assumption: all other objects must have already been created; the genome metadata is always added last
# remember, date must be in dateTime format (YYYY-MM-DDTHH:MM:SS)
# TODO: refactor code to be less clunky (so much duplication!) and have more checks for created items
# TODO: include environment as a genome source (under organism or under source???)
def create_genome(name, date = None, location = None, accession = None, bioproject = None, biosample = None, strain = None,
                  organism = None, from_host = None, from_source = None, syndrome = None, Htype = None, Otype =None ):
    create_class(name)
    g.add( (n[name], rdf.type, gfvo.Genome) )

    if date is not None:
        g.add( (n[name], n.has_isolation_date, Literal(date, datatype = XSD.dateTime)) )

    if location is not None:
        g.add ( (n[name], n.has_geographic_location, Literal(location, datatype = XSD.string)) )

    if accession is not None:
        g.add( (n[name], n.has_accession, Literal(accession, datatype = XSD.string)) )

    if bioproject is not None:
        g.add( (n[name], n.has_bioproject, Literal(bioproject, datatype = XSD.string)) )

    if biosample is not None:
        g.add( (n[name], n.has_biosample, Literal(biosample, datatype = XSD.string)) )

    if strain is not None:
        g.add( (n[name], n.has_strain, Literal(strain, datatype = XSD.string)) )

    # TODO: need to implement a way to check if organism is created, else return a failure
    if organism is not None:
        g.add( (n[name], n.is_genome_of, n[organism]) )
        g.add( (n[organism], n.has_genome, n[name]))

    # TODO: need to implement a way to check if the host and from_host entity is created, else return a failure
    # TODO: decide if this takes in from_host or host (current assumption: from_host)
    if from_host is not None:
        g.add( (n[name], n.has_isolation_attribute, n[from_host]) )
        g.add( (n[from_host], n.is_isolation_attribute_of, n[name]) )

    if from_source is not None:
        g.add( (n[name], n.has_isolation_attribute, n[from_source]) )
        g.add( (n[from_source], n.is_isolation_attribute_of, n[name]) )

    if syndrome is not None:
        g.add( (n[name], n.has_isolation_attribute, n[syndrome]) )
        g.add( (n[syndrome], n.is_isolation_attribute_of, n[name]) )

    # remember, there can be unknown serovars.
    # Htype is just ID for now
    # TODO: standardize H/O types? Make sure they ALL have H/O prefixes or keep using id #'s only
    if Htype is None:
        g.add( (n[name], n.has_Htype, n.HUnknown) )
        g.add( (n.HUnknown, n.is_Htype_of, n[name]) )
    else:
        g.add( (n[name], n.has_Htype, n["H" + str(Htype)]) )
        g.add( (n["H" + str(Htype)], n.is_Htype_of, n[name]) )

    if Otype is None:
        g.add( (n[name], n.has_Otype, n.OUnknown) )
        g.add( (n.OUnknown, n.is_Otype_of, n[name]) )
    else:
        g.add( (n[name], n.has_Otype, n["O" + str(Otype)]) )
        g.add( (n["O" + str(Otype)], n.is_Otype_of, n[name]) )


# TODO: definitely need to simplify and refactor genome uploader!!!
def create_pending_genome(name, date = None, location = None, accession = None, bioproject = None, biosample = None, strain = None,
                  organism = None, from_host = None, from_source = None, syndrome = None, Htype = None, Otype =None ):
    create_genome(name, date, location, accession, bioproject, biosample, strain, organism, from_host, from_source, syndrome, Htype, Otype)
    g.add( (n[name], rdf.type, n.pending_genome) )


def create_completed_genome(name, date = None, location = None, accession = None, bioproject = None, biosample = None, strain = None,
                  organism = None, from_host = None, from_source = None, syndrome = None, Htype = None, Otype =None ):
    create_genome(name, date, location, accession, bioproject, biosample, strain, organism, from_host, from_source, syndrome, Htype, Otype)
    g.add( (n[name], rdf.type, n.completed_genome) )


def create_isolation_syndrome(name, host = None):
    create_class(name)
    g.add( (n[name], rdf.type, n.isolation_syndrome) )

    if host is not None:
        g.add( (n[name], n.has_object, n[host]) )
        g.add( (n[host], n.is_object_of, n[name]) )


def create_from_isolation_source(name):
    create_class(name)
    g.add( (n[name], rdf.type, n.isolation_from_source))


def create_host_category(name):
    create_class(name)
    g.add( (n[name], rdf.type, n.host_category) )


# can create any number of Otypes as required (1-187 and unknown)
def create_Otype(id):
    create_class("O" + str(id))
    g.add( (n["O" + str(id)], rdf.type, n["Otype"]) )

# can create any number of Htypes as required (1-56, -, and unknown)
def create_Htype(id):
    create_class("H" + str(id))
    g.add( (n["H" + str(id)], rdf.type, n["Htype"]) )


def generate_output():
    g.serialize(destination="results.txt",format="turtle")


""" ======== TESTING ======== """

'''
create_Htype(7)
create_Otype(157)
create_Htype("-")
create_Htype("Unknown")
create_Otype("Unknown")
create_microbe("Ecoli", 562)
'''

create_host("Hsapiens", "Mammal")



