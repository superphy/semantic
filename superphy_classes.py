__author__ = 'ubiquitin'

from rdflib import Graph, Namespace, Literal, XSD
import superphySPARQL

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


class NamedIndividual(object):
    def __init__(self, name):
        self.name = name

    def rdf(self):
        g.add( (n[self.name], rdf.type, owl.NamedIndividual) )


class Organism(NamedIndividual):
    def __init__(self, name, label, scientific_name, common_name, taxonomy_id=None):
        super(Organism, self).__init__(name)
        self.label = label
        self.scientific_name = scientific_name
        self.common_name = common_name
        self.taxonomy_id = taxonomy_id


    def rdf(self):
        super(Organism, self).rdf()
        g.add( (n[self.name], rdf.type, n.Organism) )
        g.add( (n[self.name], rdfs.label, Literal(str(self.label), datatype = XSD.string)) )
        g.add( (n[self.name], n.scientific_name, Literal(str(self.scientific_name), datatype = XSD.string)) )
        g.add( (n[self.name], n.common_name, Literal(str(self.common_name), datatype = XSD.string)) )
        g.add( (n[self.name], n.has_taxonomy_id, Literal(str(self.taxonomy_id), datatype = XSD.string)) )

class Host(Organism):
    def __init__(self, name, label, scientific_name, common_name, host_category, taxonomy_id=None):
        super(Host, self).__init__(name, label, scientific_name, common_name, taxonomy_id)
        self.host_category = host_category


    def rdf(self):
        super(Host, self).rdf()
        g.add( (n[self.name], rdf.type, n.Host) )
        g.add( (n[self.name], n.has_host_category, n[self.host_category]) )
        g.add( (n[self.host_category], n.is_host_category_of, n[self.name]) )
        FromHost(self.name, self.host_category).rdf()

class Microbe(Organism):
    def rdf(self):
        super(Microbe, self).rdf()
        g.add( (n[self.name], rdf.type, n.Microbe) )

class Attribute(NamedIndividual):
    def rdf(self):
        super(Attribute, self).rdf()
        g.add( (n[self.name], rdf.type, n.Attribute) )

class HostCategory(Attribute):
    def __init__(self, name, label):
        super(HostCategory, self).__init__(name)
        self.label = label

    def rdf(self):
        super(Attribute, self).rdf()
        g.add( (n[self.name], rdf.type, n.host_category) )
        g.add( (n[self.name], rdfs.label, Literal(str(self.label), datatype = XSD.string)) )

class IsolationAttribute(Attribute):
    def rdf(self):
        super(IsolationAttribute, self).rdf()
        g.add( (n[self.name], rdf.type, n.isolation_attribute) )

class FromHost(IsolationAttribute):
    def __init__(self, host, host_category):
        self.name = "from_" + host
        super(FromHost, self).__init__(self.name)
        self.host = host
        self.host_category = host_category

    def rdf(self):
        super(FromHost, self).rdf()
        g.add( (n[self.name], rdf.type, n.isolation_from_host) )
        g.add( (n[self.name], n.has_object, n[self.host]) )
        g.add( (n[self.host], n.is_object_of, n[self.name]))
        g.add( (n[self.name], n.has_host_category, n[self.host_category]) )
        g.add( (n[self.host_category], n.is_host_category_of, n[self.name]) )

class FromSource(IsolationAttribute):
    def __init__(self, name, label, host_category):
        super(FromSource, self).__init__(name)
        self.label = label
        self.host_category = host_category

    def rdf(self):
        g.add( (n[self.name], rdf.type, n.isolation_from_source) )
        g.add( (n[self.name], rdfs.label, Literal(str(self.label), datatype = XSD.string)) )
        g.add( (n[self.name], n.has_host_category, n[self.host_category]) )
        g.add( (n[self.host_category], n.is_host_category_of, n[self.name]) )

class IsolationSyndrome(IsolationAttribute):
    def __init__(self, name, label, host_category):
        super(FromSource, self).__init__(name)
        self.label = label
        self.host_category = host_category

    def rdf(self):
        g.add( (n[self.name], rdf.type, n.isolation_syndrome) )
        g.add( (n[self.name], rdfs.label, Literal(str(self.label), datatype = XSD.string)) )
        g.add( (n[self.name], n.has_host_category, n[self.host_category]) )
        g.add( (n[self.host_category], n.is_host_category_of, n[self.name]) )

class Serotype(Attribute):
    def rdf(self):
        super(Serotype, self).rdf()
        g.add( (n[self.name], rdf.type, n.serotype) )

class Otype(Serotype):
    def __init__(self, id):
        self.name = "O" + str(id)
        super(Otype, self).__init__((self.name))

    def rdf(self):
        super(Otype, self).rdf()
        g.add( (n[self.name], rdf.type, n.Otype) )

class Htype(Serotype):
    def __init__(self, id):
        self.name = "H" + str(id)
        super(Htype, self).__init__((self.name))

    def rdf(self):
        super(Htype, self).rdf()
        g.add( (n[self.name], rdf.type, n.Htype) )

class Genome(NamedIndividual):

    def __init__(self, name, **kwargs):
        SearchParam = ["date", "location", "accession", "bioproject", "biosample", "strain", "organism", "from_host",
                       "from_source", "syndrome", "Htype", "Otype", "User"]

        super(Genome, self).__init__(name)
        self.kwargs = {key:value for key, value in kwargs.items() if key in SearchParam}

    def rdf(self):
        super(Genome, self).rdf()

        if "Htype" not in self.kwargs:
            self.Htype()
        if "Otype" not in self.kwargs:
            self.Otype()

        for key, value in self.kwargs.iteritems():
            getattr(self, key)(value)

    def date(self, date):
        [g.add( (n[self.name], n.has_isolation_date, Literal(item, datatype = XSD.date)) ) for item in date]

    def location(self, location):
        [g.add( (n[self.name], n.has_geographic_location, Literal(item, datatype = XSD.string)) ) for item in location]

    def accession(self, accession):
        [g.add( (n[self.name], n.has_accession, Literal(item, datatype = XSD.string)) ) for item in accession]

    def bioproject(self, bioproject):
        [g.add( (n[self.name], n.has_bioproject, Literal(item, datatype = XSD.string)) ) for item in bioproject]

    def biosample(self, biosample):
        [g.add( (n[self.name], n.has_biosample, Literal(item, datatype = XSD.string)) ) for item in biosample]

    def strain(self, strain):
        [g.add( (n[self.name], n.has_strain, Literal(item, datatype = XSD.string)) ) for item in strain]

    def organism(self, organism):
        g.add( (n[self.name], n.is_genome_of, n[organism]))
        g.add( (n[organism], n.has_genome, n[self.name]))

    def from_host(self, from_host):
        for item in from_host:
            node = superphySPARQL.find_from_host(item).split("#", 1)[1]
            g.add( (n[self.name], n.has_isolation_attribute, n[node]) )
            g.add( (n[node], n.is_isolation_attribute_of, n[self.name]) )

    def from_source(self, from_source):
        for item in from_source:
            node = superphySPARQL.find_source(item).split("#", 1)[1]
            g.add( (n[self.name], n.has_isolation_attribute, n[node]) )
            g.add( (n[node], n.is_isolation_attribute_of, n[self.name]) )

    def syndrome(self, syndrome):
        for item in syndrome:
            node = superphySPARQL.find_syndrome(item).split("#", 1)[1]
            g.add( (n[self.name], n.has_isolation_attribute, n[node]) )
            g.add( (n[node], n.is_isolation_attribute_of, n[self.name]) )

    def Htype(self, Htype=None):
        if Htype is None:
            g.add((n[self.name], n.has_Htype, n.HUnknown))
            g.add((n.HUnknown, n.is_Htype_of, n[self.name]))
        else:
            g.add((n[self.name], n.has_Htype, n["H" + str(Htype)]))
            g.add((n["H" + str(Htype)], n.is_Htype_of, n[self.name]))

    def Otype(self, Otype=None):
        if Otype is None:
            g.add((n[self.name], n.has_Otype, n.OUnknown))
            g.add((n.OUnknown, n.is_Otype_of, n[self.name]))
        else:
            g.add((n[self.name], n.has_Otype, n["O" + str(Otype)]))
            g.add((n["O" + str(Otype)], n.is_Otype_of, n[self.name]))

    def User(self, User):
        g.add( (n[self.name], n.is_owned_by, n[User]) )
        g.add( (n[User], n.owns, n[self.name]) )

class PendingGenome(Genome):
    def rdf(self):
        super(PendingGenome, self).rdf()
        g.add( (n[self.name], rdf.type, n.pending_genome) )

class CompletedGenome(Genome):
    def rdf(self):
        super(CompletedGenome, self).rdf()
        g.add( (n[self.name], rdf.type, n.completed_genome) )

def generate_output(destination):
    g.serialize(destination=destination,format="turtle")



""" ================================================================================================================ """

item = Host("hsapiens", "Homo sapiens (human)", "Homo sapiens", "human", "human")
item.rdf()
item = Microbe("ecoli", "Escherichia coli (E. coli)", "Escherichia coli", "E. coli")
item.rdf()

"""
list(Htype(num).rdf()  for num in xrange(1,56))
list(Otype(num).rdf() for num in xrange(1,187))
"""

genome = Genome("NC_455763", **{"date":{"2010-10-09"}, "location":{"Canada"}, "Htype":"7", "from_host":{"Homo sapiens (human)"}, "organism": "ecoli", "from_source":{"Blood"}})
print genome.kwargs
genome.rdf()



g.serialize(destination="test.txt", format="turtle")

