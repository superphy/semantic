__author__ = 'ubiquitin'

from rdflib import Graph, Namespace, Literal, XSD

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
        fh = FromHost(self.name)

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









""" ================================================================================================================ """

item = Host("hsapiens", "Homo sapiens (human)", "Homo sapiens", "human")
item.rdf()
item = Microbe("ecoli", "Escherichia coli (E. coli)", "Escherichia coli", "E. coli")
item.rdf()

g.serialize(destination="test.txt", format="turtle")
