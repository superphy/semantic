__author__ = 'ubiquitin'

from superphy_classes import Host, HostCategory, FromSource, IsolationSyndrome, Microbe, Htype, Otype, generate_output

host_categories = [("human", "Human"),
                   ("mammal", "Non-human Mammalia"),
                   ("bird", "Aves"),
                   ("env", "Environmental Sources")]

hosts = [("hsapiens", "Homo sapiens (human)", "Homo sapiens", "human", "human"),
         ("btaurus", "Bos taurus (cow)", "Bos taurus", "cow", "mammal"),
         ("sscrofa", "Sus scrofa (pig)", "Sus scrofa", "pig", "mammal"),
         ("mmusculus", "Mus musculus (mouse)", "Mus musculus", "mouse", "mammal"),
         ("oaries", "Ovis aries (sheep)", "Ovis aries", "sheep", "mammal"),
         ("ggallus", "Gallus gallus (chicken)", "Gallus gallus", "chicken", "bird"),
         ("ocuniculus", "Oryctolagus cuniculus (rabbit)", "Oryctolagus cuniculus", "rabbit", "mammal"),
         ("clupus", "Canis lupus familiaris (dog)", "Canis lupus familiaris", "dog", "mammal"),
         ("fcatus", "Felis catus (cat)", "Felis catus", "cat", "mammal"),
         ("environment", "Environmental source", "Environmental source", "environment", "env"),
         ("eferus", "Equus ferus caballus (horse)", "Equus ferus caballus", "horse", "mammal"),
         ("caegagrus", "Capra aegagrus hircus (goat)", "Capra aegagrus hircus", "goat", "mammal"),
         ("acepa", "Allium cepa (onion)", "Allium cepa", "onion", "env"),
         ("mgallopavo", "Meleagris gallopavo (turkey)", "Meleagris gallopavo", "turkey", "bird")
         ]

sources = [("stool", "Stool", "human"),
           ("urine", "Urine", "human"),
           ("colon", "Colon", "human"),
           ("ileum", "Ileum", "human"),
           ("cecum", "Cecum", "human"),
           ("intestine", "Intestine", "human"),
           ("blood", "Blood", "human"),
           ("liver", "Liver", "human"),
           ("cerebrospinal_fluid", "cerebrospinal_fluid", "human"),
           ("feces", "Feces", "mammal"),
           ("urine", "Urine", "mammal"),
           ("meat", "Meat", "mammal"),
           ("blood", "Blood", "mammal"),
           ("liver", "Liver", "mammal"),
           ("intestine", "Intestine", "mammal"),
           ("feces", "Feces", "bird"),
           ("yolk", "Yolk", "bird"),
           ("meat", "Meat", "bird"),
           ("blood", "Blood", "bird"),
           ("liver", "Liver", "bird"),
           ("veggiefood", "Vegetable-based food", "env"),
           ("meatfood", "Meat-based food", "env"),
           ("water", "Water", "env"),
           ("colon", "Colon", "mammal"),
           ("cecum", "Cecum", "mammal"),
           ("urogenital", "Urogenital system", "human"),
           ("milk", "Milk", "mammal"),
           ("soil", "Soil", "env"),
           ("marine_sediment", "Marine sediment", "env"),
           ("bronchoalveolar_lavage", "Bronchoalveolar lavage", "human"),
           ("gelatinous_edema", "gelatinous edema", "bird"),
           ("perianal", "Perianal", "human"),
           ("enteral_feeding_tube", "Enteral feeding tube", "human"),
           ("wound", "wound", "human"),
           ("wound", "wound", "mammal"),
           ("wound", "wound", "bird"),
           ("abscess", "abscess", "human"),
           ("bile", "bile", "human"),
           ("peritoneal_fluid", "Peritoneal fluid", "human"),
           ("peritoneal_fluid", "Peritoneal fluid", "mammal"),
           ]

syndromes = [("gastroenteritis", "Gastroenteritis", "human"),
             ("bloody_diarrhea", "Bloody diarrhea", "human"),
             ("hus", "Hemolytic-uremic syndrome", "human"),
             ("hc", "Hemorrhagic colitis", "human"),
             ("uti", "Urinary tract infection (cystitis)", "human"),
             ("crohns", "Crohn's Disease", "human"),
             ("uc", "Ulcerateive colitis", "human"),
             ("meningitis", "Meningitis", "human"),
             ("pneumonia", "Pneumonia", "human"),
             ("pyelonephritis", "Pyelonephritis", "human"),
             ("bacteriuria", "Bacteriuria", "human"),
             ("pneumonia", "Pneumonia", "mammal"),
             ("diarrhea", "Diarrhea", "mammal"),
             ("septicaemia", "Septicaemia", "mammal"),
             ("mastitis", "Mastitis", "mammal"),
             ("peritonitis", "Peritonitis", "mammal"),
             ("pneumonia", "Pneumonia", "bird"),
             ("diarrhea", "Diarrhea", "bird"),
             ("septicaemia", "Septicaemia", "bird"),
             ("peritonitis", "Peritonitis", "bird"),
             ("asymptomatic", "Asymptomatic", "human"),
             ("asymptomatic", "Asymptomatic", "mammal"),
             ("asymptomatic", "Asymptomatic", "bird"),
             ("bacteremia", "Bacteremia", "human"),
             ("diarrhea", "Diarrhea", "human"),
             ("septicaemia", "Septicaemia ", "human"),
             ("swollen_head_syndrome", "swollen head syndrome", "bird"),
             ("necrotizing_fasciitis", "Necrotizing fasciitis", "human"),
             ("omphalitis", "Omphalitis", "bird")
             ]

microbes = [("ecoli", "Escherichia coli (E. coli)", "Escherichia coli", "E. coli")]

for host_category in host_categories:
    name, label = host_category
    HostCategory(name, label).rdf()

for host in hosts:
    name, label, sci_name, com_name, host_category = host
    Host(name, host_category, label, sci_name, com_name).rdf()

for source in sources:
    name, label, host_category = source
    FromSource(name, label, host_category).rdf()

for syndrome in syndromes:
    name, label, host_category = syndrome
    IsolationSyndrome(name, label, host_category).rdf()

for microbe in microbes:
    name, label, sci_name, com_name = microbe
    Microbe(name, label, sci_name, com_name).rdf()

list(Htype(num).rdf()  for num in xrange(1,56))
Htype("Unknown").rdf()
Htype("-").rdf()

list(Otype(num).rdf() for num in xrange(1,187))
Otype("Unknown").rdf()


generate_output("setup.ttl")
