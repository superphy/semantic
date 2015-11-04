# semantic   [![Build Status](https://travis-ci.org/superphy/semantic.svg?branch=master)](https://travis-ci.org/superphy/semantic)
SuperPhy for the semantic web     

# Quick Start
- Install Flask
    * bash bash/install
    * ./manage.py runserver
- Start Blazegraph
    * `bash bash/start_blazegraph`
    * to stop it use `bash bash/kill_port_9999`
- Running Gulp will compile all of the LESS to CSS, convert Coffeescript to JS, combine and minify all files
- The index.html file therefore only needs to load `mithril.min.js` and `all.min.js`, where all.min.js is the combined and minified file of all the Coffeescript in the project.
- Likewise, the css files included are `bootstrap.min.css` and `all.min.css`, where all.min.cs is the LESS --> CSS, combined and minified.
- All files created in the `coffeee` and `LESS` folders get automatically included
- gulpfile.js specifies what is done when gulp is run. It can also run tests etc.
- Your IDE / editor and be setup to run Gulp whenever `build` is run, which is very convenient
- point the root for the website to the root of the repository
- the `scrypt` python module is in the `scripts` directory


# Semantic SuperPhy Overview
- Project overview is in `uml/project_overview.puml`
    * note that the GraphViz executable is called `dot`

## [Coffeescript](http://coffeescript.org/)
* Setup coffeelint to ensure we are all following the same coding style
* coffeelint.json in the root directory ensures the rules apply to all subdirectories


## [Mithril 0.2.0](http://lhorie.github.io/mithril/)
The MVC framework is called Mithril, in Javascript.
It is minimalistic and can be used more like a library than a huge framework.
It is small, fast and allows easy integration with other libraries.
Like react, it provides a virtual DOM, and only re-draws what has changed.


## [Gulp.js](http://gulpjs.com/)
Used to compile and minify everything for the site
- Install nodejs and npm
- Install all of the dependencies
    * npm install gulp
    * npm install gulp-coffee
    * npm install gutil
    * npm install gulp-uglify
    * npm install gulp-rename
    * npm install gulp-less
    * npm install gulp-minify-css

## [LESS](http://lesscss.org/features/)
- CSS pre-processor


### For this version:
Use Mithril

Use Gulp

Use LESS

Use Bootstrap.css

Use IntroJS for page introduction


Note that valid JSON requires the property names to be quoted.
[Link](http://stackoverflow.com/questions/15699196/syntaxerror-json-parse-expected-property-name-or-while-using-highcharts)
eg.
```json
    good: [{"name": "John"}, {"name": "Mary"}]
    bad: [{name: "John"}, {name: "Mary"}]
```

## [Blazegraph TripleStore](https://www.blazegraph.com/)
A standalone Blazegraph `.jar` file can be found at [Blazegraph Download](www.blazegraph.com/download/)

### `.war` file setup
Download the latest Blazegraph .war [file:](http://sourceforge.net/projects/bigdata/)

### Configure Tomcat
A persistent Blazegraph requires Apache Tomcat. On Debian, Tomcat installs to `/var/lib/tomcat8`; the webapps directory is under this root
If you have not changed any configuration files, please examine the file conf/tomcat-users.xml in your installation. That file must contain the credentials to let you use this webapp.
For example, to add the manager-gui role to a user named tomcat with a password of s3cret, add the following to the config file listed above.
```
    <role rolename="manager-gui"/>
    <user username="tomcat" password="s3cret" roles="manager-gui"/>
```
Make sure to `service tomcat8 restart` after making changes to any configuration file.
After configuring a manager in Tomcat, use the Tomcat manager (`localhost:8080/manager`) to deploy and start the .war file.
This will not work initially, the following two changes need to be made from the webapps/bigdata directory:
 1. Modify web.xml to use an absolute path for the .properties file, eg.:
```
    <param-name>propertyFile</param-name>
    <param-value>/var/lib/tomcat8/webapps/bigdata/WEB-INF/RWStore.properties</param-value>
```
2. Modify the RWStore.properties file to use an absolute path for the triplestore database itself, eg.:
    `com.bigdata.journal.AbstractJournal.file=/var/lib/tomcat8/webapps/bigdata/bigdata.jnl`

After modifying these two files, start the deployed app via the Tomcat manager page


## Using Blazegraph
Navigate to `http://localhost:8080/bigdata/#namespaces` and create a new namespace (eg. superphy).
Blazegraph supports SPARQL 1.1, so all CRUD activities can be performed via a REST API.
It is important to specify the correct format of the query in the header, otherwise the action will fail. For example, to insert some nonsense data into the triplestore in turtle format (taken from `http://www.w3.org/TeamSubmission/turtle/`), create a `test.ttl` file with the following:
```turtle
    # this is a complete turtle document
    @prefix foo: <http://example.org/ns#> .
    @prefix : <http://other.example.org/ns#> .
    foo:bar foo: : .
    :bar : foo:bar .
```
And use curl to add this via the REST API:

    curl -X POST -H 'Content-Type:application/x-turtle' --data-binary '@test.ttl' http://localhost:8080/bigdata/sparql

Note that the default namespace is 'kb' and can be accessed at

    http://localhost:8080/bigdata/sparql

To use a specific namespace, it is specified as follows:

    http://localhost:8080/bigdata/namespace/kb/sparql

For superphy, it would be:

    http://localhost:8080/bigdata/namespace/superphy/sparql

If the file was not in turtle format, or the Content-Type was wrong, the action fails. See `https://wiki.blazegraph.com/wiki/index.php/NanoSparqlServer#INSERT` for MIME types, and the REST API for Blazegraph


The following simple query will return one statement from the default KB instance:

    curl -X GET http://localhost:8080/bigdata/namespace/superphy/sparql --data-urlencode 'query=SELECT * { ?s ?p ?o } LIMIT 1' -H 'Accept:application/rdf+xml'


If you want the result set in JSON, use:

    curl -X GET http://localhost:8080/bigdata/namespace/superphy/sparql --data-urlencode 'query=SELECT * { ?s ?p ?o } LIMIT 1' -H 'Accept:application/sparql-results+json'

In order to use SPARQL syntax to perform insert or update, one must use the 'update', rather than 'query' syntax. For example, to insert the example given at `http://www.w3.org/TR/2013/REC-sparql11-update-20130321/#updateLanguage`, which is below from their 'Adding some triples to a graph' section:
```
    PREFIX dc: <http://purl.org/dc/elements/1.1/>
    INSERT DATA
    {
      <http://example/book1> dc:title "A new book" ;
                             dc:creator "A.N.Other" .
    }
```
Use the following command:

    curl -X POST http://localhost:8080/bigdata/namespace/superphy/sparql --data-urlencode 'update=PREFIX dc: <http://purl.org/dc/elements/1.1/> INSERT DATA { <http://example/book1> dc:title "A new book" ; dc:creator "A.N.Other" . }' -H 'Accept:application/sparql-results+json'


## Adding real genomic data to the Blazegraph triplestore
The Java Heap Size will need to be changed (probably) for datasets of our size -- even loading the basic GO owl file required this. To do so, locate the base tomcat directory. On Debian it is /usr/share/tomcat8/bin/ (there should be a startup.sh file in this directory). The catalina.sh file is the one invoked for Tomcat startup -- in it, it describes the following: create a new file in the same directory as catalina.sh called setenv.sh and specify the environment parameters in it. eg:
    export JAVA_OPTS="-Xmx4g"
This will give bigdata a 4GB heap space, which should be enough.

The OWL file for the Genomic Feature and Variation Ontology (GFVO) An ontology for describing genomic features and variants; in particular the contents of GFF3, GTF, GVF and VCF files, should be obatined from: (https://raw.github.com/BioInterchange/Ontologies/master/gfvo.xml) and installed


    curl -X POST -H 'Content-Type:application/rdf+xml' --data-binary '@gfvo.xml' http://localhost:8080/bigdata/namespace/superphy/sparql

GO ontology for genes in OWL (http://geneontology.org/page/download-ontology). 10 tips for using GO (www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1003343#s2). The latest version includes the ~900 PAMGO (http://pamgo.vbi.vt.edu/) terms for virulence


    curl -X POST -H 'Content-Type:application/rdf+xml' --data-binary '@go.owl' http://localhost:8080/bigdata/namespace/superphy/sparql

This takes a few minutes (Modified: 1349345 Milliseconds: 926484)

The genome should undergo the quality checks before being loaded into the triplestore. Following addition, each new analysis should query the triplestore for the sequence, and then update the triplestore with the result. This allows things to be done in parallel, rather than one monolithic analysis that can fail. Also, using GET uses caching, so their shouldn't be a huge amount of data transfer (eg. querying multiple times for the same genome sequence)

## Ontology
An ontology describes hierarchical relationships, and the properties.
To create an actual RDF instance, you can use the rdf:type  -- this is so common that the agreed upon shorthand is the single letter 'a' -- it means exactly the same thing.
See (http://www.linkeddatatools.com/introducing-rdfs-owl) for a well-explained introduction using OWL and RDF.




