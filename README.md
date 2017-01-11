# semantic   [![Build Status](https://travis-ci.org/superphy/semantic.svg?branch=master)](https://travis-ci.org/superphy/semantic)
SuperPhy for the semantic web

# Starting from a fresh install of linuxmint-17.3-cinnamon-64bit

    * sudo apt-get update && sudo apt-get upgrade -y
    * sudo apt-get install apache2 curl git libapache2-mod-wsgi libyajl2 mummer muscle python-dev python-virtualenv wget xvfb -y
    * git clone https://github.com/superphy/semantic.git ~/superphy
    * cd ~/superphy
    * python setup.py install --sys
    *
    * //to start blazegraph
    * python setup.py run
    *
    * cd app/
    * //to start superphy webapp
    * python run.py
    * //Navigate to localhost:5000, and verify the page loads, go back to the terminal and hit ctrl-c
    * //Navigate to localhost, and verify the apache-mod_wsgi server is running.

/*If any thing doesn't load, try refreshing the page again, as the first load is broken at the time of writing this.

Congratulations. You have a working version of SuperPhy. Next thing to do is add in genome data to the Blazegraph Triplestore.
    

## [Coffeescript](http://coffeescript.org/)
* Setup coffeelint to ensure we are all following the same coding style
* coffeelint.json in the root directory ensures the rules apply to all subdirectories


## [Mithril 0.2.0](http://lhorie.github.io/mithril/)
The MVC framework is called Mithril, in Javascript.
It is minimalistic and can be used more like a library than a huge framework.
It is small, fast and allows easy integration with other libraries.
Like react, it provides a virtual DOM, and only re-draws what has changed.

## [LESS](http://lesscss.org/features/)
- CSS pre-processor

## [Blazegraph TripleStore](https://www.blazegraph.com/)
A standalone Blazegraph `.jar` file can be found at [Blazegraph Download](www.blazegraph.com/download/)

## Adding real genomic data to the Blazegraph triplestore
The Java Heap Size will need to be changed (probably) for datasets of our size -- even loading the basic GO owl file required this. 
This will give blazegraph a 4GB heap space, which should be enough.

The OWL file for the Genomic Feature and Variation Ontology (GFVO) An ontology for describing genomic features and variants; in particular the contents of GFF3, GTF, GVF and VCF files, should be obatined from: (https://raw.github.com/BioInterchange/Ontologies/master/gfvo.xml) and installed


    curl -X POST -H 'Content-Type:application/rdf+xml' --data-binary '@gfvo.xml' http://localhost:8080/blazegraph/namespace/superphy/sparql

GO ontology for genes in OWL (http://geneontology.org/page/download-ontology). 10 tips for using GO (www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1003343#s2). The latest version includes the ~900 PAMGO (http://pamgo.vbi.vt.edu/) terms for virulence


    curl -X POST -H 'Content-Type:application/rdf+xml' --data-binary '@go.owl' http://localhost:8080/blazegraph/namespace/superphy/sparql

This takes a few minutes (Modified: 1349345 Milliseconds: 926484)

The genome should undergo the quality checks before being loaded into the triplestore. Following addition, each new analysis should query the triplestore for the sequence, and then update the triplestore with the result. This allows things to be done in parallel, rather than one monolithic analysis that can fail. Also, using GET uses caching, so their shouldn't be a huge amount of data transfer (eg. querying multiple times for the same genome sequence)

## Ontology
An ontology describes hierarchical relationships, and the properties.
To create an actual RDF instance, you can use the rdf:type  -- this is so common that the agreed upon shorthand is the single letter 'a' -- it means exactly the same thing.
See (http://www.linkeddatatools.com/introducing-rdfs-owl) for a well-explained introduction using OWL and RDF.

# Directory structure for semantic superphy (OUTDATED)

# Top level

* - Directory is expanded below top level

Directory.........................Description

app/*.............................Flask appication files; files for creating the website
    auth/*........................views, and forms having to do with the authentication content
    main/*........................views, and forms having to do with the main content
    static/*......................Static assets: pictures, css, js, cscript.
    templates/*...................Jinja2 templates. .HTML files.
    api_1_0/*.....................Relic from Flasky skeleton. RESTful API Blueprint
bash/*............................bash scripts for matenance and setup
migrations/.......................Relic from Flasky skeleton. Iterative deployment of the table SQL database is stored here.
ontology/.........................Linked data files (TTL,RDF,etc). These were made with 'Protege'. They provide the framework for adding linked data to the triplestore.
python_package/...................Individual packages 
    scrypt/.......................?Script?
    sparql/.......................Sparql queries & updates. sparql endpoint interface.
    uploader/.....................Genome import module
tests/............................Testing code for the entire project
tmp/..............................Temporary storage. This folder mostly ignored by git
    findhome/.....................Stranded code. This folder is not ignored by git.
uml/..............................Specification, Visualization, Construction and Documentation artifacts
venv/.............................Virtual Environment to keep the dependencies required by the project separate from the dev system.
    requirments.txt...............Requirements install file for the virtual environment

# app/

app/api_1_0/
    authentication.py.............flask.ext.httpauth authentication
    comments.py...................Relic from Flasky skeleton.
    decorators.py.................Function Decorators
    errors.py.....................Error handling
    __init__.py...................API blueprint constructor
    posts.py......................Relic from Flasky skeleton.
    users.py......................Relic from Flasky skeleton. (Could keep track of some things here)

app/auth/
    forms.py......................Authentication Form Classes
    __init__.py...................blueprint constructor
    views.py......................Authentication View Functions

app/main/
    errors.py.....................Main Error View Functions (Keep this)
    forms.py......................Relic from Flasky skeleton.
    __init__.py...................blueprint constructor
    views.py......................Relic from Flasky skeleton.

app/static/
    coffee/.......................all coffeescript files for creating the website
    css/..........................all CSS files (generated by LESS, or 3rd party)
    images/.......................all images for the site
    js/...........................generated javascript (from Coffee files) and 3rd party
    less/.........................style coding, generated CSS is in separate folder
    coffeelint.json...............coffeescript style conventions
    gulpfile.js...................config file for Gulp
    index.html....................Mithril home (includes HTML sekelton and includes minified JS and CSS)

app/templates/
    auth/.........................Authentication Webpages
        change_email.html.........^
        change_password.html......^
        email/....................Email related webpages
            change_email.html.....^
            change_email.txt......^
            confirm.html..........^
            confirm.txt...........^
            reset_password.html...^
            reset_password.txt....^
        login.html................Login Webpage
        register.html.............Register Webpage
        reset_password.html.......Password Reset Webpage
        unconfirmed.html..........Unconfirmed account Webpage
    mail/.........................Relic from Flasky skeleton.
        new_user.html.............Relic from Flasky skeleton.
        new_user.txt..............Relic from Flasky skeleton.
    403.html......................403 Error Webpage
    404.html......................404 Error Webpage
    500.html......................500 Error Webpage
    base.html.....................Base webpage for many other pages to inherit.
    _comments.html................Relic from Flasky skeleton.
    edit_post.html................Relic from Flasky skeleton.
    edit_profile.html.............Relic from Flasky skeleton?
    error_page.html...............Base webpage for many other pages to inherit.
    followers.html................Relic from Flasky skeleton.
    index.html....................Index
    _macros.html..................Relic from Flasky skeleton?
    moderate.html.................Relic from Flasky skeleton?
    post.html.....................Relic from Flasky skeleton.
    _posts.html...................Relic from Flasky skeleton.
    user.html.....................Relic from Flasky skeleton.

app/decorators.py.................Function Decorators
app/email.py......................Sending email (broken)
app/exceptions.py.................Relic from Flasky skeleton. Dummy?
app/__init.py.....................Main app init function

#
