#!/bin/bash
#Run chmod a+x init.sh to make this exucutable
#Need to have virtualenv installed
#Sets up a virtual environment, this allow you to install dependancies without breaking your build.
mkdir superphy &> /dev/null
if find superphy -maxdepth 0 -empty | read v; then 
	virtualenv superphy;
	superphy/binpip install -r /path/to/requirements.txt &> /dev/null
fi
#Downloads the blazegraph client from sourceforge. This will have to be updated when we migrate to a new version of BG.
mkdir db &> /dev/null
if ! find db/bigdata-bundled.jar | read v; then
	cd db;
	wget "http://iweb.dl.sourceforge.net/project/bigdata/bigdata/1.5.3/bigdata-bundled.jar"; 
	cd ..;
fi