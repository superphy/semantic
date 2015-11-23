#!/bin/bash

chmod a+x *.py

if ! find venv/bin/pip | read v; then
	echo Setting up Virtual Environment
	virtualenv --no-site-packages venv
fi

echo "source $TRAVIS_BUILD_DIR/bash/postactivate.sh" >> venv/bin/activate
source venv/bin/activate
pip install --upgrade pip
pip install -r venv/requirements.txt

if ! find venv/bin/node | read v; then
    nodeenv -p --prebuilt --requirements=venv/npm-requirements.txt
fi

deactivate

#Getting the graph db jar file from remote server.
mkdir db &> /dev/null
if ! find db/bigdata-bundled.jar | read v; then

	if find ../db/bigdata-bundled.jar | read v;
	then
		echo Copying Blazegraph jar file...
		cp ../db/bigdata-bundled.jar db/bigdata-bundled.jar
	else
		echo Downloading Blazegraph jar file...
		cd db;
		wget "http://iweb.dl.sourceforge.net/project/bigdata/bigdata/1.5.3/bigdata-bundled.jar";
		cd ..;
		mkdir ../db
		cp db/bigdata-bundled.jar ../db/bigdata-bundled.jar
	fi
fi

#Setting up sqlite server for user auth
if ! find data-dev.sqlite | read v; then
	echo Setting up SQL server
	./manage.py db upgrade &> /dev/null
fi

# Add environment setup scripts to activate
BASHDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo Adding postactivate script to venv/bin/activate
echo "# Superphy Environment setup" >> venv/bin/activate
echo "source $BASHDIR/postactivate.sh" >> venv/bin/activate


#Run related scripts
bash bash/start_blazegraph
bash bash/blastinstall

echo Finished
echo """$ bash bash/run""" to run the server
exit 0