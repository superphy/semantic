#!/bin/bash

chmod a+x *.py

if ! find venv/bin/pip | read v; then
	echo Setting up Virtual Environment
	virtualenv --no-site-packages venv

    # Add environment setup scripts to activate
    BASHDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    echo Adding postactivate script to venv/bin/activate
    echo "# Superphy Environment setup" >> venv/bin/activate
    echo "source $BASHDIR/postactivate.sh" >> venv/bin/activate
fi

source venv/bin/activate
pip install --upgrade pip
pip install -r venv/requirements.txt

if ! find venv/bin/node | read v; then
    nodeenv -p --prebuilt --requirements=venv/npm-requirements.txt
fi

deactivate

: <<'END'
#Setting up sqlite server for user auth
if ! find data-dev.sqlite | read v; then
	echo Setting up SQL server
	./manage.py db upgrade &> /dev/null
fi
echo sqlite server setup complete!

#Getting the graph db jar file from remote server.
mkdir db
if ! [ -s db/bigdata-bundled.jar ]; then
    rm db/bigdata-bundled.jar
	echo Downloading Blazegraph jar file...
	cd db;
	wget "http://iweb.dl.sourceforge.net/project/bigdata/bigdata/1.5.3/bigdata-bundled.jar";
	cd ..;
fi
echo Blazegraph setup complete!

#Getting the BLAST+ gzip file from remote server.
mkdir blast
if ! find blast/ncbi*/ | read v; then
	echo Downloading NCBI BLAST+ jar file...
	cd blast;
	wget -rc -nd -A "*x64-linux.tar.gz" "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
	tar zxvpf *x64-linux.tar.gz
	cd ..;
fi
echo BLAST+ setup complete!

#Run related scripts
bash bash/start_blazegraph
END
echo Finished
echo """$ bash bash/run""" to run the server
exit 0