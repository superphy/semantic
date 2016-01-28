#!/bin/bash

#Config file
if ! find ../config/superphy.cfg | read v; then
	echo "Creating ../config ..."
	mkdir ../config
	IP="$(ifconfig | grep -A 1 'eth0' | tail -1 | cut -d ':' -f 2 | cut -d ' ' -f 1)"
	#IP="localhost"
	PORT="9000"
	echo "[rdf]" >> ../config/superphy.cfg
	echo "url = http://"$IP":"$PORT"/blazegraph/namespace/superphy/sparql" >> ../config/superphy.cfg
fi

#sudo installs
read -r -p "Install sudo packages to your dev system? [Y/N] " response
case $response in
    [yY][eE][sS]|[yY])
		#Put sudo packages that install to the development machine here. There should be no npm or pip here, because that is installed to the virtual environment below.
		sudo apt-get install python-virtualenv
		sudo apt-get install xvfb
		sudo apt-get install libyajl1
		sudo apt-get install libyajl2
        ;;
    *)
        echo "Not installing sudo packages."
        ;;
esac

#deploy apache
bash bash/deploy_apache.sh superphy

#Virtualenvironment install
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

#Getting the BLAST+ gzip file from remote server.
mkdir venv/lib/blast &> /dev/null || true
if ! find venv/lib/blast/ncbi*/ | read v; then

	if find ../blast/ncbi*/ | read v;
	then
		echo Copying NCBI BLAST+ jar file...
		cp -r ../blast/ncbi*/ venv/lib/blast/
	else
		echo Downloading NCBI BLAST+ jar file...
		cd venv/lib/blast;
	    wget -rc -nd -A "*x64-linux.tar.gz" "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
		tar zxvpf *x64-linux.tar.gz

		cd ../../..;
		pwd
		read -r -p "Do you want to backup the NCBI BLAST+ client in ../blast ? [Y/N] " response
		case $response in
		    [yY][eE][sS]|[yY])
			mkdir ../blast || true #Redundant?
			cp -r venv/lib/blast/ncbi*/ ../blast/
			;;*)
		esac
	fi
fi

#Downloading Graph database
cd database
	bash scripts/start.sh
cd ..

echo Finished
echo """$ bash scripts/run.sh""" to run the server
exit 0