#!/bin/bash
source bash/config
chmod a+x *.py

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
read -r -p "Deploy apache? [Y/N] " response
case $response in 
   [yY][eE][sS]|[yY])
	bash bash/deploy_apache.sh superphy
	;;
   *)
	echo "Not deploying apache."
	;;
esac


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
mkdir blast &> /dev/null || true
if ! find blast/ncbi*/ | read v; then

	if find ../blast/ncbi*/ | read v;
	then
		echo Copying NCBI BLAST+ jar file...
		cp -r ../blast/ncbi*/ blast/
	else
		echo Downloading NCBI BLAST+ jar file...
		cd blast;
	    wget -rc -nd -A "*x64-linux.tar.gz" "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
		tar zxvpf *x64-linux.tar.gz

		cd ..;
		read -r -p "Do you want to backup the NCBI BLAST+ client in ../blast ? [Y/N] " response
		case $response in
		    [yY][eE][sS]|[yY])
			mkdir ../blast || true #Redundant?
			cp -r blast/ncbi*/ ../blast/
			;;*)
		esac
	fi
fi

#Setting up sqlite server for user auth
if ! find data-dev.sqlite | read v; then
	echo Setting up SQL server
	./manage.py db upgrade &> /dev/null
fi

echo Finished
echo """$ bash bash/run.sh""" to run the server
exit 0