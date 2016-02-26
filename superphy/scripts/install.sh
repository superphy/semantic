#!/bin/bash

#Deploy apache
#echo "Deploy apache"
#INDEX=${PWD}/superphy/src/app/mithril
#APACHE=/var/www/html/superphy

#echo $INDEX
#echo $APACHE

#[ -L "$APACHE" ] && sudo rm $APACHE #&& echo "Removing existing symlink at ${APACHE} ..."
#sudo ln  $INDEX $APACHE -sf

#Install sudo packages
if [ "$(whoami)" == "root" ]; then
	echo "Install sudo packages"
	for package in "python-virtualenv" "xvfb" "libyajl2" "wget" "MUMmer" "muscle"; do
		PKG_OK=$(dpkg-query -W --showformat='${Status}' $package|grep "install ok installed")
		if [ "" == "$PKG_OK" ]; then
		  echo "No $package Setting up $package"
		  sudo apt-get --force-yes --yes install $package
		fi
	done
fi
#Initialize venv
echo "Initialize venv"
virtualenv --no-site-packages venv
BASHDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#What does this do?
echo Adding postactivate script to venv/bin/activate
echo "# Superphy Environment setup" >> venv/bin/activate
echo "source $BASHDIR/postactivate.sh" >> venv/bin/activate
source venv/bin/activate
pip install --upgrade pip
pip install -r venv/requirements.txt
nodeenv -p --prebuilt --requirements=venv/npm-requirements.txt
deactivate

if ! find superphy/database/blazegraph/blazegraph-bundled.jar | read v; then
	if find ../blazegraph-bundled.jar | read v; then
		cp ../blazegraph-bundled.jar superphy/database/blazegraph/blazegraph-bundled.jar
	fi
fi

mkdir blast &> /dev/null
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
		mkdir ../blast
		cp -r blast/ncbi*/ ../blast/
	fi
fi

echo ""
echo "Install complete!"