#!bin/bash
source scripts/config
if find blazegraph/${JARNAME} | read v; then
	rm blazegraph/${JARNAME}
fi
cd blazegraph
wget "http://sourceforge.net/projects/bigdata/files/latest/download/${JARNAME}";
cd ..
echo Blazegraph has been downloaded!