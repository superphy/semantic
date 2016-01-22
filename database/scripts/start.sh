#!bin/bash
source scripts/config

#!/bin/bash
#starts up Blazegraph

bash scripts/stop.sh <&- 1>/dev/null

echo Starting Blazegraph...


cd blazegraph
if find ${JARNAME} | read v; then
	/usr/bin/java -ea -Xmx4g -server -cp ${JARNAME} com.bigdata.rdf.sail.webapp.NanoSparqlServer $PORT $NAMESPACE $PROPERTIES_FILE 2>&1 &
	response=$(curl -s -i http://localhost:${PORT}/blazegraph/ | grep "200 OK")
	while [[ $response = "" ]]; do
		echo "Blazegraph not initialized"
		sleep 1
		response=$(curl -s -i http://localhost:${PORT}/blazegraph/ | grep "200 OK")
	done
	echo Blazegraph has been started!
else
	echo Cannot find ${JARNAME}. Use 'Init.sh' to download jar file.
fi
cd ..
exit 0