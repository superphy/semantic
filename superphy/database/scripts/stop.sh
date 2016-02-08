#!bin/bash
source scripts/config
result=$(ps ax | grep ${PORT} | grep java | awk '{print $1}')
while [ ! -z $result ]; do
	result=$(ps ax | grep ${PORT} | grep java | awk '{print $1}')
	if [ ! -z "$result" ]; then
		{
			kill -9 $result
		} || {
			echo "Sudo killing existing blazegraph..."
			sudo kill -9 $result
		}
		echo "Blazegraph has been terminated!"
	else
		echo "No Blazeraph running on port $PORT..."
	fi
done