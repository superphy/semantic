#!/bin/bash
#This removes the old symlink from apache,
#and adds a new one; the new one points to the working 
#index.html that is in the folder that the script is run on.

INDEX_LOCATION="/app/static/"

index="$(pwd)"${INDEX_LOCATION}
apache="/var/www/html/superphy"

[ -L $apache ] && sudo rm $apache && echo "Removing existing symlink at ${apache} ..."
echo Making symlink from ${index}, to ${apache} ...

sudo ln  $index $apache -s

echo Apache Deployed!