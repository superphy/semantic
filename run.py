#!superphy/bin/python
import subprocess
from app import app
subprocess.call(['./bash/init.sh'])
subprocess.call(['./bash/startdb.sh'])
app.run(debug=True)