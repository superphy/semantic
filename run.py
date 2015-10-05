#!superphy/bin/python
import subprocess
from app import app
subprocess.call(['./init.sh'])
subprocess.call(['./startdb.sh'])
app.run(debug=True)