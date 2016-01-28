#!venv/bin/python
"""
	Author: Bryce Drew

	If this file is run with the system binary 'python', then it will take the steps nessesary to run the program with the virtual environment.
	- This includes downloading and initializing the virtual environment, and then running the program again from the virtual environment.
"""
import os, sys, subprocess

def verify_install(verbose = True):
	#If the file not running from venv/bin/python, install venv and then run from venv/bin/python
	if str(sys.argv[0]) != "./%s" % (os.path.basename(__file__)):
		sudo_install(verbose=False)
		#Check if venv binaries are installed
		if not (os.path.exists(os.path.dirname(os.path.abspath(__file__)) + "/venv/bin")):
			if verbose: print "Setting up Virtual Environment"
			cmd = ["virtualenv --no-site-packages venv;",
				"source venv/bin/activate;",
				"pip install --upgrade pip",
				"pip install -r venv/requirements.txt",
				"nodeenv -p --prebuilt --requirements=venv/npm-requirements.txt",
				"deactivate"]
			subprocess.call(cmd, shell=True)
		#Run this file using the new venv, and 
		args = ""
		if (len(sys.argv) >= 2):
			args = sys.argv[1]
			for arg in sys.argv[2:]:
				args = args + " " + arg
		os.system("chmod +x %s" % (os.path.basename(__file__)))
		os.system("./%s %s" % (os.path.basename(__file__), args))
		exit(0)

def sudo_install(verbose = False):
	"""
		This function is for installing the ubuntu packages nessessary to run the system. It CAN run the program as sudo.
		It is currently unclear weather sudo permissions continue to be with the user after it is re-run with the virtual environment
		python.
	"""
	import apt #This is not part of the venv, so it won't be imported unless we are running python from the the system binaries.
	packages = ["python-virtualenv", "xvfb", "libyajl2"]

	for package in packages:
		cache = apt.cache.Cache()
		if cache[package].is_installed:
			if verbose: print "%s is already installed!" % package
		else:
			if os.geteuid() != 0:
				if verbose: print "In order to install nessesary packages, %s has been given sudo privleadges." % os.path.basename(__file__)
				args = ['sudo', sys.executable] + sys.argv + [os.environ]
				os.execlpe('sudo', *args)
			package = cache[package]
			package.mark_install()
			if verbose: print "Installing %s..." % package
			cache.commit()
	if verbose: print "Success! Ubuntu packages are installed!"


verify_install()

import time 
while True:
    print "This prints once every 5 seconds"
    time.sleep(5)