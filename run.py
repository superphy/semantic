#!venv/bin/python
"""
  Author: Bryce Drew

  If this file is run with the system binary 'python', then it will take the steps nessesary to run the program with the virtual environment.
  - This includes downloading and initializing the virtual environment, and then running the program again from the virtual environment.
"""
import os, sys, subprocess
import socket, fcntl, struct

def verify_install(verbose = False):
  #If the file not running from venv/bin/python, install venv and then run from venv/bin/python
  if str(sys.argv[0]) != "./%s" % (os.path.basename(__file__)) or ((len(sys.argv) > 1) and sys.argv[1] != "install"):
    """
      Sudo Packages:
      This is for installing the ubuntu packages nessessary to run the system. It CAN run the program as sudo.
      It is currently unclear weather sudo permissions continue to be with the user after it is re-run with the virtual environment
      python.
    """
    if verbose: print "Verifying Install..."
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

    """
      Virtualenv:
      This initializes and downloads the virtualenvironment.
    """

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

    """
      Run with venv:
      This starts a new
    """
    #Run this file using the new venv, and 
    args = ""
    if (len(sys.argv) >= 2):
      args = sys.argv[1]
      for arg in sys.argv[2:]:
        args = args + " " + arg
    os.system("chmod +x %s" % (os.path.basename(__file__)))
    os.system("./%s %s" % (os.path.basename(__file__), args))
    exit(0)
  """
  Deploy Apache
  """
  #src = os.path.dirname(os.path.realpath(__file__))+"/superphy/src/main"
  #dst = os.path.realpath("/var/www/html/superphy")
  #if src != dst

    #print os.path.realpath("/var/www/html/superphy")
    #print os.path.dirname(os.path.realpath(__file__))

def get_ip_address(ifname):
  s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
  return socket.inet_ntoa(fcntl.ioctl(
    s.fileno(),
    0x8915,  # SIOCGIFADDR
    struct.pack('256s', ifname[:15])
  )[20:24])

def import_config(verbose = True):
  #Import Environment Variables
  #   for example: [rdf] url = localhost/blazegraph/namespace/superphy/sparql
  #   becomes:     ['SUPERPHY_RDF_URL', ' localhost/blazegraph/namespace/superphy/sparql']
  if not os.path.exists("config/superphy.cfg"):
    #Copy of config file from default
    try : os.system("cp config/default.cfg config/superphy.cfg")
    except :
      print sys.exc_info()
      sys.exit(0)
  if verbose: print('Setting config to environment...')
  for line in open('config/superphy.cfg'):
    #Ignore # lines
    if line.startswith("#"): continue
    #Seperate and Filter the line into a key-value pair
    env = line\
      .replace("$IP_ADDRESS", get_ip_address('eth0'))\
      .strip()\
      .split('=')
    #Remove appended and prepended whitespace from value
    env[1] = env[1].strip() 
    #Seperate the different headers of our key into a list
    keys = env[0]\
      .strip()\
      .translate(None, '[ ')\
      .split(']')
    #All our env keys start with header 'SUPERPHY'
    env[0] = "SUPERPHY"
    #Aggregate the keys into one named key variable
    for key in keys:
      #Append the 
      env[0] = env[0] + "_" + key.upper()
    #Commit the key value pair as an os environment variable
    if len(env) == 2:
      os.environ[env[0]] = env[1]
      if verbose: print "Environment Variable Uploaded: '" + env[0] + "' = " + env[1]
  #Shows all environment varaibles
  #for key in os.environ.keys():
  #    print "%30s %s" % (key,os.environ[key])

def start_blazegraph(verbose = False):
  #Start blazegraph if we are using a local host
  if get_ip_address("eth0") in os.environ["SUPERPHY_RDF_URL"]:
    if verbose: print "Hosting blazegraph on local server..."
    append = "> /dev/null"
    if verbose: append = ""
    os.system("cd superphy/database; bash scripts/start.sh %s" % append)

def gulp():
  cmd = ["source venv/bin/activate;",
  "",
  "deactivate"]
  subprocess.call(cmd, shell=True)
  os.system("cd superphy/src/main/mithril; ../../../../venv/bin/gulp")

verify_install()
import_config()
start_blazegraph()
gulp()


import time 
while True:
  print "This prints once every 5 seconds"
  time.sleep(5)