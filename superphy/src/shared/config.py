#!/usr/bin/python

import os, subprocess

def install():
    os.system("bash superphy/scripts/install.sh")
    apache = os.path.join('/var/www/html', 'superphy')
    index = os.path.join(os.getcwd(), "superphy/src/app/mithril")
    if os.path.realpath(apache) != index:
        subprocess.call("sudo rm '%s'; sudo ln -s '%s' '%s'" % (apache, index, apache), shell=True)
        print "Apache is now deployed!"
    else:
      print "Apache is already deployed!"

#This function belongs in superphy/... package
def start_database(default="development"):
    os.system("cd superphy/database; bash scripts/start.sh %s" % default)

#This function belongs in superphy/... package
def import_env(verbose = True):
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


#This function belongs in superphy/... package
def get_ip_address(ifname):
    import socket, fcntl, struct
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    return socket.inet_ntoa(fcntl.ioctl(
        s.fileno(),
        0x8915,  # SIOCGIFADDR
        struct.pack('256s', ifname[:15])
    )[20:24])
