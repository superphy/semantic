"""
This is the file that should be used to build and run superphy.
    Parses std input for program

    TODO: verbose is only showing debug statments. We need to tie it in with
        the subprocesses.
"""
# pylint: disable=line-too-long

import argparse
import subprocess
import os

def run(args):
    """
    Compile and run the code
    """
    if not args.no_compile:
        subprocess.call("bash app/SuperPhy/static/compile.sh", shell=True)
    if not args.no_restart:
        subprocess.call("bash superphy/database/scripts/start.sh", shell=True)
        subprocess.call("sudo /etc/init.d/apache2 reload", shell=True)
    if args.pyserver:
        #Run Python server
        #Remember you don't want this running all the time. This is a
        #   security hazard if you allow port 5000 traffic.
        subprocess.call("/usr/bin/python app/run.py", shell=True)



GIT_CONF_FILE = "development_virtualhost.conf"
APACHE = "apache2"
APACHE_CONF_DESTINATION = "/etc/{apache}/sites-available".format(apache=APACHE)
APACHE_SYMLINK_DESTINATION = '/var/www/SuperPhy'
NETWORK_FACEING_SUBDIRECTORY = 'app'
VIRTUALENV_FOLDER = "venv"
NPM_REQUIRMENTS = "npm-requirements.txt"
BOWER_COMPONENTS_DESTINATION = "app/SuperPhy/static/js/bower_components/mithril-components"

def install(args):
    """
    Download, install, initialize, etc.

    You should already have system packages installed on your build.
    """
    if args.sys:
        args.symlink = True
        #Connect this repo with apache
        subprocess.call("sudo a2enmod wsgi", shell=True)
        subprocess.call("sudo cp $(pwd)/{conf} {conf_destination}/000-default.conf".format(conf=GIT_CONF_FILE, conf_destination=APACHE_CONF_DESTINATION), shell=True)
        subprocess.call("sudo service {apache} reload".format(apache=APACHE), shell=True)

    #This is a very large download. If you aren't going to be uploading data, don't bother downloading it.
    #subprocess.call("mkdir blast", shell=True)
    #subprocess.call("if ! find blast/ncbi*/ | read v; then cd blast; wget -rc -nd -A "*x64-linux.tar.gz" "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"; tar zxvpf *x64-linux.tar.gz; cd ..; fi", shell=True

    if args.symlink:
        if args.verbose:
            print "SYMLINKING"
        subprocess.call("sudo rm -f {0}".format(APACHE_SYMLINK_DESTINATION), shell=True)
        subprocess.call("sudo ln -s {0} {1}".format(os.path.join(os.getcwd(), NETWORK_FACEING_SUBDIRECTORY), APACHE_SYMLINK_DESTINATION), shell=True)

    #Install virtualenv binaries

    if args.verbose:
        print "Initializing virtualenv"
    subprocess.call("virtualenv --no-site-packages {venv}".format(venv=VIRTUALENV_FOLDER), shell=True)

    if args.verbose:
        print "Upgrading virtualenv pip"
    subprocess.call("{venv}/bin/pip install --upgrade pip".format(venv=VIRTUALENV_FOLDER), shell=True)

    if args.verbose:
        print "Installing Pip requirments into {venv}".format(venv=VIRTUALENV_FOLDER)
    subprocess.call("{venv}/bin/pip install -r venv/requirements.txt".format(venv=VIRTUALENV_FOLDER), shell=True)

    if args.verbose:
        print "Installing Npm requirments into virtualenv"
    subprocess.call("cd {venv}/lib; ../bin/nodeenv -p --prebuilt --requirements=../{npm}".format(npm=NPM_REQUIRMENTS, venv=VIRTUALENV_FOLDER), shell=True)


    #Refactor This into a submodule of our project?
    if args.verbose:
        print "Grabbing mithril components from github"
    subprocess.call("git clone --depth=1 https://github.com/eddyystop/mithril-components.git {dest};".format(dest=BOWER_COMPONENTS_DESTINATION), shell=True)

    #For some reason, even apache doesn't work if you don't run the pyserver beforehand.
    if args.verbose:
        print "Restarting and running pyserver. (Press Ctrl-c to exit):"

    #subprocess.call("python {this} run --pyserver {verbose}".format(this=__file__, verbose="--verbose" if args.verbose else ""), shell=True)


def pull(args):
    """
    Grab any changes to the git master branch. Pull the development changes
    into production.
    """
    pass

def test(args):
    """
    Test the code. (Need to be in the venv context to run this)
    """
    #Start up a new Blazegraph instance, but make the DB files reference non-production.
    subprocess.call("bash superphy/database/scripts/start.sh testing", shell=True)

    subprocess.call("venv/bin/nosetests -vv --exe", shell=True)
    

def parse_args():
    """
    Parse the command-line arguments
    """
    #Define the acceptable arguements for the program:

    # create the top-level parser
    parser = argparse.ArgumentParser(prog='SuperPhy')
    subparsers = parser.add_subparsers(help='help for subcommand', dest='method')

    # Run subparser
    run_parser = subparsers.add_parser('run', help='Compile, and restart servers.')
    run_parser.add_argument("-p", "--pyserver", help="Run the python development server.", action="store_true")
    run_parser.add_argument("--no-compile", help="Don't recompile code", action="store_true")
    run_parser.add_argument("--no-restart", help="Don't restart servers", action="store_true")
    run_parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")

    # Install subparser
    install_parser = subparsers.add_parser('install', help='Download, install, initialize, etc.')
    install_parser.add_argument('--sys', help='Install system packages. (Includes --symlink)', action="store_true")
    install_parser.add_argument('--symlink', help='Point a symlink from apache to SuperPhy', action="store_true")
    install_parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")

    # Pull subparser
    pull_parser = subparsers.add_parser('pull', help='Pull the development changes into production.')
    pull_parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")

    # Test subparser
    test_parser = subparsers.add_parser('test', help='Test the code')
    test_parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")

    #Parse the arguements defined above:

    args = parser.parse_args()

    if args.verbose:
        print args.method

    return args

def import_env(verbose=True):
    """
    Import Environment Variables from configuration file.
     for example: [rdf] url = localhost/blazegraph/namespace/superphy/sparql
     becomes:     ['SUPERPHY_RDF_URL', ' localhost/blazegraph/namespace/superphy/sparql']
    """

    DEFAULT_CONFIG_DIRECTORY = "config/default.cfg"
    ACTIVE_CONFIG_DIRECTORY = "config/active.cfg"

    #Create active configuration file.
    if not os.path.exists("{active}".format(active=ACTIVE_CONFIG_DIRECTORY)):
        #Copy of config file from default
        subprocess.call("cp {default} {active}".format(active=ACTIVE_CONFIG_DIRECTORY, default=DEFAULT_CONFIG_DIRECTORY), shell=True)
    #Setting the Environment Variables
    for line in open("{active}".format(active=ACTIVE_CONFIG_DIRECTORY)):
        #Ignore lines that are commented out.
        if line.startswith("#"):
            continue
        #Seperate and Filter the line into a key-value pair
        key, value = line.strip().split('=')

        #Sort the key into each tag
        key = key.strip().translate(None, '[ ').split(']')

        #Add the header 'SUPERPHY' to the front of all keys
        if key[0] is not "SUPERPHY":
            key.insert(0, "SUPERPHY")

        #Aggregate the keys into one string.
        for i, tag in enumerate(key):
            key[i] = tag.upper()
        key = '_'.join(key)

        #Remove Whitespace from value
        value = value.strip()

        #Set the environment variables.
        os.environ[key] = value

    #Show all environment varaibles
    if verbose:
        for key in os.environ.keys():
            if "SUPERPHY" in key:
                print key, os.environ[key]




if __name__ == "__main__":

    OPTIONS = {
        'run': run,
        'install': install,
        'pull': pull,
        'test': test
    }

    ARGS = parse_args()
    import_env(ARGS.verbose)
    OPTIONS[ARGS.method](ARGS)
