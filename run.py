"""
This is the file that should be used to build and run superphy.
    Parses std input for program
"""
#pylint: skip-file

import argparse
import subprocess
import os

def parse_args():

    def run(args):
        if not args.no_compile:
            subprocess.call("bash var_www_SuperPhy/SuperPhy/static/compile.sh", shell=True)
        if args.pyserver:
            #Run Python server
            subprocess.call("/usr/bin/python var_www_SuperPhy/run.py", shell=True)
        if not args.no_restart:
            subprocess.call("bash superphy/database/scripts/start.sh", shell=True)
            subprocess.call("sudo /etc/init.d/apache2 reload", shell=True)

    def install(args):
        """
        Download, install, initialize, etc.
        """
        if args.sys:
            args.symlink = True
            #Setup fresh image to use SuperPhy
            #subprocess.call("sudo apt-get update && sudo apt-get upgrade -y", shell=True)
            subprocess.call("sudo apt-get install apache2 curl libapache2-mod-wsgi libyajl2 MUMmer muscle python-dev python-virtualenv wget xvfb -y", shell=True)
            subprocess.call("sudo a2enmod wsgi", shell=True)
            subprocess.call("sudo cp $(pwd)/development_virtualhost.conf /etc/apache2/sites-available/000-default.conf)", shell=True)
            subprocess.call("sudo service apache2 reload", shell=True)

            #This is a very large download. If you aren't going to be uploading data, don't bother downloading it.
            #subprocess.call("mkdir blast", shell=True) 
            #subprocess.call("if ! find blast/ncbi*/ | read v; then cd blast; wget -rc -nd -A "*x64-linux.tar.gz" "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"; tar zxvpf *x64-linux.tar.gz; cd ..; fi", shell=True)

        if args.symlink:
            print "SYMLINKING"
            apache = os.path.join('/var/www/', 'SuperPhy')
            print "APACHE: %s" %apache
            project = os.path.join(os.getcwd(), 'var_www_SuperPhy')
            subprocess.call("sudo rm -f %s" % apache, shell=True)
            subprocess.call("sudo ln -s %s %s" % (project, apache), shell=True)
        
        #Install virtualenv binaries

        subprocess.call("virtualenv --no-site-packages venv", shell=True)
        subprocess.call("venv/bin/pip install --upgrade pip", shell=True)
        subprocess.call("venv/bin/pip install -r venv/requirements.txt", shell=True)
        subprocess.call("venv/bin/nodeenv -p --prebuilt --requirements=venv/npm-requirements.txt", shell=True)
        subprocess.call("git clone --depth=1 https://github.com/eddyystop/mithril-components.git var_www_SuperPhy/SuperPhy/static/js/bower_components/mithril-components;", shell=True)

        #For some reason, even apache doesn't work if you don't run the pyserver beforehand.
        subprocess.call("python run.py run --pyserver", shell=True)

    def pull(args):
        """
        Grab any changes to the git master branch. Pull the development changes
        into production.
        """
        pass

    def test(args):
        pass

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
    install_parser.add_argument('--sys', help='help for b', action="store_true")
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
        print args

    print args.method

    # switch for which method is being used.
    options = {
        'run': run,
        'install': install,
        'pull': pull,
        'test': test
    }

    options[args.method](args)

if __name__ == "__main__":
    parse_args()