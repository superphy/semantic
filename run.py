#!venv/bin/python
"""
  Author: Bryce Drew
"""

import os
import sys
import subprocess

sys.path.append(os.getcwd()+"/venv/lib/python2.7/site-packages")

import superphy

from superphy.shared import config
config.import_env()

#Debug allows the execution of arbitrary code. Do not use it with production
def run():
    """
    .
    """
    config.start_database()
    os.system("cd var_www_SuperPhy/SuperPhy/static; bash compile.sh; cd ../../..")
    os.system("sudo /etc/init.d/apache2 reload")
    #from superphy.app import create_app
    #app = create_app(os.getenv('FLASK_CONFIG') or 'default')
    #app.run(host='0.0.0.0', debug=True, use_reloader=False)

def install():
    """
    .
    """
    config.start_database()
    config.install()
    exit()

def uploader():
    """
    .
    """
    superphy.upload.foo()
    exit()

def shell():
    """
    .
    """
    import code
    code.interact(local=dict(globals(), **locals()))

def test():
    """
    .
    """
    config.start_database("testing")
    from superphy.upload import main as upload
    upload.init()
    subprocess.call(
        "for f in superphy/src/*; do echo $f; nosetests $f -vv --exe; done",
        shell=True
    )

OPTIONS = {
    "install"   :   install,
    "run"       :   run,
    "upload"    :   uploader,
    "sparql"    :   superphy.shared.config.start_database,
    "shell"     :   shell,
    "test"      :   test
}

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        if sys.argv[1] in OPTIONS:
            OPTIONS[sys.argv[1]]()
        else:
            print "OPTIONS:"
            for key in OPTIONS:
                print key
    else:
        #Default
        OPTIONS['run']()
