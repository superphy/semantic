#!venv/bin/python
"""
  Author: Bryce Drew
"""

import os, sys, subprocess

sys.path.append(os.getcwd()+"/venv/lib/python2.7/site-packages")

import superphy

from superphy.shared import config
config.import_env()

#Debug allows the execution of arbitrary code. Do not use it with production
def run():
    config.start_database()
    os.system("cd superphy/src/app; bash gulp.sh; cd ../../..")
    from superphy.app import create_app
    app = create_app(os.getenv('FLASK_CONFIG') or 'default')
    app.run(host='0.0.0.0', debug=True, use_reloader=False)

def install():
    config.start_database()
    config.install()
    exit()

def upload():
    superphy.upload.foo()
    exit()

def shell():
    import code
    code.interact(local=dict(globals(), **locals()))

def test():
    config.start_database("testing")
    from superphy.upload import main as upload
    upload.init()
    subprocess.call("for f in superphy/src/*; do echo $f;  nosetests $f -vv --exe; done ", shell=True)

options = {
    "install"   :   install,
    "run"       :   run,
    "upload"    :   upload,
    "sparql"    :   superphy.shared.config.start_database,
    "shell"     :   shell,
    "test"      :   test
}

if __name__ == '__main__':
    if (len(sys.argv) >= 2):
        if (sys.argv[1] in options):
            options[sys.argv[1]]()
        else:
            print "Options:"
            for key in options:
                print key
    else:
        #Default
       options['run']()