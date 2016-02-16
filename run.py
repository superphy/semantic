#!venv/bin/python
"""
  Author: Bryce Drew
"""

import os, sys

try:
    import superphy
except:
    if not os.path.isfile(os.getcwd()+"venv/bin/python"):
        sys.path.append(os.getcwd()+"/venv/lib/python2.7/site-packages")
        import superphy
        superphy.config.install()
        exit()
    else:
        exit()

sys.path.append(os.getcwd()+"/superphy/src/main")

#Debug allows the execution of arbitrary code. Do not use it with production
def run():
    superphy.config.import_env()
    #superphy.config.start_database()
    os.system("cd superphy/src/main; bash gulp.sh; cd ../../..")
    from app import create_app
    app = create_app(os.getenv('FLASK_CONFIG') or 'default')
    app.run(host='0.0.0.0', debug=True, use_reloader=False)

def install():
    superphy.config.install()
    superphy.config.start_database()
    exit()

def upload():
  #superphy.config.import_env()
  #superphy.config.start_database()
    superphy.upload.foo()
    exit()

def shell():
    import code
    code.interact(local=dict(globals(), **locals()))

options = {
    "install"   :   install,
    "run"       :   run,
    "upload"    :   upload,
    "sparql" : superphy.config.start_database,
    "shell"     :   shell
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