#!/usr/bin/env python
# -*- coding: UTF-8 -*-

if __name__ == "__main__":
    import subprocess
    from os import listdir

    for filename in listdir('tmp'):
        print 'working on ' + filename
        try:
            subprocess.call('python insert.py -i tmp/' + filename)
        except:
            print 'call failed'
