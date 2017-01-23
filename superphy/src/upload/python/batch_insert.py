#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import subprocess

from multiprocessing import Pool, cpu_count
from os import listdir

def batch_call(filename):
    subprocess.call(['./insert.py', '-i', 'tmp/' + filename])
if __name__ == "__main__":

    p = Pool(cpu_count())
    p.map(batch_call, listdir('tmp'))

    print '***ALL DONE***'
