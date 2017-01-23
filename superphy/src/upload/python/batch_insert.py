#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import subprocess

from multiprocessing import Pool, cpu_count
from os import listdir
from time import time

def batch_call(filename):
    subprocess.call(['./insert.py', '-i', 'tmp/' + filename])
if __name__ == "__main__":

    start = time()
    print 'Starting batch insert at: ', start

    p = Pool(cpu_count())
    p.map(batch_call, listdir('tmp'))

    print '***ALL DONE***'
    print 'Completed at: ', time()
    print 'Elapsed: ', time() - s
