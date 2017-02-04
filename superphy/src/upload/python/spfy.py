import logging
import time
import os

# Redis Queue
from redis import Redis
from rq import Queue

# other libraries for rdflib
from rdflib import Graph

# our own slightly more general stuff
from insert import insert
from turtle_grapher import generate_output
from turtle_utils import generate_uri as gu, generate_hash

# for various features we add
from savvy import savvy  # serotype/amr/vf

# the only ONE time for global variables
# when naming queues, make sure you actually set a worker to listen to that queue
# we use the high priority queue for things that should be immediately returned to the user
high = Queue('high', connection=Redis())
low = Queue('low', connection=Redis(), default_timeout=600)

def blob_savvy(args_dict):
    '''
    Handles savvy.py's pipeline.
    '''
    # run the much faster vf and serotyping separately of amr
    vf_s = high.enqueue(savvy, dict(args_dict.items() + {'disable_amr': True}))
    amr = low.enqueue(savvy, dict(args_dict.items() + {'disable_vf':True,'disable_serotype':True}))

def monitor():
    '''
    Meant to run until all jobs are finished. Monitors queues and adds completed graphs to Blazegraph.
    '''
    while (high.job_ids not None) or (low.job_ids not None):
        queued_job_ids_high = high.job_ids
        for job_id in queued_job_ids_high:
            job_high = high.fetch_job(job_id_high)
            if job_high.is_finished():
                insert(job_high.result)

        queued_job_ids_low = low.job_ids
        for job_id in queued_job_ids_low:
            job_low = low.fetch_job(job_id_low)
            if job_low.is_finished():
                insert(job_low.result)

def spfyids_single(args_dict):
    from settings import database

    # this is temporary, TODO: include a spqarql query to the db
    uriIsolate = gu(':spfy' + database['count'])

    uriGenome = gu(':' + generate_hash(args_dict['i']))

    args_dict[uris] = {uriIsolate:uriGenome}

    return args_dict

def spfyids_directory(args_dict):
    '''
    TODO: make the database count actually work
    This is meant to preallocate spfyIDs
    -note may have problems with files that fail (gaps in id range)
    '''
    from settings import database
    files = os.listdir(args_dict['i'])
    count = database['count']
    uris = {}
    for f in files:
        uris[gu(':spfy' + count]=gu(':' +generate_hash(f))
    count=count + len(files)
    args_dict['uris'] = uris
    #TODO: write-out

    return args_dict

def spfy(args_dict):
    '''
    '''
    # check if a directory was passed or a just a single file
    # updates args_dict with appropriate rdflib.URIRef's
    if os.path.isdir(args_dict['i']):
        args_dict = spfyids_directory(args_dict)
    else:
        args_dict=spfyids_single(args_dict)

    print 'Starting savvy call'
    logging.info('Starting savvy call...')
    sav = high.enqueue(savvy, args_dict)
    logging.info(sav.id)
    time.sleep(180)
    graph = sav.result

    logging.info('uploading to blazegraph')
    print "Uploading to Blazegraph"
    print insert(graph)
    print 'uploaded wooot!'


if __name__ == "__main__":
    import argparse

    from ConfigParser import SafeConfigParser

    # parsing cli-input
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        help="FASTA file or directory",
        required=True
    )
    parser.add_argument(
        "--disable-serotype",
        help="Disables use of the Serotyper. Serotyper is triggered by default.",
        action="store_true"
    )
    parser.add_argument(
        "--disable-vf",
        help="Disables use of ECTyper to get associated Virulence Factors. VFs are computed by default.",
        action="store_true"
    )
    parser.add_argument(
        "--disable-amr",
        help="Disables use of RGI to get Antimicrobial Resistance Factors.  AMR genes are computed by default.",
        action="store_true"
    )

    args = parser.parse_args()
    # we make a dictionary from the cli-inputs and add are uris to it
    # mainly used for when a func needs a lot of the args
    args_dict = vars(args)

    # starting logging
    #TODO: move this to global and see it if breaks
    logging.basicConfig(
        filename='outputs/spfy' + __name__ +
        args_dict['i'] + '.log',
        level=logging.INFO
    )

    spfy(args_dict)
