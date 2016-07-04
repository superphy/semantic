#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
A module used to initialize and set up the superphy namespace in Blazegraph
"""

import gc

from SuperPhy.models.upload import _sparql
from SuperPhy.models.upload.blazegraph_upload import BlazegraphUploader
from SuperPhy.models.upload.blazegraph_setup import BlazegraphSetup

#response = BlazegraphUploader().create_namespace()

def init():
    """
    .
    """
    BlazegraphSetup().setup_curated_data()
    BlazegraphUploader().upload_all_ontologies()
    gc.collect()
    _sparql.delete_blank_nodes()
