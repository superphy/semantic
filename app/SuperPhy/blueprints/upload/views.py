"""
.
"""

from flask import jsonify

from SuperPhy.blueprints.upload import upload

from SuperPhy.models.upload.metadata_upload import MetadataUploader, Metadata, GenomeMetadataUploader, GenomeMetadata, GeneMetadataUploader, GeneMetadata
from SuperPhy.models.upload._utils import generate_path

@upload.route('/', methods=['GET', 'POST'])
def upload_file():
    """
    We also want to be changing the filename somehow, and keep track of who
    is uploading what
    """
    case = MetadataUploader('samples/2_genome.json')
    metadata = Metadata("JHNI00000000")



    return jsonify({})
