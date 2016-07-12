from flask import jsonify, request
from SuperPhy.models import sparql

from SuperPhy.blueprints.upload import upload

def validate_parameters(key, value):
    default = "ERROR: NO ERROR CHECKER YET. VALUE: %s" % value
    switcher = {
        "Asymptomatic": default,
        "Bacteremia": default,
        "Bacteriuria": default,
        "Bloody diarrhea": default,
        "Crohn's Disease": default,
        "Diarrhea": default,
        "Gastroenteritis": default,
        "Hemolytic-uremic syndrome": default,
        "Hemorrhagic colitis": default,
        "Mastitis": default,
        "Meningitis": default,
        "Necrotizing fasciitis": default,
        "Omphalitis": default,
        "Peritonitis": default,
        "Pneumonia": default,
        "Pyelonephritis": default,
        "Septicaemia": default,
        "Ulcerateive colitis": default,
        "Urinary tract infection (cystitis)": default,
        "additional_notes": default,
        "city": default,
        "data_of_isolation": default,
        "external_db_identifier_accession_id": default,
        "external_db_identifier_database": default,
        "external_db_identifier_version_number": default,
        "file": default,
        "genome_alias": default,
        "genome_description": default,
        "genome_name": default,
        "genome_sequence_status": default,
        "group_label": default,
        "isolation_host": default,
        "isolation_source": default,
        "keywords": default,
        "molecular_type": default,
        "owner": default,
        "privacy_settings": default,
        "province/state": default,
        "pubmed_id": default,
        "serotype": default,
        "strain": default,
        "swollen head syndrome": default
    }
    if key in switcher:
        return switcher[key]
    else:
        return default

@upload.route('/genome', methods=['POST', 'GET'])
def genome():
    """
    Endpoint for posting genomes.

    Todo:
    add security
    do something with the data
    validate the data
    """

    #Fix items that can't be posted due to special characters
    value = request.json.get("Urinary tract infection")
    if value != '':
        request.json["Urinary tract infection (cystitis)"] = value

    value = request.json.get("Crohns Disease")
    if value != '':
        request.json["Crohn's Disease"] = value

    #Get parameters from request object
    params = dict((key, "") for key in [
        'file', 'genome_name', 'strain', 'serotype', 'isolation_host',
        'isolation_source', 'data_of_isolation',
        'province/state', 'city', 'privacy_settings',
        'group_label', 'additional_notes', 'genome_description',
        'keywords', 'owner', 'genome_alias', 'external_db_identifier_database',
        'external_db_identifier_accession_id',
        'external_db_identifier_version_number', 'pubmed_id',
        'genome_sequence_status', 'molecular_type'
    ])
    params.update(dict((key, "") for key in sparql.get_all_syndromes()))

    errors = {'status': 'error'}
    for key in params:
        params[key] = request.json.get(key)
        error = validate_parameters(key, params[key])
        if error:
            errors[key] = error

    #Show all errors
    if errors:
        return jsonify(errors)

    #Add data to uploader processes

    #Return success notification
    return jsonify({'status':'ok', 'params': params})

'''
curl -i -X POST -H "Content-Type: application/json" -d '{ "Asymptomatic": "", "Bacteremia": "", "Bacteriuria": "", "Bloody diarrhea": "", "Crohn's Disease": "", "Diarrhea": "", "Gastroenteritis": "", "Hemolytic-uremic syndrome": "", "Hemorrhagic colitis": "", "Mastitis": "", "Meningitis": "", "Necrotizing fasciitis": "", "Omphalitis": "", "Peritonitis": "", "Pneumonia": "", "Pyelonephritis": "", "Septicaemia": "", "Ulcerateive colitis": "", "Urinary tract infection (cystitis)": "", "additional_notes": "", "city": "", "data_of_isolation": "", "external_db_identifier_accession_id": "", "external_db_identifier_database": "", "external_db_identifier_version_number": "", "file": "", "genome_alias": "", "genome_description": "", "genome_name": "", "genome_sequence_status": "", "group_label": "", "isolation_host": "", "isolation_source": "", "keywords": "", "molecular_type": "", "owner": "", "privacy_settings": "", "province/state": "", "pubmed_id": "", "serotype": "", "strain": "", "swollen head syndrome": "" }' http://127.0.0.1:5000/upload/genome
'''
