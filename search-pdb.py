import os
import requests
from requests import Request, Session
from hashlib import sha3_256

url = 'https://www.ebi.ac.uk/pdbe/search/pdb/select'

properties = [
    "pdb_id",
    "citation_title",
    "citation_authors",
    "title",
    "experimental_method",
    "entry_authors",
    "pubmed_id",
    "citation_year",
    "journal",
    "organism_scientific_name",
    "assembly_composition",
    "interacting_ligands",
    "tax_id",
    "resolution",
    "status",
    "release_date",
    "prefered_assembly_id",
    "entry_author_list",
    "entry_organism_scientific_name",
    "data_quality",
    "model_quality",
    "experiment_data_available",
    "deposition_date",
    "release_year",
    "molecule_type",
    "pfam_name",
    "uniprot_coverage",
    "compound_id",
    "bound_compound_id",
    "modified_compound_id",
    "uniprot_accession_best",
    "interacting_uniprot_accession",
]

params = {
    "group": "true",
    "group.field": "pdb_id",
    "group.ngroups": "true",
    "json.nl": "map",
    "fl": ",".join(properties),
    "rows": 10000,
    "start": 0,
    # Comment out parameters with special characters
    #"sort":"phmmer(e_value) asc",
    #"xjoin_phmmer.fl": "*",
    #"xjoin_phmmer": "true",    
    # "q":"*:*",
    # "fq":b"{!xjoin}xjoin_phmmer",
    # "group.limit":100,
    # "wt":"json"
}

def get_checksum(seq):
    buffer = seq.encode() + b"\n"
    checksum = sha3_256(buffer).digest().hex()
    return checksum


def detect_homologs(seq):
    myparams = params.copy()
    
    # Does not work, special characters
    #myparams["xjoin_phmmer.external.sequence"] = seq    
    #request = requests.get(url, params=myparams)  

    # Do it the hard way...
    s = Session()
    request = Request('GET', url, params=myparams)
    prepped = request.prepare()
    prepped.url += "&sort=phmmer(e_value)%20asc&xjoin_phmmer.fl=*&xjoin_phmmer=true"
    prepped.url += "&xjoin_phmmer.external.sequence=" + seq
    prepped.url += "&q=*:*&fq={!xjoin}xjoin_phmmer&group.limit=100&wt=json"

    response = s.send(prepped)
    return response.json()


def detect_all_homologs(targets, result_dir):
    result = {}
    for target_name, target in targets.items():
        if target["status"] == "rejected":
            continue
        curr_result = {}
        for entity_name, entity in target.items():
            if not isinstance(entity, dict):
                continue
            if entity["classification"] != "protein":
                continue
            print(entity_name, file=sys.stderr)
            sequence = entity["sequence"]
            checksum = get_checksum(sequence)
            homologs_file = os.path.join(result_dir, checksum + ".json")
            if os.path.exists(homologs_file):
                homologs = json.load(open(homologs_file))
            else: 
                return result ###
                homologs = detect_homologs(sequence)
                json.dump(homologs, open(homologs_file, "w"), sort_keys=False, indent=2)
            curr_result[entity_name] = entity.copy()
            curr_result[entity_name]["homologs"] = homologs_file
        result[target_name] = curr_result
    return result

if __name__ == "__main__":
    import json, sys
    filtered_targets = json.load(open(sys.argv[1]))
    result_dir = sys.argv[2]
    result = detect_all_homologs(filtered_targets, result_dir)
    print(json.dumps(result, sort_keys=True, indent=2))