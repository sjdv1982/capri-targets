import json

def analyze_all_searches(targets):
    result = {}
    for target_name, target in targets.items():
        result[target_name] = {}
        for entity_name, entity in target.items():            
            homologs_file = entity["homologs"]
            homologs = json.load(open(homologs_file))
            matches = []
            for match in homologs["xjoin_phmmer"]["external"]:
                mymatch = {}
                m = match["doc"]
                if m["e_value"] > 0.05:
                    continue
                mymatch["pdb"] = match["joinId"]
                mymatch["pdb_chain"] = m["target"]
                mymatch["query_region"] = m["query_sequence_start"], m["query_sequence_end"]
                mymatch["target_region"] = m["target_sequence_start"], m["target_sequence_end"]
                mymatch["seqid"] = m["identity_percent"]
                matches.append(json.dumps(mymatch))
            result[target_name][entity_name] = entity.copy()
            result[target_name][entity_name]["matches"] = matches
    return result


if __name__ == "__main__":
    import json, sys
    searched_pdb = json.load(open(sys.argv[1]))
    result = analyze_all_searches(searched_pdb)
    print(json.dumps(result, sort_keys=True, indent=2))
