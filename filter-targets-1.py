def filter_targets(classified_entities):
    result = classified_entities.copy()
    impossible = (
        "nucleic acid",
        "unnatural peptide",
        "unnatural protein",
        "large peptide"
    )
    for target_name, target in classified_entities.items():
        nprot = 0
        npep = 0
        rejection = None
        for entity_name, entity in target.items():
            for impos in reversed(impossible):
                classification = entity["classification"]
                if classification == impos:
                    rejection = "Target contains " + classification
                    break
            if rejection is not None:
                break

        if rejection is None:
            for entity_name, entity in target.items():
                classification = entity["classification"]
                if classification == "cofactor":
                    pass
                elif classification == "protein":
                    nprot += 1
                elif classification == "small peptide":
                    npep += 1
            else:
                if nprot < 2:
                    if npep == 0:
                        rejection = "Only one protein partner"
                    else:
                        rejection = "Small peptide-protein complex"
        if rejection is None:
            result[target_name]["status"] = "potential"
        else:    
            result[target_name]["status"] = "rejected"
            result[target_name]["rejection"] = rejection
    return result


if __name__ == "__main__":
    import json, sys
    classified_entities = json.load(open(sys.argv[1]))
    result = filter_targets(classified_entities)
    print(json.dumps(result, sort_keys=True, indent=2))
