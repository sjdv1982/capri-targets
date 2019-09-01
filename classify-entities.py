import re

non_aa = re.compile('\([^)]*\)')

def adjust_sequence(sequence_original):
    """Adjusts a sequence, replacing all non-amino acids
    Input: a sequence where all non-amino acids are between ()
    Example: ACD(MSE)DQD(HYP)D
    
    Output: the sequence with the following modifications:
    - (MSE) replaced by M
    - (XXX) replace by H, where XXX is  one of HIE, HIP, HSD, HSE
    - All other (...) replaced by X
    """
    repl = {
        "MSE": "M",
        "HIE": "H",
        "HIP": "H",
        "HSD": "H",
        "HSE": "H"
    }
    sequence = sequence_original
    for source, target in repl.items():
        sequence = sequence.replace("(" + source + ")", target)
    sequence = non_aa.sub("X", sequence)
    return sequence


def classify_entity(entity):
    result = {}
    sequence_original = entity["molecule_sequence"]
    result["sequence_original"] = sequence_original
    sequence = adjust_sequence(sequence_original)
    result["sequence"] = sequence
    if len(sequence) == 0:
        return
    nmod = sequence.count("X")
    pct_mod = nmod / len(sequence) * 100
    if "U" in sequence:
        pct_mod = 100
    if len(sequence) < 6:
        if pct_mod > 50:
            classification = "cofactor"
        elif pct_mod > 20:
            classification = "unnatural peptide"
        else:
            classification = "small peptide"
    elif len(sequence) <= 10:
        if pct_mod >= 50:
            classification = "nucleic acid"
        elif pct_mod >= 20:
            classification = "unnatural peptide"
        else:
            classification = "small peptide"
    elif len(sequence) <= 30:
        if pct_mod >= 50:
            classification = "nucleic acid"
        elif pct_mod >= 20:
            classification = "unnatural peptide"
        else:
            classification = "large peptide"
    else:
        if pct_mod >= 50:
            classification = "nucleic acid"
        elif pct_mod >= 10:
            classification = "unnatural protein"
        else:
            classification = "protein"
    result["classification"] = classification
    return result

def classify_on_hold_targets(on_hold_targets):
    result = {}
    for target in on_hold_targets["grouped"]["pdb_id"]["groups"]:
        name = target["groupValue"]
        entities = {}
        ok = True
        for entity in target["doclist"]["docs"]:
            entity_result = classify_entity(entity)
            entities[entity["entry_entity"]] = entity_result
        result[name] = entities
    return result

if __name__ == "__main__":
    import json, sys
    on_hold_targets = json.load(open(sys.argv[1]))
    result = classify_on_hold_targets(on_hold_targets)
    print(json.dumps(result, sort_keys=True, indent=2))
