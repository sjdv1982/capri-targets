def filter_targets(targets):
    result = {}
    for target_name, target in targets.items():
        if target["classification"] != "rejected":
            result[target_name] = target
    return result
    
if __name__ == "__main__":
    import json, sys
    classified_targets = json.load(open(sys.argv[1]))
    result = filter_targets(classified_targets)
    print(json.dumps(result, sort_keys=True, indent=2))
