import sys
import json
import copy
import numpy as np

class Template:
    # contains all the matches from a single PDB
    def __init__(self, pdb, seqlens):
        self.pdb = pdb
        if seqlens is not None:
            mask = {}
            for entity, seqlen in seqlens.items():
                mask[entity] = np.zeros(seqlen, bool)
            self.mask = mask

    def add_region(self, entity, region):
        self.mask[entity][region[0]-1:region[1]] = 1

    def is_trivial(self):
        """Assesses if template covers at least 20 residues in every entity"""
        for entity in self.mask:
            m = self.mask[entity]
            if not contiguous_region(m, 20, "any"):
                print("Not trivial", self.pdb, entity, np.where(np.diff(m)), m)
                return False
        return True

    def redundant(self, other):
        """Assesses if template is redundant with another
        -1: other fully contains this
        0: non-redundant
        1: this fully contains other"""
        m0 = self.mask
        m1 = other.mask
        delta1, delta2 = False, False
        for entity in m0:
            mm0, mm1 = m0[entity], m1[entity]
            assert len(mm0) == len(mm1)
            d = mm0.astype(int) - mm1.astype(int)
            curr_delta1 = np.any(d>0)
            curr_delta2 = np.any(d<0)
            if curr_delta1:
                delta1 = True
            if curr_delta2:
                delta2 = True
        if not delta2:
            return 1
        elif not delta1:
            return -1
        else:
            return 0

def contiguous_region(mask, minlength, where):
    assert where in ("middle", "termini", "any")
    emask = np.zeros(len(mask)+2, int)
    emask[1:-1] = mask
    dif = np.diff(emask)
    starts = np.where(dif==1)[0]
    ends = np.where(dif==-1)[0] - 1
    for start, end in zip(starts,ends):
        if end - start + 1 < minlength:
            continue
        terminus = (start == 0 or end == len(mask)-1)
        if where == "middle" and terminus:
            continue
        if where == "termini" and not terminus:
            continue
        return True
    return False

def remove_small_regions(mask, minlength):
    result = copy.deepcopy(mask)
    emask = np.zeros(len(mask)+2, int)
    emask[1:-1] = mask
    dif = np.diff(emask)
    starts = np.where(dif==1)[0]
    ends = np.where(dif==-1)[0] - 1
    for start, end in zip(starts,ends):
        if end - start + 1 < minlength:
            result[start:end+1] = 0
    return result

def detect_coverage(template_groups):
    mask = copy.deepcopy(template_groups[0].mask)
    for group in template_groups:
        for entity in mask:
            mask[entity] |= group.mask[entity]
    result = {}
    for entity in mask:
        mask2 = remove_small_regions(mask[entity], 20)
        if not mask2.sum():
            result[entity] = "bad"
            continue
        inv_mask = ~mask2
        big_termini_gaps = contiguous_region(inv_mask, 51, "termini")
        middle_gaps = contiguous_region(inv_mask, 1, "middle")
        big_middle_gaps = contiguous_region(inv_mask, 21, "middle")
        if big_middle_gaps:
            result[entity] = "bad"
        elif big_termini_gaps:
            if middle_gaps:
                result[entity] = "multi-domain"
            else:
                result[entity] = "contiguous"
        else:
            if middle_gaps:
                result[entity] = "multi-domain"
            else:
                result[entity] = "complete"
    return result


def get_templates_and_groups(target, seqid_threshold):
    templates = {}
    seqlens = {}
    for entity_name, entity in target.items(): 
        seqlens[entity_name] = len(entity["sequence"])
    for entity_name, entity in target.items(): 
        for match_json in entity["matches"]:
            match = json.loads(match_json)
            pdb = match["pdb"]
            if match["seqid"] < seqid_threshold:
                continue
            if pdb not in templates:
                templates[pdb] = Template(pdb, seqlens)
            
            templates[pdb].add_region(entity_name, match["query_region"])
    
    for pdb1, t1 in list(templates.items()):
        if pdb1 not in templates:
            continue
        for pdb2, t2 in list(templates.items()):
            if t1 is t2:
                continue
            if pdb2 not in templates:
                continue
            redundant = t1.redundant(t2)
            if redundant == 1:
                #print("REDUNDANT1:", pdb1, "over", pdb2)
                templates.pop(pdb2)
            elif redundant == -1:
                #print("REDUNDANT2:", pdb2, "over", pdb1)
                templates.pop(pdb1)
                break        

    template_pdbs = list(templates.keys())
    if not len(templates):
        return [], []
    first_template = templates[template_pdbs[0]]
    group1 = Template([first_template.pdb], None)
    group1.mask = copy.deepcopy(first_template.mask)
    template_groups = [group1]
    for t_pdb in template_pdbs[1:]:
        t = templates[t_pdb]
        for g in template_groups:
            for entity in target:
                m1 = t.mask[entity]
                m2 = g.mask[entity]
                combined = m1 & m2
                if not contiguous_region(combined, 20, where="any"):
                    break
            else:
                g.pdb.append(t_pdb)
                for entity in target:
                    m1 = t.mask[entity]
                    m2 = g.mask[entity]
                    m1 |= m2
                break
    
    return templates, template_groups

def filter_all_targets(targets):
    result = {}    
    for target_name, target in targets.items():
        print(target_name, file=sys.stderr)
        templates, template_groups = get_templates_and_groups(target, seqid_threshold=0)
        coverage = detect_coverage(template_groups)
        result["coverage"] = coverage
        any_bad = any(["bad" in coverage.values()])
        any_multi_domain = any(["multi_domain" in coverage.values()])
        close_templates, close_template_groups = get_templates_and_groups(target, seqid_threshold=40)
        close_coverage = detect_coverage(close_template_groups)
        print("CLOSE COV", close_coverage)
        close_template_coverage = not any(["bad" in close_coverage.values()]) \
            and not any(["multi_domain" in close_coverage.values()])
        
        my_result = {"entities": list(target.keys())}
        my_result["classification"] = None
        if any_bad:
            my_result["classification"] = "rejected"
            my_result["rejection"] = "Unbound coverage is bad for one or more entities"
        elif any_multi_domain:                
            my_result["classification"] = "rejected"
            my_result["rejection"] = "Unbound coverage is multi-domain for one or more entities"
        else:
            trivial_pdbs = []
            trivial_group = None            
            for pdb, t in list(close_templates.items()):
                if t.is_trivial():
                    trivial_pdbs.append(pdb)
            if len(close_template_groups) == 1 and close_template_coverage:
                assert close_template_groups[0].is_trivial()                
            for group in close_template_groups:
                if group.is_trivial():
                    trivial_group = group
                    break

            if len(trivial_pdbs):
                my_result["classification"] = "trivial"
                if len(trivial_pdbs) == 1:
                    my_result["rejection"] = "Trivial template: " + trivial_pdbs[0]
                else:
                    my_result["rejection"] = "Trivial templates: " + str(trivial_pdbs)
            elif trivial_group is not None:
                my_result["classification"] = "trivial"
                my_result["rejection"] = "Trivial template group: " + str(trivial_group.pdb)
        
        if my_result["classification"] is None:
            if len(template_groups) == 1:
                my_result["classification"] = "template-based"
            elif len(template_groups) == 2:
                my_result["classification"] = "ab initio"
            else:
                my_result["classification"] = "multi-body"

        result[target_name] = my_result
    return result


if __name__ == "__main__":
    analyzed_searches = json.load(open(sys.argv[1]))
    result = filter_all_targets(analyzed_searches)
    print(json.dumps(result, sort_keys=True, indent=2))
