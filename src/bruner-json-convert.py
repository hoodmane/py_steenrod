import json

def readFile(input_file):
    f = open(input_file, 'r')
    contents = f.read()
    f.close()
    return contents
    
def writeFile(output_file, contents):
    f = open(output_file, 'w')
    f.write(contents)
    f.close()    

def brunerFileToJSONFile(input_file, output_file):
    bruner_str = readFile(input_file)    
    module = parseBrunerString(bruner_str)
    json_str = json.dumps(module)
    writeFile(output_file, json_str)

def JSONFileToBrunerFile(input_file, output_file):
    json_str = readFile(input_file)
    module = json.loads(json_str)
    bruner_str = moduleToBrunerString(module)
    writeFile(output_file, bruner_str)

def parseBrunerString(module):
    lines = [s for s in module.split("\n") if s]
    num_gens = int(lines[0])
    gens =  [int(s) for s in lines[1].split(" ") if s]
    gens_dict = {}
    idx = 0
    for g in gens:
        gens_dict[str(idx)] = g
        idx += 1
    actions = []
    for l in lines[2:]:
        l = [int(s) for s in l.split(" ") if s]
        act = {"Sq" : l[1], "input" : str(l[0]), "output" : []}
        for e in l[3:]:
            act["output"].append([1, str(e)])
        actions.append(act)
    return {"generators" : gens_dict, "actions" : actions }
    
def moduleToBrunerString(obj):
    gens = obj["generators"]
    actions = obj["actions"]
    num_gens = len(gens)
    gen_degrees = gens.values()
    # We need a map (index in generator list) --> (generator name) because Bruner modules only refer to elts by their index.
    l = [(gens[k], k) for k in gens]
    l.sort()
    gen_name_to_idx = { str(ind) : x[1] for ind, x in enumerate(l)}
    output_lines = [str(num_gens)]
    output_lines += [" ".join(map(str,gen_degrees))]
    for act in actions:
        action_list = [gen_name_to_idx[act["input"]]]
        action_list += [str(act["Sq"])]
        action_list += [str(len(act["output"]))] 
        action_list += [str(gen_name_to_idx[x[1]]) for x in act["output"]]
        output_lines += [ " ".join(action_list) ]
    return "\n".join(output_lines)
    
        
        
