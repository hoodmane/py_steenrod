import json
with open('output_data/S3-100-external.json') as external_json:  
    external_stems = json.load(external_json)
    with open('output_data/AdemAlgebra3_100.json') as our_json:  
        our_stems = json.load(our_json)
        for i in range(len(our_stems)):
            for j in range(len(our_stems[i])):
                assert external_stems[i][j] == our_stems[i][j]