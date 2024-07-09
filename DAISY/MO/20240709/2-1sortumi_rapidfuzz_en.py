import pandas as pd
from rapidfuzz import process, fuzz
import sys


input_file = sys.argv[1]
output_file = sys.argv[2]


df = pd.read_csv(input_file, sep='\t')

def find_similar_groups(umi_list):
    groups = []
    for umi in umi_list:
        found = False
        for group in groups:
            if fuzz.ratio(umi, group[0]) >= 80:  # 0.8 * 100
                group.append(umi)
                found = True
                break
        if not found:
            groups.append([umi])
    return groups

results = []

for tenx, group in df.groupby('tenx'):
    umi_list = group['umi'].tolist()
    static = group['static'].tolist()
    mutate = group['mutate'].tolist()
    
    similar_groups = find_similar_groups(umi_list)
    
    for similar_group in similar_groups:
        count = len(similar_group)
        representative_umi = similar_group[0]
        indices = [umi_list.index(umi) for umi in similar_group]
        
        results.append({
            'tenx': tenx,
            'umi': representative_umi,
            'static': static[indices[0]],
            'mutate': mutate[indices[0]],
            'overlapumi': count
        })

result_df = pd.DataFrame(results)

result_df.to_csv(output_file, sep='\t', index=False)