import pandas as pd
from rapidfuzz import process, fuzz
import sys


input_file = sys.argv[1]
output_file = sys.argv[2]


df = pd.read_csv(input_file, sep='\t')


def find_similar_groups(static_list, mutate_list):
    groups = []
    for i, static in enumerate(static_list):
        found = False
        for group in groups:
            for j, group_static in enumerate(group['static']):
                if fuzz.ratio(static, group_static) >= 90:  # 0.9 * 100
                    if fuzz.ratio(mutate_list[i], group['mutate'][j]) >= 90:
                        group['static'].append(static)
                        group['mutate'].append(mutate_list[i])
                        found = True
                        break
            if found:
                break
        if not found:
            groups.append({'static': [static], 'mutate': [mutate_list[i]]})
    return groups


results = []


for tenx, group in df.groupby('tenx'):
    static_list = group['static'].tolist()
    mutate_list = group['mutate'].tolist()
    

    similar_groups = find_similar_groups(static_list, mutate_list)
    

    for similar_group in similar_groups:
        count = len(similar_group['static'])
        representative_static = similar_group['static'][0]
        
        overlapmutate = sum(1 for _ in similar_group['mutate'])
        
        results.append({
            'tenx': tenx,
            'static': representative_static,
            'mutate': similar_group['mutate'][0],
            'overlapmutate': overlapmutate
        })

result_df = pd.DataFrame(results)

result_df.to_csv(output_file, sep='\t', index=False)