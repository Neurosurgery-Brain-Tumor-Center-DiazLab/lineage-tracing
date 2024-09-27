from ete3 import Tree, TreeStyle
import os
os.chdir("/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/results/")

# Load the tree
tree = Tree("/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate/SB28_DAISY_DAY0_sub2_A3_S122_L007/BCall.nexus.treefile")


# Create a TreeStyle object for circular layout
ts = TreeStyle()
ts.mode = "c"
ts.show_leaf_name = False
ts.show_branch_length = False
ts.show_branch_support = True
# ts.scale = 0

# Render the tree
tree.render("all_circular_tree.png", w=2000, tree_style=ts, dpi = 300)


import pandas as pd
import matplotlib.pyplot as plt
import squarify

import pandas as pd

# Sample dataframe
data = {
    'Static Barcode': ['GAATTCAAAT', 'TTGTGCGTAT', 'TGACAGTCTT', 'TTGGAGTTCC', 'TCCTGCAGCA', 'TCGGTTTTAT', 
                       'GGAGAGTCTC', 'CGCATTTTAT', 'TCACCTACTT', 'ATTTCAGTTC', 'AAATTACTCA'],
    'Mutable Barcode': [
        'TTTGCTACCTATTACTAGGACAAGTGGAGGCGTTGAGAATGTCTCGGAATGTGCCGGAAATTTGCATCGTATTACTAGGACAATTCGAGTCCTTGAGGAAGTCTCTCGATGTGCCGGAAA',
        'TTTGCTACCTATTACTAGGACAAGTGGAGGCGTTGAGAATGTTTCGGAATGTGCCGGAAATTTGCATCGTATTACTAGGACAATTCGAGTCCTTGAGGAAGTCTCTCGATGTGCCGGAAA',
        'TTTGCTACCTATTACTAGGACAAGTGGAGGCGTTGAGAATGTCTCGGAATGTGCCGGAAATTTGCATCGTATTACTAGGACAATTCGAGTCCTTGAGGAAGTCTCTCGATGTGCCGGAAA',
        'TTTGCTACCTATTACTAGGACAAGTGGAGGCGTTGAGAATGTCTCGGAATGTGCCGGAAATTTGCATCGTATTACTAGGACAATTCGAGTCCTTGAGGAAGTCTCTCGATGTGCCGGAAA',
        'TTTGCTACCTATTACTAGGACAAGTGGAGGCGTTGAGAATGTCTCGGAATGTGCCGGAAATTTGCATCGTATTACTAGGACAATTCGAGTCCTTGAGGAAGTCTCTCGATGTGCCGGAAA',
        'TTTGCTACCTATTACTAGGACAAGTGGAGGCGTTGAGAATGTCTCGGAATGTGCCGGAAATTTGCATCGTATTACTAGGACAATTCGAGTCCTTGAGGAAGTCTCTCGATGTGCCGGAAA',
        'TTTGCTACCTATTACTAGGACAAGTGGAGGCGTTGAGAATGTCTCGGAATGTGCCGGAAATTTGCATCGTATTATTAGGACAATTCGAGTCCTTGAGGAAGTCTCTCGATGTGCCGGAAA',
        'TTTGCTACCTATTACTAGGACAAGTGGAGGCGTTGAGAATGTCTCGGAATGTGCCGGAAATTTGCATCGTATTACTAGGACAATTCGAGTCCTTGAGGAAGTCTCTCGATGTGCCGGAAA',
        'TTTGCTACCTATTGGCGTATTGTGCCGGAAATTTGCATCGTATTACTAGGCCTTCTCTCGATGTGCCGGAAA',
        'TTTGCTACCTATTACTAGGACAAGTGGAGGCGTTGGAATGTGCCGGAAATTTGCATCGTATTACTAGGACAATTCGAGTCCTTGAGGAAGTCTCTCGATGTGCCGGAAA',
        'TTTGCTACCTATTACTAGGACAAGTGGAGGCGTTGAGAATGTCTCGGAATGTGCCGGAAATTTGCATCGTATTACTAGGACAATTCGAGTCCTTGAGGAAGTCTCTCGATGTGCCGGAAA'
    ]
}

df = pd.DataFrame(data)

# Calculate sizes (e.g., lengths of mutable barcodes)
df['Size'] = df['Mutable Barcode'].apply(len)

import matplotlib.pyplot as plt
import squarify

# Normalize sizes
sizes = df['Size'].values
labels = df['Static Barcode'].values

squarify.plot(sizes=sizes, label=labels, alpha=0.8)
plt.title('Circle Packing of Static Barcodes')
plt.axis('off')

# Display plot
plt.show()


import numpy as np

def circle_packing(sizes, labels, ax):
    sizes = np.sqrt(sizes)  # Convert areas to radii
    positions = []
    for radius in sizes:
        while True:
            x, y = np.random.uniform(radius, 1-radius, size=2)
            if all(np.sqrt((x - px)**2 + (y - py)**2) > radius + pr for px, py, pr in positions):
                positions.append((x, y, radius))
                break

    for (x, y, radius), label in zip(positions, labels):
        circle = plt.Circle((x, y), radius, alpha=0.6)
        ax.add_patch(circle)
        plt.text(x, y, label, ha='center', va='center', fontsize=8)

# Plot circle packing
fig, ax = plt.subplots(figsize=(8, 8))
circle_packing(sizes, labels, ax)
ax.set_aspect('equal')
ax.axis('off')
plt.show()


