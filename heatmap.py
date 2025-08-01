import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

csvfile = "mutation_matrix.csv"
df = pd.read_csv(csvfile,  index_col=0)


selected = ["BRP_24_G_13", "BRP_6_H_11", "BRP_2_H_17", "BRP_6_G_7",
            "BRP_1_H_20", "BRP_1_G_3", "BRP_2_H_7", "BRP_6_G_18",
            "BRP_6_H_18", "BRP_1_H_31", "BRP_1_H_15", "BRP_1_G_23",
            "BRP_2_G_18", "BRP_6_H_8", "BRP_24_H_20", "BRP_6_G_14",
            "BRP_24_H_16", "BRP_24_H_7", "BRP_2_G_14", "BRP_24_G_17",
            "BRP_24_G_6", "BRP_2_H_13", "BRP_2_G_25", "BRP_1_G_19"]

df_sel = df[selected]
# Sort columns alphanumerically
sorted_columns = sorted(df_sel.columns, key=str)
df_sel = df_sel[sorted_columns]



plt.figure(figsize=(15, 12))
sns.heatmap(df_sel,
            cmap="viridis",
            norm=plt.Normalize(vmin=df_sel.min().min(), vmax=df_sel.max().max()),
            cbar_kws={'label': 'Mutation Count'},
            yticklabels=True)

plt.title("Mutation Context Heatmap — Selected Samples")
plt.xlabel("Sample")
plt.ylabel("Mutation Context")
plt.tight_layout()
#plt.show()
plt.savefig("heatmap.png", dpi=300)
