import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# load data
file = "./INPUT1_DESeq2_Results/multilevel/cold GFP invsout - Rbp1 invsout_EVERYTHING.log"
df = pd.read_csv(file, sep="\t")

# prepare df for plotting
LfcCutoff = 2.5
padjCutoff = 0.001
df["log10padj"] = -np.log10(df["padj"]) # compute log10(padj) for the y-axis
df["significant"] = (df["log2FoldChange"] >= LfcCutoff) & (df["padj"] <= padjCutoff)   # differentiate significant and non-significant for color
df_sig = df.loc[df.significant]                                             
df_nonsig = df.loc[~df.significant]
annotate = df[df["Gene Name"] == "psaA_motif"]

#print(annotate)

# plotting
lfc_sig = df_sig["log2FoldChange"]                                                     # differentiate significant and non-significant for color
log10padj_sig = df_sig["log10padj"]
lfc_nonsig = df_nonsig["log2FoldChange"]
log10padj_nonsig = df_nonsig["log10padj"]

n_anno = ["psaA", "rimO"]
x_anno = [-0.54826, 5.837211]
y_anno = [0.297443, 128.877960]

#print(x_anno)

fig, ax = plt.subplots()
ax.scatter(lfc_sig, log10padj_sig, c="#07bcc5", alpha=0.35)                           # alpha is for transparency
ax.scatter(lfc_nonsig, log10padj_nonsig, c="grey", alpha=0.35)                        # alpha is for transparency
ax.axhline(y=-np.log10(padjCutoff), color='black', linestyle='dashed')                # add horizontal and vertical line to show cutoffs
ax.axvline(x=LfcCutoff, color='black', linestyle='dashed')
ax.axvline(x=-LfcCutoff, color='black', linestyle='dashed')
for i, txt in enumerate(n_anno):
    ax.annotate(txt,                                                                  # generate annotation from n_anno, x_anno and y_anno
                xy=(x_anno[i], y_anno[i]), 
                xytext=(x_anno[i]+0.7, y_anno[i]+30), 
                arrowprops=dict(facecolor='black', width=0.0001, headwidth=0.0001))
ax.set_xlabel("log2 fold change")
ax.set_ylabel("-log10(adjusted p-value)")
plt.savefig("./vulcano_plots/VulcanoPlot_DESeq2.jpg", dpi=300)
#plt.show()