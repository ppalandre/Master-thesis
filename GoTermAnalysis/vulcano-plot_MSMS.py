import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 26}) # define font size for whole plot

# load data
file = "./INPUT1_DESeq2_Results/msms/output_Rbp1-coimmunoprecipitation-lcmsms.xlsx"
df = pd.read_excel(file, sheet_name=2)
df = df[["Protein IDs", "log2_ratio_Rbp1_cold", "pValue_Rbp1_cold", "log10_pValue_Rbp1_cold", "significant_Rbp1_cold"]]  # log2_ratio_Rbp1_cold is LFQ
df = df.dropna()

# prepare df for plotting
LfcCutoff = 2.5
padjCutoff = 0.05
#df["log10padj"] = -np.log10(df["padj"]) # compute log10(padj) for the y-axis
df_sig = df.loc[(df["log2_ratio_Rbp1_cold"] >= LfcCutoff) & (df["pValue_Rbp1_cold"] <= padjCutoff)]                                            
df_nonsig = df.loc[(df["log2_ratio_Rbp1_cold"] < LfcCutoff) | (df["pValue_Rbp1_cold"] > padjCutoff)] 
annotate = df[df["Protein IDs"] == "FLAG_sfGFP"] #rnhA: Q55801;A0A068N4X2, ssb1: Q55499;A0A068MX81, Rbp1_FLAG;Q57014;A0A068MYQ9, FLAG_sfGFP

print(annotate)

# plotting
lfc_sig = df_sig["log2_ratio_Rbp1_cold"]                                                     # differentiate significant and non-significant for color
log10padj_sig = df_sig["log10_pValue_Rbp1_cold"]
lfc_nonsig = df_nonsig["log2_ratio_Rbp1_cold"]
log10padj_nonsig = df_nonsig["log10_pValue_Rbp1_cold"]

n_anno = ["rnhA", "ssb", "Rbp1-FLAG", "GFP-FLAG"]
x_anno = [3.075292, 3.536258, 6.566367, -3.06033]
y_anno = [6.670715, 8.305219, 10.547889, 3.998427]

fig, ax = plt.subplots(figsize=(12,10))
ax.scatter(lfc_sig, log10padj_sig, c="orange", alpha=0.5)                             # alpha is for transparency
ax.scatter(lfc_nonsig, log10padj_nonsig, c="grey", alpha=0.35)                        # alpha is for transparency
ax.axhline(y=-np.log10(padjCutoff), color='black', linestyle='dashed')                # add horizontal and vertical line to show cutoffs
ax.axvline(x=LfcCutoff, color='black', linestyle='dashed')
ax.axvline(x=-LfcCutoff, color='black', linestyle='dashed')
for i, txt in enumerate(n_anno):
    ax.annotate(txt,                                                                  # generate annotation from n_anno, x_anno and y_anno
                xy=(x_anno[i], y_anno[i]), 
                xytext=(x_anno[i]-0.4, y_anno[i]+0.2), 
                arrowprops=dict(arrowstyle="-", facecolor='white', alpha=0))          # alpha=0 makes the arrow totally transparent
ax.set_xlabel("log2 fold change")
ax.set_ylabel("-log10(adjusted p-value)")
plt.savefig("./vulcano_plots/VulcanoPlot_MSMS.jpg", dpi=300)
#plt.show()