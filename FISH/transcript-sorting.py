# sort the transcripts of the DESeq2 analysis to find the most abundant and long ones, i. e. the most suitable for FISH

# packages needed
import pandas as pd

# load data
DESeq2_file = "../GoTermAnalysis/INPUT1_DESeq2_Results/multilevel/coldGFPinvsout-Rbp1invsout_EVERYTHING.txt"
mastertable = "../GOTermAnalysis/INPUT2_GO-mapping/GO-Data-Phillipp.txt"

DESeq2Results = pd.read_csv(DESeq2_file, sep="\t")
mapping_df = pd.read_csv(mastertable, sep="\t")


# ...
mapping_df = mapping_df[["Entry", "UniProt", "Protein name", "Protein Sequence lengths", "Fasta Headers", "end", "start", "GOBP", "GOCC", "GOMF"]]
mapping_df = mapping_df.rename({"Entry": "Gene Name"}, axis=1)

merged_df = pd.merge(DESeq2Results, mapping_df, on="Gene Name")
merged_df["rounded LFC"] = round(merged_df["log2FoldChange"], 1)
merged_df["gene length"] = merged_df["end"]-merged_df["start"]
merged_df = merged_df.sort_values(by = ['rounded LFC', 'Protein Sequence lengths'], ascending = [False, False], na_position = 'last')

merged_df.to_csv("FISH-targets.txt", sep="\t")

print(merged_df)
