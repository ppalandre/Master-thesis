# summary table occurences: number of occurences of each sequence identifier in whole genome, cold RBP1 and standard RBP1
# summary percentages: percentages of occurences with a positive/negative LFC for each sequence identifyier in cold RBP1


# packages needed
import pandas as pd

# load, format and filter relevant files: annotation file of all Synechocystis sp. strain SCC 6803 transcripts and their sequence identifyer, and DESeq2 results
annotation_file = "./summary-table/20210217_syne_onlyUnique_withFeat.gff3"
DESeq2_standard = "./INPUT1_DESeq2_Results/multilevel/stdGFPinvsout-Rbp1invsout_EVERYTHING.log"
DESeq2_cold = "./INPUT1_DESeq2_Results/multilevel/coldGFPinvsout-Rbp1invsout_EVERYTHING.log"

## format annotation file so it has the same sequence identifyer names as the DESeq2 result files
formatted_file = pd.read_csv(annotation_file, sep="\t")["Sequence Identifyer"]
for identifyer in formatted_file:
    if identifyer == "3UTR":
        formatted_file = formatted_file.replace(identifyer, "3'UTR")
    if identifyer == "5UTR":
        formatted_file = formatted_file.replace(identifyer, "5'UTR")
    if identifyer == "misc":
        formatted_file = formatted_file.replace(identifyer, "misc_RNA")
formatted_file.to_csv("./summary-table/20210217_syne_onlyUnique_withFeat_FORMATTED.gff3", sep="\t")

annotation_file = "./summary-table/20210217_syne_onlyUnique_withFeat_FORMATTED.gff3"

## define LFC and padj cutoff values and filter DESeq2 result files using them to get enriched transcripts
LfcCutoff = 0
padjCutoff = 0.01

filtered_standard = pd.read_csv(DESeq2_standard, sep="\t") # standard conditions
filtered_standard = filtered_standard.loc[(filtered_standard["log2FoldChange"] > LfcCutoff) & (filtered_standard['padj'] <= padjCutoff)]
filtered_standard.to_csv("./INPUT1_DESeq2_Results/multilevel/stdGFPinvsout-Rbp1invsout_FILTERED.log", sep="\t")
filtered_standard = "./INPUT1_DESeq2_Results/multilevel/stdGFPinvsout-Rbp1invsout_FILTERED.log"

filtered_cold = pd.read_csv(DESeq2_cold, sep="\t") # cold conditions
filtered_cold = filtered_cold.loc[(filtered_cold["log2FoldChange"] > LfcCutoff) & (filtered_cold['padj'] <= padjCutoff)]
filtered_cold.to_csv("./INPUT1_DESeq2_Results/multilevel/coldGFPinvsout-Rbp1invsout_FILTERED.log", sep="\t")
filtered_cold = "./INPUT1_DESeq2_Results/multilevel/coldGFPinvsout-Rbp1invsout_FILTERED.log"

## final files
all_input_files = [annotation_file, filtered_standard, filtered_cold]


############################################################################################################################

# generate summary dict and populate it with a list of all existing sequence identifyers
summary = {"Sequence Identifyer": []}

for file in all_input_files:
    df = pd.read_csv(file, sep="\t")
    locus_categories = df["Sequence Identifyer"].unique()
    ## use set to unique sequence identifyers, use * to get one set rather than a set of sets
    summary["Sequence Identifyer"] = set([*summary["Sequence Identifyer"], *list(locus_categories)]) 

# populate summary dict with the number of occurence of each Sequence Identifyer in each file
for file in all_input_files:
    df = pd.read_csv(file, sep="\t")
    summary[f"{file}"] = [] # key is the filename
    ## iterate over each type of sequence identifyer: CDS, 5'UTR, 3'UTR, ncRNA, ect.
    for sequence_identifyer in summary["Sequence Identifyer"]:
        if sequence_identifyer in df["Sequence Identifyer"].values:
            summary[f"{file}"].append(df["Sequence Identifyer"].value_counts()[sequence_identifyer])
        else:
            summary[f"{file}"].append(0)

# format and save the final dataframe
## rename columns (current names are the filenames)
summary["Total"] = summary.pop("./summary-table/20210217_syne_onlyUnique_withFeat_FORMATTED.gff3")
summary["Standard"] = summary.pop("./INPUT1_DESeq2_Results/multilevel/stdGFPinvsout-Rbp1invsout_FILTERED.log")
summary["Cold"] = summary.pop("./INPUT1_DESeq2_Results/multilevel/coldGFPinvsout-Rbp1invsout_FILTERED.log")
## transform the set back to a list is needed in order to make a dataframe out of the dict
summary["Sequence Identifyer"] = list(summary["Sequence Identifyer"])
## make and format dataframe
summary = pd.DataFrame(summary)
summary = summary.replace({"rubbish": "other"})
summary = summary.set_index("Sequence Identifyer")
summary = summary.sort_values('Total', ascending=False)
## save dataframe
summary.to_csv("./summary-table/summary-table.log", sep="\t")


#print(summary)

############################################################################################################################

# summary percentages

df = pd.read_csv(DESeq2_cold, sep = "\t")
SeqIds = list(df["Sequence Identifyer"].unique())

percent_summary = {"Sequence Identifyer": SeqIds, 
                   "Total (nb)": [], 
                   "Up (nb)": [], 
                   "Down (nb)": [], 
                   "Up (%)": [], 
                   "Down (%)": []}

LfcThreshold = 1

for id in SeqIds:
    df_filtered = df[df["Sequence Identifyer"] == id]
    total = df_filtered.shape[0]
    percent_summary["Total (nb)"].append(total)
    up = df_filtered[df_filtered["log2FoldChange"] >= LfcThreshold].shape[0]
    percent_summary["Up (nb)"].append(up)
    down = df_filtered[df_filtered["log2FoldChange"] <= -LfcThreshold].shape[0]
    percent_summary["Down (nb)"].append(down)
    percent_up = (up/total)*100
    percent_summary["Up (%)"].append(round(percent_up, 1))
    percent_down = (down/total)*100
    percent_summary["Down (%)"].append(round(percent_down, 1))

percent_summary = pd.DataFrame(percent_summary)
percent_summary = percent_summary.replace({"rubbish": "other"})
percent_summary = percent_summary.set_index("Sequence Identifyer")
percent_summary = percent_summary.sort_values('Total (nb)', ascending=False)
## save dataframe
percent_summary.to_csv("./summary-table/sPercentSummary-table.log", sep="\t")

print(percent_summary)