# packages needed
import pandas as pd
import os 
import re
import numpy as np

###############################################################################################################################

# define functions needed

def MSFileFormatting(MSMSdf):
    """
    filter and format the df to make it work as an input for the GO enrichement script

    returns a tuple with the formatted MSMS df with all its columns [0] (to keep a trace of what was done)
    and the one with just the columns needed for GO enrichement [1]
    """

    # filtering according to set cutoff or to already performed cutoff
    LfcCutoff = 2.5
    padjCutoff = 0.05

    MSMSdf = MSMSdf.dropna()
    #MSMSdf = MSMSdf.loc[(MSMSdf["significant"]) == True] # cutoff already performed by MS facility
    MSMSdf = MSMSdf.loc[(MSMSdf["log2FoldChange"] >= LfcCutoff) 
                                          & (MSMSdf["padj"] <= padjCutoff)] # choice of own cutoff

    # move Fasta headers column to the back as it contains interesting information but is not used for computing GO enrichement
    last_column = MSMSdf.pop("Fasta headers")                             
    MSMSdf.insert(MSMSdf.shape[1], "Fasta headers", last_column)

    # only keep the Uniprot ID using RegEx: create new column, populate it with RegEx matches and add it to the dataframe
    uniprot_ids = []
    for prot_id in MSMSdf["Protein IDs"]:
        prot_id = re.findall("[OPQ][0-9][A-Z0-9][A-Z0-9][A-Z0-9][0-9]", prot_id)
        uniprot_ids.append(prot_id)       
    MSMSdf["Uniprot IDs"] = uniprot_ids
    
    # get rid of rows without Uniprot IDs
    test = []
    for i in MSMSdf["Uniprot IDs"]:
        if len(i) == 0:
            test.append(True)
        else:
            test.append(False)
    MSMSdf["Test"] = test
    MSMSdf = MSMSdf.loc[MSMSdf["Test"] == False]
    
    # final formatting
    MSMSdf = MSMSdf.explode("Uniprot IDs") # separate rows with multiple Uniprot IDs in the single Uniprot IDs (one ID per row)
    MSMSdf_GOinput = MSMSdf[["Uniprot IDs", "log2FoldChange", "padj"]]
    MSMSdf_GOinput['Index'] = MSMSdf_GOinput.loc[:, 'Uniprot IDs']
    MSMSdf_GOinput = MSMSdf_GOinput.set_index("Index")

    return (MSMSdf, MSMSdf_GOinput)

###############################################################################################################################

# actual code

## go to right working directory
os.chdir("./INPUT1_DESeq2_Results/msms")

## load data, select relevant columns (one dataset each for standard and cold) and rename the columns to have same names for both datasets to use in a loop
MSMSdf = pd.read_excel("output_Rbp1-coimmunoprecipitation-lcmsms.xlsx", sheet_name=2)
MSMSdf_standard = MSMSdf[["Protein IDs", "Fasta headers", 
                 "log2_ratio_Rbp1_regular", "pValue_Rbp1_regular", "significant_Rbp1_regular"]]
MSMSdf_standard = MSMSdf_standard.rename({"log2_ratio_Rbp1_regular": "log2FoldChange", 
                                          "pValue_Rbp1_regular": "padj", 
                                          "significant_Rbp1_regular": "significant"}, axis=1)
MSMSdf_cold = MSMSdf[["Protein IDs", "Fasta headers", 
                 "log2_ratio_Rbp1_cold", "pValue_Rbp1_cold", "significant_Rbp1_cold"]]
MSMSdf_cold = MSMSdf_cold.rename({"log2_ratio_Rbp1_cold": "log2FoldChange", 
                                          "pValue_Rbp1_cold": "padj", 
                                          "significant_Rbp1_cold": "significant"}, axis=1)

## apply function MSFileFormatting to data
MSMSdf_standard = MSFileFormatting(MSMSdf_standard)[0]
MSMSdf_standard_GOinput = MSFileFormatting(MSMSdf_standard)[1]
MSMSdf_cold = MSFileFormatting(MSMSdf_cold)[0]
MSMSdf_cold_GOinput = MSFileFormatting(MSMSdf_standard)[1]

## save the results
#MSMSdf_standard_GOinput.to_csv("msms_standard.txt", sep="\t")
#MSMSdf_cold_GOinput.to_csv("msms_cold.txt", sep="\t")

# replace UniProt IDs by corresponding gene locus tag
mapping_df = pd.read_csv("../../INPUT2_GO-mapping/output_symbols_PHILLIPP.log", sep="\t")
#mapping_df = mapping_df.dropna()
mapping_df["UniProt"] = mapping_df["UniProt"].str.split(";")
mapping_df = mapping_df.explode("UniProt")
mapping_df = mapping_df.set_index("UniProt")
mapping_df = mapping_df.to_dict()

for i in MSMSdf_cold_GOinput["Uniprot IDs"]:
    try:
        MSMSdf_cold_GOinput = MSMSdf_cold_GOinput.replace(i, mapping_df["locus_tag"][i])
    except:
        MSMSdf_cold_GOinput = MSMSdf_cold_GOinput.replace(i, np.nan)

MSMSdf_cold_GOinput = MSMSdf_cold_GOinput.sort_values("log2FoldChange", ascending=False)

MSMSdf_cold_GOinput.to_csv("test.log", sep="\t")
