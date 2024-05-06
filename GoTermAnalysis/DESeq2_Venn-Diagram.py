# packages needed
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles
import os


# define functions

def Filtering(filename: str, LfcCutoff: float=2.5, padjCutoff: float=0.05) -> pd.DataFrame:
    """
    takes a (DESeq2 results) df as an input and returns the same df filtered after a given LFC and padj value

    filename: name of the DESeq2 file we are working with
    LfcCutoff: log2 fold change cutoff, default is value used by Andreas
    padjCutoff: p adjusted cutoff, default is value used by Andreas

    positive LFC means that the second component is enriched as compared to the first one,
    e. g. Rbp1 in the dataframe "coldGFPinvsout-Rbp1invsout"
    """
    
    df = pd.read_csv(filename, sep="\t")
    if LfcCutoff >= 0:
        df = df.loc[(df["log2FoldChange"] >= LfcCutoff) & (df["padj"] <= padjCutoff)]
    elif LfcCutoff < 0:
        df = df.loc[(df["log2FoldChange"] <= LfcCutoff) & (df["padj"] <= padjCutoff)]

    return df

def EnrichementCount(df: pd.DataFrame) -> int:
    """
    returns the number of enriched transcripts in a DESeq2 file
    """
    count = df.shape[0]

    return count


def dfCompare(df1: pd.DataFrame, df2: pd.DataFrame, on: str="Gene Name") -> tuple:
    """
    returns the number of common entries (default: Gene Names) in two columns of a (DESeq2 results) df
    """
    common_elements = df1.merge(df2, left_on=on, right_on=on, how="inner")[on].unique()

    return (len(common_elements), common_elements)


# code

os.chdir("./summary-table")

## input files
ColdEnriched_Rbp1_file = "../INPUT1_DESeq2_Results/singlelevel/Rbp1-input-cold-vs-std_EVERYTHING.log"
ColdEnriched_GFP_file = "../INPUT1_DESeq2_Results/singlelevel/GFP-input-cold-vs-std_EVERYTHING.log"
Rbp1Bound_file = "../INPUT1_DESeq2_Results/multilevel/coldGFPinvsout-Rbp1invsout_EVERYTHING.log"

input_files = [ColdEnriched_Rbp1_file, ColdEnriched_GFP_file, Rbp1Bound_file]

## file to store the results
summary_file = open(f"summary_VennDiagram_LfcRbp1is-25_LfcGfpis-25.log", "w")

## dfs we are working with
Standard_Rbp1_df = Filtering(ColdEnriched_Rbp1_file, LfcCutoff=2.5)
Standard_GFP_df = Filtering(ColdEnriched_GFP_file, LfcCutoff=2.5)
ColdEnriched_Rbp1_df = Filtering(ColdEnriched_Rbp1_file, LfcCutoff=-2.5)
ColdEnriched_GFP_df = Filtering(ColdEnriched_GFP_file, LfcCutoff=-2.5)
Rbp1Bound_df = Filtering(Rbp1Bound_file, LfcCutoff=2.5)

## enrichement in each file
standard_genes_Rbp1 = EnrichementCount(Standard_Rbp1_df)
summary_file.write(f"Number of genes in standard conditions (Rbp1): {standard_genes_Rbp1}\n")

standard_genes_GFP = EnrichementCount(Standard_GFP_df)
summary_file.write(f"Number of genes in standard conditions (GFP): {standard_genes_GFP}\n")

cold_enriched_genes_Rbp1 = EnrichementCount(ColdEnriched_Rbp1_df)
summary_file.write(f"Number of cold enriched genes (Rbp1): {cold_enriched_genes_Rbp1}\n")

cold_enriched_genes_GFP = EnrichementCount(ColdEnriched_GFP_df)
summary_file.write(f"Number of cold enriched genes (GFP): {cold_enriched_genes_GFP}\n")

Rbp1_bound_genes = EnrichementCount(Rbp1Bound_df)
summary_file.write(f"Number of genes in {Rbp1Bound_file.split("/")[-1]}: {Rbp1_bound_genes}\n")

## genes enriched in Rbp1 and GFP (control) samples under standard conditions
genes_in_Rbp1_and_GFP_STD = dfCompare(Standard_Rbp1_df, Standard_GFP_df)
summary_file.write(f"Number of genes common to {ColdEnriched_Rbp1_file.split("/")[-1]} and {ColdEnriched_GFP_file.split("/")[-1]} under standard conditions: {genes_in_Rbp1_and_GFP_STD[0]}\n")

## genes enriched in Rbp1 and GFP (control) samples under cold conditions
genes_in_Rbp1_and_GFP_COLD = dfCompare(ColdEnriched_Rbp1_df, ColdEnriched_GFP_df)
#df = ColdEnriched_Rbp1_df
#df[df["Gene Name"].isin(genes_in_Rbp1_and_GFP_COLD[1])].to_csv("genes_in_Rbp1_and_GFP.log", sep="\t") # list of genes in Rbp1 and GFP
summary_file.write(f"Number of genes common to {ColdEnriched_Rbp1_file.split("/")[-1]} and {ColdEnriched_GFP_file.split("/")[-1]} in cold conditions: {genes_in_Rbp1_and_GFP_COLD[0]}\n")

## genes enriched in cold conditions and bound to Rbp1
genes_in_Rbp1_and_Rbp1Bound = dfCompare(ColdEnriched_Rbp1_df, Rbp1Bound_df)
summary_file.write(f"Number of genes common to {ColdEnriched_Rbp1_file.split("/")[-1]} and {Rbp1Bound_file.split("/")[-1]}: {genes_in_Rbp1_and_Rbp1Bound[0]}\n")

summary_file.close()

## Venn diagrams
v_ColdBound = venn2(subsets = (cold_enriched_genes_Rbp1-genes_in_Rbp1_and_Rbp1Bound[0], Rbp1_bound_genes-genes_in_Rbp1_and_Rbp1Bound[0], genes_in_Rbp1_and_Rbp1Bound[0]), 
                    set_labels = ('Upregulated in cold', 'Enriched in Rbp1 pulldown'),
                    set_colors=('b','r'),
                    normalize_to=2,
                    alpha=0.5)
c_ColdBound = venn2_circles(subsets = (cold_enriched_genes_Rbp1-genes_in_Rbp1_and_Rbp1Bound[0], Rbp1_bound_genes-genes_in_Rbp1_and_Rbp1Bound[0], genes_in_Rbp1_and_Rbp1Bound[0]),
                            linestyle='dashed', 
                            linewidth=1, 
                            color='black',
                            normalize_to=2,
                            alpha=0.5)
plt.savefig("VennDiagram_ColdBound.jpg")
plt.close()

v_Rbp1GFP_COLD = venn2(subsets = (cold_enriched_genes_Rbp1-genes_in_Rbp1_and_GFP_COLD[0], cold_enriched_genes_GFP-genes_in_Rbp1_and_GFP_COLD[0], genes_in_Rbp1_and_GFP_COLD[0]), 
                  set_labels = ('Upregulated in cold (Rbp1 samples)', 'Upregulated in cold (GFP controls)'),
                  set_colors=('b','g'),
                  normalize_to=2,
                  alpha=0.5)
c_Rbp1GFP_COLD = venn2_circles(subsets = (cold_enriched_genes_Rbp1-genes_in_Rbp1_and_GFP_COLD[0], cold_enriched_genes_GFP-genes_in_Rbp1_and_GFP_COLD[0], genes_in_Rbp1_and_GFP_COLD[0]),
                          linestyle='dashed', 
                          linewidth=1, 
                          color='black',
                          normalize_to=2,
                          alpha=0.5)
plt.savefig("VennDiagram_Rbp1GFP_COLD.jpg")
plt.close()

v_Rbp1GFP_STD = venn2(subsets = (standard_genes_Rbp1-genes_in_Rbp1_and_GFP_STD[0], standard_genes_GFP-genes_in_Rbp1_and_GFP_STD[0], genes_in_Rbp1_and_GFP_STD[0]), 
                  set_labels = ('Standard transcriptome (Rbp1 samples)', 'Standard transcriptome (GFP controls)'),
                  set_colors=('b','g'),
                  normalize_to=2,
                  alpha=0.5)
c_Rbp1GFP_STD = venn2_circles(subsets = (standard_genes_Rbp1-genes_in_Rbp1_and_GFP_STD[0], standard_genes_GFP-genes_in_Rbp1_and_GFP_STD[0], genes_in_Rbp1_and_GFP_STD[0]),
                          linestyle='dashed', 
                          linewidth=1, 
                          color='black',
                          normalize_to=2,
                          alpha=0.5)
plt.savefig("VennDiagram_Rbp1GFP_STD.jpg")
plt.close()