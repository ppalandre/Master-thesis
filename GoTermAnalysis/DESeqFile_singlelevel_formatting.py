import pandas as pd

# singlelevel or multilevel
level = "singlelevel"

# summary of length of files before and after dropping na values
summary_file = open(f"./INPUT1_DESeq2_Results/{level}/summary_before-vs-after-na-drop.txt", "w")

print("Working on file: ")

for sheet_number in range(12):
    # read each sheet of the data (range of the number of sheets in the Excel file)
    df_raw = pd.read_excel(f"./INPUT1_DESeq2_Results/{level}/DESeq2_new_summary.xlsx", sheet_name=sheet_number) 
    file_title = list(df_raw)[0].split("Name ")[1]                               # file title is the sheet name
    print(file_title)                                                            # makes sure each sheet gets read once

    # filter the dataframe
    columns_to_keep = [list(df_raw)[0], list(df_raw)[2], list(df_raw)[6]] 
    df_working = df_raw[columns_to_keep]                                         # keep only the columns needed for GO terms analysis: Gene Name, Log2FoldChange and padj    
    #df_working = df_working.dropna()                                            # get rid of missing values (not necessary, R script also checks for missing values)                                                            
    """
    Note on p-values set to NA: some values in the results table can be set to NA for one of the following reasons:

    If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p value and adjusted p value will all be set to NA.
    If a row contains a sample with an extreme count outlier then the p value and adjusted p value will be set to NA. These outlier counts are detected by Cook’s distance. Customization of this outlier filtering and description of functionality for replacement of outlier counts and refitting is described below
    If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p value will be set to NA. Description and customization of independent filtering is described below

    from: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
    """

    # format the dataframe
    df_working = df_working.rename({list(df_working)[1]: "log2FoldChange", 
                                    list(df_working)[2]: "padj"}, axis=1)        # gets rid of the "" in the column names
    df_working = df_working.rename({list(df_working)[0]: "Gene Name"}, axis=1)   # makes column name readable by GO Enrichment R script
    #for gene_name in df_working["Gene Name"]:
    #    df_working = df_working.replace(gene_name, gene_name.split("\"")[1])    # gets rid of the "" in the gene names          
    for gene_name in df_working["Gene Name"]:
        if "-" not in gene_name:
            df_working = df_working.replace(gene_name, "other-"+gene_name)     # adds prefix to all entries that do not have one
    # splits after the first -, puts what is left (suffix) of it in "Gene Name" and what is right (prefix) in "Sequence Identifyer"
    df_working[["Sequence Identifyer", "Gene Name"]] = df_working["Gene Name"].str.split("-", n=1, expand=True)
    for gene_name in df_working["Gene Name"]:
        df_working = df_working.replace(gene_name, gene_name.split("-")[0])      # gets rid of further suffixes to isolate gene name
    df_working = df_working.set_index("Gene Name")
    df_working = df_working.dropna() 

    # last filtering to get rid of the "Sequence Identifyer" column and optionally further filter the dataframe
    df_with_everything = df_working[["log2FoldChange", "padj"]]
    df_with_everything = df_working
    df_CDS = df_working[df_working["Sequence Identifyer"] == "CDS"]
    df_CDS = df_CDS[["log2FoldChange", "padj"]]
    df_CDS_and_UTRs = df_working[(df_working["Sequence Identifyer"] == "CDS") | (df_working["Sequence Identifyer"] == "5'UTR") | (df_working["Sequence Identifyer"] == "3'UTR")]
    df_CDS_and_UTRs = df_CDS_and_UTRs[["log2FoldChange", "padj"]]
    
    # save formatted and filtered dataframes to tsv files for use in GO Enrichment R script
    df_with_everything.to_csv(f"./INPUT1_DESeq2_Results/{level}/{file_title}_EVERYTHING.txt", sep="\t")
    df_CDS.to_csv(f"./INPUT1_DESeq2_Results/{level}/{file_title}_CDS.txt", sep="\t")
    df_CDS_and_UTRs.to_csv(f"./INPUT1_DESeq2_Results/{level}/{file_title}_CDS_and_UTRs.txt", sep="\t")


    summary_file.write(f"{file_title}\n")
    summary_file.write(f"before dropping na: {df_raw.shape[0]}\n")
    summary_file.write(f"after dropping na: {df_working.shape[0]}\n")
    summary_file.write("#####\n")

summary_file.close()