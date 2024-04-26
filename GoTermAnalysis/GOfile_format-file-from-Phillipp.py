import pandas as pd

df_raw = pd.read_csv("./INPUT2_GO-mapping/GO-Data-Phillipp.log", sep="\t")
df_raw = df_raw.rename({"Entry": "locus_tag"}, axis=1)

# generate output_go_terms files, mapping the locus_tags to the GO terms, 
df_BP = df_raw[["locus_tag", "GOBP"]]                                   # BP: Biological Processes 
df_CC = df_raw[["locus_tag", "GOCC"]]                                   # CC: Cellular Components
df_MF = df_raw[["locus_tag", "GOMF"]]                                   # MF: Molecular Functions
data_all_GOcategories = {"GOBP": df_BP, "GOCC": df_CC, "GOMF": df_MF}
for GO_category in data_all_GOcategories:
    df = data_all_GOcategories[GO_category]                             # define the df we are working with to avoid 
                                                                        # using data_all_GOcategories[GO_category] (too lengthy)
    df = df.rename({f"{GO_category}": "GOTerm"}, axis=1)               # let's see if it's needed by the R script
                                                                        # GO category is the same as the column name we want to rename
    df = df[~df["locus_tag"].isna()]                                    # get rid of all entries that do not have a locus_tag
    df = df[~df["GOTerm"].isna()]                               # gets rid of all entries without a GO annotation
    df["GOTerm"] = df["GOTerm"].str.split(";")          # explode needs a list, GO_category column currently contains a string 
                                                                        # enumerating GO terms separated by ;
    df = df.explode("locus_tag").explode("GOTerm")              # one GO term per row -> duplicating the locus_tags that have multiple GO terms 
    df.to_csv(f"./INPUT2_GO-mapping/output_{GO_category}_terms_PHILLIPP.log", sep="\t", index=False)
    

# generate output_symbols files, mapping the locus_tags to the UniProt entries
df_UniProt = df_raw[["locus_tag", "UniProt"]]
#df_UniProt = df_UniProt.rename({"UniProt": "Entry"}, axis=1) # let's see if it's needed by the R script
df_UniProt.to_csv("./INPUT2_GO-mapping/output_symbols_PHILLIPP.log", sep="\t", index=False)

### !!!!! MISSING: need to make each entry CDS-locus_tag-gene
### and duplicate it as 5'UTR and 3'UTR ???