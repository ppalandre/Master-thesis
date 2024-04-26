# packages needed
import pandas as pd
import os

summary = {"Filename": [], 
           "BP occurences": [], 
           "CC occurences": [], 
           "MF occurences": [], 
           "ALL occurences": []}                                # empty dict to store results

data_directory = "./OUTPUT1_GO-enrichement-RAW"
dir_name = os.chdir(data_directory)                             # get current working directory

for file in os.listdir(dir_name):                               # for each file, count the number of BP, CC and MF occurences and store in summary dict

    df = pd.read_csv(file, sep="\t")

    if df.isin(["BP"]).any().any():
        BP_occurences = df["ONTOLOGY"].value_counts()["BP"]
    else:
        BP_occurences = 0
    if df.isin(["CC"]).any().any():
        CC_occurences = df["ONTOLOGY"].value_counts()["CC"]
    else:
        CC_occurences = 0
    if df.isin(["MF"]).any().any():
        MF_occurences = df["ONTOLOGY"].value_counts()["MF"]
    else:
        MF_occurences = 0
    ALL_occurences = df.shape[0]

    summary["Filename"].append(file)
    summary["BP occurences"].append(BP_occurences)
    summary["CC occurences"].append(CC_occurences)
    summary["MF occurences"].append(MF_occurences)
    summary["ALL occurences"].append(ALL_occurences)

summary = pd.DataFrame(summary)                                 # make dataframe out of summary dict and store as file
summary = summary.set_index("Filename")

summary.to_csv("summary.log", sep="\t")
print(summary)
