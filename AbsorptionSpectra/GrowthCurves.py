import pandas as pd 
import numpy as np
import os
import absorption as abs
import samples

#################################################################################################################################################################

working_directory = "./growth_curves"

# load data
dir_name = os.chdir(working_directory) # get current working directory 
rawdata = {}
for file in os.listdir(dir_name): 
    if ".csv" in file:
        rawdata[file] = pd.read_csv(file, sep="\t")

#################################################################################################################################################################

# biological replicates
ToPlot = ["A750_20degrees-75uE-KW8.csv", "A750_15degrees-140uE-KW9.csv"]

#################################################################################################################################################################

# technical replicates GFP
SampleList = ["WT", "delta-rbp1 #4", "delta-rbp1 #6", "RBP1-GFP #4", "RBP1-GFP #11"]
file = "A750_20degrees_75uE_KW11.csv"                                            
formatted_data = {file: rawdata[file]}
# no normalization but statistics function works with "_Normalized" columns
for sample in formatted_data[file]:
    formatted_data[file]["new_col"] = formatted_data[file][sample].apply(lambda row: abs.Normalize(row, formatted_data[file], sample, norm="nope")) # applies the function Normalize to each row of the dataframe and stores it in new column
    formatted_data[file] = formatted_data[file].rename(columns={"new_col": f"{sample}_Normalized"}) 
# statistics
formatted_data[file] = abs.Statistics(formatted_data[file], SampleList, NbOfReps=3)
# plotting
SampleList = ["WT", "delta-rbp1 #4", "RBP1-GFP #4"] # we don't want to plot all samples
fig, ax = abs.GrowthCurveTechRep(formatted_data[file], SampleList, samples.samples)
fig.savefig(f"{file.split(".")[0]}.jpg", transparent=True, dpi=300)


# technical replicates FLAG
SampleList = ["WT", "delta-rbp1 #6", "RBP1-FLAG #6.6"]
file = "A750_20degrees_75uE_FLAG_KW16.csv"
formatted_data = {file: rawdata[file]}
# no normalization but statistics function works with "_Normalized" columns
for sample in formatted_data[file]:
    formatted_data[file]["new_col"] = formatted_data[file][sample].apply(lambda row: abs.Normalize(row, formatted_data[file], sample, norm="nope")) # applies the function Normalize to each row of the dataframe and stores it in new column
    formatted_data[file] = formatted_data[file].rename(columns={"new_col": f"{sample}_Normalized"}) 
# statistics
formatted_data[file] = abs.Statistics(formatted_data[file], SampleList, NbOfReps=3)
# plotting
fig, ax = abs.GrowthCurveTechRep(formatted_data[file], SampleList, samples.samples)
fig.savefig(f"{file.split(".")[0]}.jpg", transparent=True, dpi=300)

#################################################################################################################################################################

# everything in 1 plot

## functions for statistics
def ci(stddev, n):
    return 1.96 * stddev / np.sqrt(n)

def ciUpper(mean, ci):
    return mean + ci

def ciLower(mean, ci):
    return mean - ci

## load data
fileFLAG = "A750_20degrees_75uE_FLAG_KW16.csv"
fileGFP = "A750_20degrees_75uE_KW11.csv"
rawdata[fileFLAG] = rawdata[fileFLAG].set_index("time (h)")
rawdata[fileGFP] = rawdata[fileGFP].set_index("time (h)")
formatted_data = pd.concat([rawdata[fileFLAG], rawdata[fileGFP]], axis=1)
formatted_data = formatted_data.drop(["new_col"], axis=1)
for sample in formatted_data:
    formatted_data = formatted_data.rename(columns={sample: sample.split("_")[0].split(" #")[0]}) # drom the _1/2/3 and the #{strain number} in sample names

## statistics: WT: 6 replicates, delta-rbp1: 9 replicates, Rbp1-GFP: 6 replicates, Rbp1-FLAG: 3 replicates
samps = {"WT": [6,6,6,6], "delta-rbp1": [9,9,9,9], "RBP1-GFP": [6,6,6,6], "RBP1-FLAG": [3,3,3,3]} # list of 4 times the number of reps to use in map for ci calculation
stat = {"time (h)": []}
for timepoint in formatted_data.index.to_list():
    stat["time (h)"].append(timepoint)
    for s in samps:
        df = formatted_data[s]
        stat[f"{s}_mean"] = list(df.mean(axis=1)) # compute mean value of each row
        stat[f"{s}_stdev"] = list(df.std(axis=1)) # compute mean value of each row
        stat[f"{s}_ci"] = list(map(ci, stat[f"{s}_stdev"], samps[s]))
        stat[f"{s}_ci-lower"] = list(map(ciUpper, stat[f"{s}_mean"], stat[f"{s}_ci"]))
        stat[f"{s}_ci-upper"] = list(map(ciLower, stat[f"{s}_mean"], stat[f"{s}_ci"]))
stat = pd.DataFrame(stat)

## plotting
SampleList = ["WT", "delta-rbp1", "RBP1-FLAG", "RBP1-GFP"]
fig, ax = abs.GrowthCurveTechRep(stat, SampleList, samples.samples)
fig.savefig(f"growth-curves_ALL.jpg", transparent=True, dpi=300)

print(stat)