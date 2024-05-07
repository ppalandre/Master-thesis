import pandas as pd 
import os
import absorption as abs
import samples


working_directory = "./growth_curves"

# load data
dir_name = os.chdir(working_directory) # get current working directory 
rawdata = {}
for file in os.listdir(dir_name): 
    if ".csv" in file:
        rawdata[file] = pd.read_csv(file, sep="\t")

# biological replicates
ToPlot = ["A750_20degrees-75uE-KW8.csv", "A750_15degrees-140uE-KW9.csv"]

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