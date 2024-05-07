import pandas as pd 
import matplotlib as plt
import os
import absorption as abs
import samples

#working_directory = "./Chlorophyll_20degrees_75uE_FLAG_KW16" #FLAG
working_directory = "./Chlorophyll_20degrees_75uE_KW11" #GFP

# load data

## chlorophyll extracts spectra
dir_name = os.chdir(working_directory) # get current working directory 
rawdata = {}
for file in os.listdir(dir_name): 
    if ".csv" in file:
        rawdata[file] = pd.read_csv(file, sep=";")

## OD_750 of the cell cultures
a750_df = working_directory.split("./Chlorophyll_")[1]
a750_df = f"../growth_curves/A750_{a750_df}.csv"
a750_df = pd.read_csv(a750_df, sep="\t")


#samps = ["WT", "delta-rbp1 #6", "RBP1-FLAG #6.6"] #FLAG
samps = ["WT", "delta-rbp1 #4", "delta-rbp1 #6", "RBP1-GFP #4", "RBP1-GFP #11"] #GFP
NbOfReps = 3

data = rawdata.copy()
#file = "20240417_Chlorophyll extraction 0h ColdShock 20degrees_75uE_FLAG.csv"
#data = {file: pd.read_csv(file, sep=";")}

# formatting
formatted_data = {}
for filename in data:    
    formatted_data[filename] = abs.FormatDataframe(rawdata[filename])

# create dict to store extracted values in 
chlorophyll_peaks = {"time (h)": []}
for s in samps:
    for rep in range(NbOfReps):
        chlorophyll_peaks[f"{s}_{rep+1}"] = []

# get absolute quantities of chlorophyll a (in µM/OD) for each timepoint and each replicate of each sample
for filename in formatted_data:
    timepoint = int(filename.split("extraction ")[1].split("h")[0])                     # extract timepoint from the filename as an integer
    a750_working_df = a750_df.loc[a750_df['time (h)'] == timepoint]
    chlorophyll_peaks["time (h)"].append(timepoint)
    for s in samps:
        for rep in range(NbOfReps):
            chlAmax = formatted_data[filename][f"{s}_{rep+1}"].loc[650.0:700.0].max()   # get the maximum of the chlorophyll a peak
            chlA_uM = round(chlAmax * 14, 3)                                            # calculate the absolute quantity of chlorophyll a (in µM) from this value
            a750 = float(a750_working_df[f"{s}_{rep+1}"])                               # normalize on OD_750 of the corresponding culture
            chlA_uM = round(chlA_uM / a750, 3)
            chlorophyll_peaks[f"{s}_{rep+1}"].append(chlA_uM)

# statistics and plot
chlorophyll_peaks = pd.DataFrame.from_dict(chlorophyll_peaks)
chlorophyll_peaks = abs.Statistics(chlorophyll_peaks, samps, NbOfReps)
samps = ["WT", "delta-rbp1 #4", "RBP1-GFP #4"] # we don't want to plot all samples for GFP
fig, ax = abs.GrowthCurveTechRep(chlorophyll_peaks, samps, samples.samples)
ax.set_ylabel("Chloropyll a (µM/OD)")
ax.set_ylim(0,50)
ax.legend(loc="upper right")
fig.savefig("Absolute chlorophyll a content.jpg", dpi=300)

chlorophyll_peaks.to_csv("extracted values.log", sep="\t")   