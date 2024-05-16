import pandas as pd 
import numpy as np
import matplotlib as plt
import os
import absorption as abs
import samples

##############################################################################################################################################################

#working_directory = "./20degrees_75uE_FLAG_KW16" #FLAG
#working_directory = "./20degrees_75uE_KW11"       #GFP

# load data
#dir_name = os.chdir(working_directory) # get current working directory 
#rawdata = {}
#for file in os.listdir(dir_name): 
 #   if ".csv" in file:
  #      rawdata[file] = pd.read_csv(file, sep=";")

# one plot per timepoint containing all samples with technical replicates, FLAG sample

#samps = ["WT", "delta-rbp1 #6", "RBP1-FLAG #6.6"] #FLAG
#samps = ["WT", "delta-rbp1 #4", "RBP1-GFP #4"]     #GFP
#NbOfReps = 3

#formatted_data = abs.ReshapeData(rawdata, samps, norm="750")

#for filename in formatted_data:
 #   fig, ax = abs.Plotting3(formatted_data[filename], samps, samples.samples)
 #   title = filename.split("spectrum ")[1].split("ColdShock")[0]  # extract info for title from filename (timepoint), the [] are because split returns a list of two strings   
 #   title = title[:1].upper() + title[1:]                         # make first letter of title uppercase
 #   ax.set_title(title)                                           # generate title     
 #   fig.savefig(f"../spectra/{filename.split("spectrum")[1].split(".")[0]}.jpg", transparent=True, dpi=300)
    #fig.close()

##############################################################################################################################################################

# plot both FLAG and GFP data in the same plot, 0h

## load data
fileFLAG = "./20degrees_75uE_FLAG_KW16/20240417_Absorption spectrum 0h ColdShock 20degrees 75uE_FLAG.csv"
fileGFP = "./20degrees_75uE_KW11/20240313_Absorption spectrum 0h ColdShock 20degrees 75muE.csv"
rawdataFLAG =  pd.read_csv(fileFLAG, sep=";")
rawdataFLAG = abs.FormatDataframe(rawdataFLAG)
rawdataGFP =  pd.read_csv(fileGFP, sep=";")
rawdataGFP = abs.FormatDataframe(rawdataGFP)
formatted_data = pd.concat([rawdataFLAG, rawdataGFP], axis=1)
## normalization (here we directly replace the value in the column by the normalized value without generating a new column to store the normalized values)
def Norm(A, A750):                          # define function which for each absorption A returns the absorption normalized after the absorption at 750 nm, A750
    return (A - A750) / (A750)
A750 = list(formatted_data.loc[750.0])      # get list of all A750
for ind in formatted_data.index:            # for each row of the dataframe, get the absorption A, then apply the function Norm, then replace the row by the result of the normalization
    A = list(formatted_data.loc[ind])
    A = list(map(Norm, A, A750))
    formatted_data.loc[ind] = A
## renaming columns so all replicates have the same name: drop the _1/2/3 and the #{strain number} in sample names
for sample in formatted_data:
    formatted_data = formatted_data.rename(columns={sample: sample.split("_")[0].split(" #")[0]})

## statistics
samps = {"WT": 6, "delta-rbp1": 9, "RBP1-GFP": 6, "RBP1-FLAG": 3}         
for sample in samps:
    temp_df = formatted_data.filter(regex=sample, axis=1)
    formatted_data[f"{sample}_mean"] = temp_df.apply(np.mean, axis=1)
    formatted_data[f"{sample}_stdev"] = temp_df.apply(np.std, axis=1)
    # confidence interval: https://www.pythoncharts.com/python/line-chart-with-confidence-interval/
    formatted_data[f"{sample}_ci"] = 1.96 * formatted_data[f"{sample}_stdev"] / np.sqrt(samps[sample]) # get number of replicates from samps dictionary
    formatted_data[f"{sample}_ci-lower"] = formatted_data[f"{sample}_mean"] - formatted_data[f"{sample}_ci"]
    formatted_data[f"{sample}_ci-upper"] = formatted_data[f"{sample}_mean"] + formatted_data[f"{sample}_ci"]

## plotting
fig, ax = abs.Plotting3(formatted_data, samps, samples.samples)
#title = filename.split("spectrum ")[1].split("ColdShock")[0]  # extract info for title from filename (timepoint), the [] are because split returns a list of two strings   
#title = title[:1].upper() + title[1:]                         # make first letter of title uppercase
#ax.set_title(title)                                           # generate title     
fig.savefig("./spectra/spectrum_0h_ALL.jpg", transparent=True, dpi=300)
#fig.close()

print(formatted_data)