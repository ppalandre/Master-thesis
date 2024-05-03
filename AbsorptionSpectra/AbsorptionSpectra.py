import pandas as pd 
import matplotlib as plt
import os
import absorption as abs
import samples


working_directory = "./20degrees_75uE_FLAG_KW16" #FLAG
#working_directory = "./20degrees_75uE_KW11"       #GFP

# load data
dir_name = os.chdir(working_directory) # get current working directory 
rawdata = {}
for file in os.listdir(dir_name): 
    if ".csv" in file:
        rawdata[file] = pd.read_csv(file, sep=";")


# one plot per timepoint containing all samples with technical replicates, FLAG sample

samps = ["WT", "delta-rbp1 #6", "RBP1-FLAG #6.6"] #FLAG
#samps = ["WT", "delta-rbp1 #4", "RBP1-GFP #4"]     #GFP
NbOfReps = 3

formatted_data = abs.ReshapeData(rawdata, samps, norm="750")

for filename in formatted_data:
    fig, ax = abs.Plotting3(formatted_data[filename], filename, samps, samples.samples)
    title = filename.split("spectrum ")[1].split("ColdShock")[0]  # extract info for title from filename (timepoint), the [] are because split returns a list of two strings   
    title = title[:1].upper() + title[1:]                         # make first letter of title uppercase
    ax.set_title(title)                                           # generate title     
    fig.savefig(f"{filename.split("spectrum")[1].split(".")[0]}.jpg", transparent=True, dpi=300)
    #fig.close()
