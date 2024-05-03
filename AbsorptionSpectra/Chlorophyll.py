import pandas as pd 
import matplotlib as plt
import os
import absorption as abs
import samples

working_directory = "./Chlorophyll_20degrees_75uE_FLAG_KW16" #FLAG

# load data
dir_name = os.chdir(working_directory) # get current working directory 
rawdata = {}
for file in os.listdir(dir_name): 
    if ".csv" in file:
        rawdata[file] = pd.read_csv(file, sep=";")


samps = ["WT", "delta-rbp1 #6", "RBP1-FLAG #6.6"] #FLAG
NbOfReps = 3

#data = rawdata.copy()
file = "20240417_Chlorophyll extraction 0h ColdShock 20degrees_75uE_FLAG.csv"
data = {file: pd.read_csv(file, sep=";")}

formatted_data = {}
for filename in data:
    # formatting
    formatted_data[filename] = abs.FormatDataframe(data[filename])
    
    # normalization
    ## get normalization information from growth curves data (we normalize to the cellular density)
    a750 = working_directory.split("./Chlorophyll_")[1]
    a750 = f"../growth_curves/A750_{a750}.csv"
    a750 = pd.read_csv(a750, sep="\t")
    timepoint = int(file.split("Chlorophyll extraction ")[1].split("h")[0]) # to filter out a750 to get the data corresponding to the timepoint of the data we are plotting (extracted from the filename)
    a750 = a750.loc[a750['time (h)'] == timepoint]
    a750 = a750.replace(timepoint, 751.0)
    a750 = a750.rename({"time (h)": "Wavelength (nm)"}, axis=1) # to have the same structure as our data to concatenate a750 with our data
    a750 = a750.set_index("Wavelength (nm)")
    formatted_data[filename] = pd.concat([formatted_data[filename], a750])
    for sample in formatted_data[filename]:
        formatted_data[filename]["new_col"] = formatted_data[filename][sample].apply(lambda row: abs.Normalize(row, formatted_data[filename], sample, norm="custom")) # applies the function Normalize to each row of the dataframe and stores it in new column
        formatted_data[filename] = formatted_data[filename].rename(columns={"new_col": f"{sample}_Normalized"}) 
    formatted_data[filename] = formatted_data[filename].drop([751.0])
    # statistics
    formatted_data[filename] = abs.Statistics(formatted_data[filename], samps, NbOfReps)

for filename in formatted_data:
    fig, ax = abs.Plotting3(formatted_data[filename], samps, samples.samples)
    title = filename.split("extraction ")[1].split("ColdShock")[0]  # extract info for title from filename (timepoint), the [] are because split returns a list of two strings   
    title = title[:1].upper() + title[1:]                         # make first letter of title uppercase
    ax.set_title(title)                                           # generate title     
    fig.savefig(f"{filename.split("extraction ")[1].split(".")[0]}.jpg", transparent=True, dpi=300)

print(formatted_data[file])