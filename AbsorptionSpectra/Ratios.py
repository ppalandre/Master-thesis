import pandas as pd 
import seaborn as sns
import matplotlib as plt
import os
import absorption as abs
import samples

samps = ["WT", "delta-rbp1 #6", "RBP1-FLAG #6.6"] #FLAG
samps = ["WT", "delta-rbp1 #4", "RBP1-GFP #4"]     #GFP
NbOfReps = 3

# load data
working_directory = "./20degrees_75uE_FLAG_KW16" #FLAG
dir_name = os.chdir(working_directory) # get current working directory 
rawdata = {}
for file in os.listdir(dir_name): 
    if ".csv" in file:
        rawdata[file] = pd.read_csv(file, sep=";")
formatted_data = {}
for file in rawdata:
    formatted_data[file] = abs.FormatDataframe(rawdata[file]) # no normalization and no statistics needed to get ratios -> makes the code quicker

working_directory = "../20degrees_75uE_KW11"       #GFP
dir_name = os.chdir(working_directory) # get current working directory 
rawdata = {}
for file in os.listdir(dir_name): 
    if ".csv" in file:
        rawdata[file] = pd.read_csv(file, sep=";")
for file in rawdata:
    formatted_data[file] = abs.FormatDataframe(rawdata[file]) # no normalization and no statistics needed to get ratios -> makes the code quicker


# get phycobillin/chlorophyll a ratios for each sample

working_directory = "../ratios"    
dir_name = os.chdir(working_directory)

ratios = {"Timepoint": [],
         "Sample": [],
         "Ratio": []}

for file in formatted_data:
    for s in formatted_data[file]:
        # extract timepoints from filenames
        timepoint = int(file.split("spectrum ")[1].split("h")[0]) 
        ratios["Timepoint"].append(timepoint)
        # get sample name without replicate number and get rid of strain number (e. g. Rbp1-FLAG instead of Rbp1-FLAG #6.6)
        sample = s.split("_")[0] 
        if "#" in sample:
            sample = sample.split(" #")[0]
        ratios["Sample"].append(sample)
        # extract phycobillin/chlorophyll a ratios from spectra
        Abs_phycobillisome = formatted_data[file][s].loc[600.0:650.0].max()
        Abs_chlorophylla = formatted_data[file][s].loc[650.0:700.0].max()
        ratio = round(Abs_phycobillisome / Abs_chlorophylla, 3)
        ratios["Ratio"].append(ratio)

ratios = pd.DataFrame(ratios)
print(ratios)

# plot the ratios

sns.set(style="whitegrid", font_scale=1.2)
sns_plot = sns.barplot(x="Timepoint", y="Ratio", hue="Sample", data=ratios, capsize=.1, errorbar="sd", palette=['grey', '#c43e96', "#317ec2", "#5cb56f"])
sns_plot.set_ylim(0.9, 1.05)
#sns_plot.legend(labels=['WT', 'delta rbp1', 'Rbp1-FLAG', "Rbp1-GFP"])
sns.move_legend(
    sns_plot, "lower center",
    bbox_to_anchor=(.5, 1), ncol=4, title=None, frameon=False,
)
#sns.swarmplot(x="Timepoint", y="Ratio", data=ratios, color="0", alpha=0.35)
fig = sns_plot.get_figure()
fig.savefig("ratios.jpg", transparent=True, dpi=300)