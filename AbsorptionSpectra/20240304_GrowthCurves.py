# packages needed

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import os

#############################################################################################################################

# functions defined

def GrowthCurveBiolRep(df: pd.DataFrame):
    """
    plots the growth curve of experiments performed with biological replicates
    """

    df["WT"].plot(color="black", label="WT", marker="o")
    df["DualTag"].plot(color="grey", label="GFP", marker="o")
    df["delta-rbp1 #2"].plot(color="#d952e7", linestyle='dashed', label=r"$\Delta \it{rbp1}$ #2", marker="o")
    df["delta-rbp1 #4"].plot(color="#ff00d0", linestyle='dashed', label=r"$\Delta \it{rbp1}$ #4", marker="o")
    df["delta-rbp1 #6"].plot(color="#a11b65", linestyle='dashed', label=r"$\Delta \it{rbp1}$ #6", marker="o")
    df["RBP1-GFP #4"].plot(color="#0c8b30", linestyle='dotted', label="RBP1-GFP #4", marker="o")
    df["RBP1-GFP #11"].plot(color="#124c02", linestyle='dotted', label="RBP1-GFP #11", marker="o")
    df["RBP1-GFP #12"].plot(color="#5cb56f", linestyle='dotted', label="RBP1-GFP #12", marker="o")
    timespan = df.index.tolist()
    plt.xlim(min(timespan),max(timespan)) 
    plt.ylim(0,8.5)                            
    plt.xlabel("Time (h)")                                 
    plt.ylabel(r"$A_{750 nm}$")                                 
    plt.legend(fontsize="8", loc="upper left")                                           

    return plt


def GrowthCurveTechRep(df: pd.DataFrame, samples: list[str], NbOfReps: int, tag: str):
    """
    plots the growth curve of experiments performed with biological replicates

    df: df containing the raw data
    samples: list of the samples (sample 1, sample 2, ..., sampleN)
    NbOfReps: number of replicates for each sample (sample1_1, sample1_2, ...sample1_n, sample2_1, ..., sampleN_n)
    tag: GFP or FLAG, for determining the color of the plots
    """

    # statistics: compute mean, standard deviation and confidence interval
    for sample in samples:
        triplicates = []
        for rep in range(NbOfReps):
            triplicates.append(f"{sample}_{rep+1}")
        df_new = df[triplicates]
        df[f"{sample}_mean"] = df_new.apply(np.mean, axis=1)
        df[f"{sample}_stdev"] = df_new.apply(np.std, axis=1)
        # confidence interval: https://www.pythoncharts.com/python/line-chart-with-confidence-interval/
        df[f"{sample}_ci"] = 1.96 * df[f"{sample}_stdev"] / np.sqrt(NbOfReps) # where does 1.96 come from???
        df[f"{sample}_ci-lower"] = df[f"{sample}_mean"] - df[f"{sample}_ci"]
        df[f"{sample}_ci-upper"] = df[f"{sample}_mean"] + df[f"{sample}_ci"]

    if tag == "GFP":
        # plot the mean with the confidence interval
        fig, ax = plt.subplots()
        x = df['time (h)']
        ## WT
        ax.plot(x, df['WT_mean'], color="black", label="WT")
        ax.fill_between(x, df['WT_ci-lower'], df['WT_ci-upper'], color='black', alpha=.15)
        ## delta-rbp1 #4
        ax.plot(x, df["delta-rbp1 #4_mean"], color="#ff00d0", label=r"$\Delta \it{rbp1}$ #4")
        ax.fill_between(x, df['delta-rbp1 #4_ci-lower'], df['delta-rbp1 #4_ci-upper'], color='#ff00d0', alpha=.15)
        ## delta-rbp1 #6
        ax.plot(x, df["delta-rbp1 #6_mean"], color="#a11b65", label=r"$\Delta \it{rbp1}$ #6")
        ax.fill_between(x, df['delta-rbp1 #6_ci-lower'], df['delta-rbp1 #6_ci-upper'], color='#a11b65', alpha=.15)
        ## RBP1-GFP #4
        ax.plot(x, df["RBP1-GFP #4_mean"], color="#0c8b30", label="RBP1-GFP #4")
        ax.fill_between(x, df['RBP1-GFP #4_ci-lower'], df['RBP1-GFP #4_ci-upper'], color='#0c8b30', alpha=.15)
        ## RBP1-GFP #11
        ax.plot(x, df["RBP1-GFP #11_mean"], color="#124c02", label="RBP1-GFP #11")
        ax.fill_between(x, df['RBP1-GFP #11_ci-lower'], df['RBP1-GFP #11_ci-upper'], color='#124c02', alpha=.15)
        ## format the plot
        timespan = df["time (h)"]
        ax.set_xlim(min(timespan),max(timespan))  
        ax.set_ylim(0,8.5)
        ax.set_xlabel("Time (h)")                                 
        ax.set_ylabel(r"$A_{750 nm}$")
        ax.set_title('Growth Curve')
        ax.legend(fontsize="8", loc="upper left")

    elif tag == "FLAG":
        # plot the mean with the confidence interval
        fig, ax = plt.subplots()
        x = df['time (h)']
        ## WT
        ax.plot(x, df['WT_mean'], color="black", label="WT")
        ax.fill_between(x, df['WT_ci-lower'], df['WT_ci-upper'], color='black', alpha=.15)
        ## delta-rbp1 #6
        ax.plot(x, df["delta-rbp1 #6_mean"], color="#ff00d0", label=r"$\Delta \it{rbp1}$ #6")
        ax.fill_between(x, df['delta-rbp1 #6_ci-lower'], df['delta-rbp1 #6_ci-upper'], color='#ff00d0', alpha=.15)
        ## RBP1-FLAG #6.6
        ax.plot(x, df["RBP1-FLAG #6.6_mean"], color="blue", label="RBP1-FLAG #6.6")
        ax.fill_between(x, df['RBP1-FLAG #6.6_ci-lower'], df['RBP1-FLAG #6.6_ci-upper'], color='blue', alpha=.15)
        ## format the plot
        timespan = df["time (h)"]
        ax.set_xlim(min(timespan),max(timespan)) 
        ax.set_ylim(0,8.5)
        ax.set_xlabel("Time (h)")                                 
        ax.set_ylabel(r"$A_{750 nm}$")
        ax.set_title('Growth Curve')
        ax.legend(fontsize="8", loc="upper left")

    return fig



#############################################################################################################################

# actual code

## set working directory
os.chdir("./growth curves")

## load data for biological replicates
filename20 = "A750_20degrees-75uE-KW8.csv"
filename15 = "A750_15degrees-140uE-KW9.csv"

for filename in [filename20, filename15]:
    df = pd.read_csv(filename, sep="\t")
    df = df.set_index("time (h)")

    ## plot data and save figure
    plt = GrowthCurveBiolRep(df)
    plt.title(f"Growth curve {filename}")                                               
    plt.savefig(f"GrowthCurve_{filename}.jpg", dpi=300)
    plt.close()

## load data for technical replicates
filenameTechRep_GFP = "A750_20degrees-75uE-KW11.csv"
filenameTechRep_FLAG = "A750_20degrees-75uE-FLAG-KW16.csv"
df_GFP = pd.read_csv(filenameTechRep_GFP, sep="\t")
df_FLAG = pd.read_csv(filenameTechRep_FLAG, sep="\t")

## plot the GFP data and save figure
samples_GFP = ["WT", "delta-rbp1 #4", "delta-rbp1 #6", "RBP1-GFP #4", "RBP1-GFP #11"]
NbOfReps_GFP = 3
fig = GrowthCurveTechRep(df_GFP, samples_GFP, NbOfReps_GFP, tag="GFP")
fig.savefig(f"GrowthCurve_{filenameTechRep_GFP}.jpg", dpi=300)

## plot the FLAG data and save figure
samples_FLAG = ["WT", "delta-rbp1 #6", "RBP1-FLAG #6.6"]
NbOfReps_FLAG = 3
fig = GrowthCurveTechRep(df_FLAG, samples_FLAG, NbOfReps_FLAG, tag="FLAG")
fig.savefig(f"GrowthCurve_{filenameTechRep_FLAG}.jpg", dpi=300)