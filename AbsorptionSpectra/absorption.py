import pandas as pd 
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
import numpy as np

# define functions needed
def FormatDataframe(df):
    """
    formats dataframe df to obtain the absorption values for each wavelength for each sample
    """
    columns_to_drop = ["No.", "Type", "Date/Time", "Note"]
    df_new = df.copy()                                       # create copy of the dataframe we're working with
    df_new = df_new.drop(columns_to_drop, axis=1)            # remove useless columns  
    df_new = df_new.rename(columns={"Name": "Wavelength"})   # rename first column as Wavelength
    df_new = df_new.set_index("Wavelength")                  # change index needed because currently, column containing wavelengths is the index
    df_new = df_new.transpose()                              # switch rows and columns
    df_new = df_new.rename(columns={"DualTag": "GFP"})   
    df_new = df_new[:-1]                                     # remove last row which is not relevant
    df_new = df_new.replace({',': '.'}, regex=True)          # replace all commas by points to be able to convert values into float
    df_new = df_new.astype(float)                            # convert values into float
    
    # change the index value from string to float
    indexes = list(df_new.index)                         # get list from all indexes
    new_indexes = {}                                     # create dict to store old and new indexes as key-value pairs
    for i in range(len(indexes)):                        # for each of the old indexes:
        new = indexes[i].replace(",", ".")               # replace the , by a . in the string
        new = float(new)                                 # then transform the string into a float
        new_indexes[indexes[i]] = new                    # finally add the new old index-new index pair to the dictionary
    df_new.rename(index=new_indexes, inplace=True)       # replace old indexes by new indexes 
    
    return df_new


def Normalize(Abs, Spectrum, Sample, norm="750"):
    """
    normalization of the absorbtion spectra after the phycobillisome 
    absorption maximum and the absorption at 750 nm for the Spectrum
    of a given Sample

    type of normalization can be chosen:
        norm="750" (default): normalization over 750 nm peak (quantity of cells)
        norm="PCA": normalization over 750 nm peak and phycobillin peak
        norm="nope": no normalization
        norm="custom": custom normalization (used for normalization of chlorophyll spectra)
                       which are normalized against the cell density (used to plot the growth curves)
    """
    if norm=="nope":
        return Abs                                                      # no normalization

    else:     
        Abs_750nm = Spectrum[Sample][750.0]                              # get absorbance at 750 nm for the sample
        Abs_phycobillisome = Spectrum[Sample].loc[600.0:650.0].max()     # get maximum of the phycobillisome peak (600-650 nm) 

        if norm=="750":
            return (Abs - Abs_750nm) / (Abs_750nm)                          # normalization over 750 nm peak (quantity of cells)
        elif norm=="PCA":
            return (Abs - Abs_750nm) / (Abs_phycobillisome - Abs_750nm)     # normalization over 750 nm peak and phycobillin peak
        elif norm=="custom":
            Abs_custom = Spectrum[Sample][751.0]
            return (Abs - Abs_custom) / (Abs_custom)
    

def Statistics(df, samples: list[str], NbOfReps: int=3):
    """
    computing the mean value and confidence interval for a set of (technical) replicates

    df: df containing the raw data
    samples: list of the samples (sample 1, sample 2, ..., sampleN)
    NbOfReps: number of replicates for each sample (sample1_1, sample1_2, ...sample1_n, sample2_1, ..., sampleN_n)
    """
    for sample in samples:
            triplicates = []
            for rep in range(NbOfReps):
                #triplicates.append(f"{sample}_{rep+1}_Normalized") # for absorption spectra
                triplicates.append(f"{sample}_{rep+1}")             # for absolute chlorophyll quantities
            df_new = df[triplicates]
            df[f"{sample}_mean"] = df_new.apply(np.mean, axis=1)
            df[f"{sample}_stdev"] = df_new.apply(np.std, axis=1)
            # confidence interval: https://www.pythoncharts.com/python/line-chart-with-confidence-interval/
            df[f"{sample}_ci"] = 1.96 * df[f"{sample}_stdev"] / np.sqrt(NbOfReps) # where does 1.96 come from???
            df[f"{sample}_ci-lower"] = df[f"{sample}_mean"] - df[f"{sample}_ci"]
            df[f"{sample}_ci-upper"] = df[f"{sample}_mean"] + df[f"{sample}_ci"]
    
    return df


def ReshapeData(rawdata:dict, samples: list[str], norm="750", NbOfReps: int=3):
    """
    gets the data in the final shape we need it to work with by formatting and normalizing it
    
    rawdata: dict with file names as keys and files read in as pandas dataframes as values
    type of normalization can be chosen:
        norm="750" (default): ormalization over 750 nm peak (quantity of cells)
        norm="PCA": normalization over 750 nm peak and phycobillin peak
        norm="nope": no normalization
        norm="custom": custom normalization (used for normalization of chlorophyll spectra)
                       which are normalized against the cell density (used to plot the growth curves)
    samples: list of the samples (sample 1, sample 2, ..., sampleN)
    NbOfReps: number of replicates for each sample (sample1_1, sample1_2, ...sample1_n, sample2_1, ..., sampleN_n)
    """

    data = rawdata.copy()

    formatted_data = {}
    for filename in data:
        # formatting
        formatted_data[filename] = FormatDataframe(data[filename])
        # normalization
        for sample in formatted_data[filename]:
            formatted_data[filename]["new_col"] = formatted_data[filename][sample].apply(lambda row: Normalize(row, formatted_data[filename], sample, norm)) # applies the function Normalize to each row of the dataframe and stores it in new column
            formatted_data[filename] = formatted_data[filename].rename(columns={"new_col": f"{sample}_Normalized"}) 
        # statistics
        formatted_data[filename] = Statistics(formatted_data[filename], samples, NbOfReps)
    return formatted_data


def Ratio(Spectrum, Sample):
    """
    returns the ratio of the maximum of the phycobillisome peak (600-650 nm)
    and the maximum of the chlorophyll a peak (651-750 nm) for the Spectrum
    of a given Sample
    """
    Abs_phycobillisome = Spectrum[Sample].loc[600.0:650.0].max()  # get maximum of the phycobillisome peak (600-650 nm)
    Abs_chlorophyllA = Spectrum[Sample].loc[651.0:750.0].max()    # get maximum of the chlorophyll a peak (651-750 nm)
    
    # calculation of ratio PCA/ChlA according to J. Myers, J.-R. Graham, R. T. Wang Plant Physiology (1980)
    # the correction factors take into account the fact that the chlorophyll a and phycobillisome peaks are not fully separated
    ratio = (1.0162*Abs_phycobillisome-0.2612*Abs_chlorophyllA)/(1.0162*Abs_chlorophyllA-0.063*Abs_phycobillisome)

    return ratio 


def Plotting(df, filename, color=1):

    """
    one plot per timepoint containing all samples
    
    df: dataframe to be plotted
    filename: name of the file the dataframe comes from
    color: default value of 1 colors the mutants in the same color and the complementants in the same color
            works only with normalized data
    """
    for s in df:
        df[f"{sample}_Normalized"].plot(color=SampleFormatting[s]["color"], 
                                        label=SampleFormatting[s]["label"], 
                                        linestyle=SampleFormatting[s]["linestyle"])
    
    # plot formatting                      
    plt.xlim(400,750)                                             # define limits of x axis to 400-750 nm
    plt.xlabel("Wavelength (nm)")                                 # label x-Axis as "Wavelength (nm)"
    plt.ylabel("Absorbance (AU)")                                 # label y-Axis as "Absorbance (AU)"
    #plt.legend(fontsize="8", loc="upper center")                  # generate legend
    title = filename.split("spectrum ")[1].split("ColdShock")[0]  # extract info for title from filename (timepoint), the [] are because split returns a list of two strings   
    title = title[:1].upper() + title[1:]                         # make first letter of title uppercase
    plt.title(title)                                              # generate title     
    
    return plt


def Plotting2(Data, Sample, Normalized=1):
    """
    one plot per sample containing all timepoints

    As a default, the data is normalized using the function Normalize, 
    set Normalized=0 to work with non-normalized data
    """
    for filename in Data:
        if Normalized: 
            Data[filename] = Data[filename].filter(regex='_Normalized', axis=1)   # only filters if the data was normalized
            legend = filename.split("spectrum ")[1].split("ColdShock")[0]
            legend = legend[:1].upper() + legend[1:]
            Data[filename][f"{Sample}"].plot(label=f"{legend}")
        else:
            legend = filename.split("spectrum ")[1].split("ColdShock")[0]
            legend = legend[:1].upper() + legend[1:]
            Data[filename][Sample].plot(label=f"{legend}")
    plt.xlim(400,750)                                       # define limits of x axis to 400-750 nm
    plt.xlabel("Wavelength (nm)")                           # label x-Axis as "Wavelength (nm)"
    plt.ylabel("Absorbance (AU)")                           # label y-Axis as "Absorbance (AU)"
    plt.legend(fontsize="8", loc="upper center")             # generate legend
    plt.title(Sample)                                       # generate title
    
    return plt


def Plotting3(df: pd.DataFrame, SampleList: list[str], SampleFormatting: dict):
    """
    plotting one plot per timepoint containing all samples with technical replicates

    df: df containing the raw data
    filename: name of the file the df comes from
    SampleList: list of all samples we want to plot
    SampleFormatting: dict mapping every sample to a color, linestyle and lable for plot formatting
    """

    # plot the mean with the confidence interval
    df['Wavelength (nm)'] = df.index
    fig, ax = plt.subplots()
    x = df['Wavelength (nm)']
    for s in SampleList: 
        ax.plot(x, df[f'{s}_mean'], color=SampleFormatting[s]["color"], label=SampleFormatting[s]["label"])
        ax.fill_between(x, df[f'{s}_ci-lower'], df[f'{s}_ci-upper'], color=SampleFormatting[s]["color"], alpha=.15)
    ## format the plot
    ax.set_xlim(400,750) 
    #ax.set_ylim(ymin=0, ymax=1.1) # to give all plots the same y axis
    ax.set_xlabel("Wavelength (nm)")                                 
    ax.set_ylabel(r"$A_{750 nm}$")
    #ax.set_xticklabels(x, fontsize=12)
    #ax.set_yticklabels(fontsize=12)
    ax.legend(fontsize="12", loc="upper right")

    return fig, ax



def GrowthCurveBiolRep(df, SampleFormatting: dict):
    """
    plots the growth curve of experiments performed with biological replicates
    """

    for s in df:
        df[s].plot(color=SampleFormatting[s]["color"], 
                    label=SampleFormatting[s]["label"], 
                    linestyle=SampleFormatting[s]["linestyle"],
                    marker='o')
        
    timespan = df.index.tolist()
    plt.xlim(min(timespan),max(timespan)) 
    plt.ylim(0,8.5)                            
    plt.xlabel("Time (h)")                                 
    plt.ylabel(r"$A_{750 nm}$")                                 
    plt.legend(fontsize="8", loc="upper left")                                           

    return plt


def GrowthCurveTechRep(df, SampleList, SampleFormatting: dict):
    """
    plots the growth curve of experiments performed with biological replicates

    """
    fig, ax = plt.subplots()
    x = df['time (h)']
    for s in SampleList: 
        ax.plot(x, df[f'{s}_mean'], color=SampleFormatting[s]["color"], label=SampleFormatting[s]["label"], marker= 'o')
        ax.fill_between(x, df[f'{s}_ci-lower'], df[f'{s}_ci-upper'], color=SampleFormatting[s]["color"], alpha=.15)

    ## format the plot
    timespan = df["time (h)"]
    ax.set_xlim(min(timespan),max(timespan))  
    ax.set_ylim(0,8.5)
    ax.set_xlabel("Time (h)")                                 
    ax.set_ylabel(r"$A_{750 nm}$")
    ax.set_title('Growth Curve')
    ax.legend(loc="upper left")

    return fig, ax