import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t
import os
import samples

working_directory = "./growth_curves"
os.chdir(working_directory)
file = "A750_20degrees-75uE-FLAG-KW16.csv"
df = pd.read_csv(file, sep="\t")

samps = ["WT", "delta-rbp1 #6", "RBP1-FLAG #6.6"]
NbOfReps = 3

fig, ax = plt.subplots()
x = df['time (h)']

for s in samps:
    triplicates = []
    for rep in range(NbOfReps):
        triplicates.append(f"{s}_{rep+1}")
    df_new = df[triplicates]
    df[f"{s}_mean"] = df_new.apply(np.mean, axis=1)
    df[f"{s}_stdev"] = df_new.apply(np.std, axis=1)
    # critical value: https://blog.finxter.com/5-best-ways-to-plot-95-confidence-interval-error-bars-with-python-pandas-dataframes-and-matplotlib/
    df[f"{s}_t_critical_up"] = t.ppf(q=0.975, df=NbOfReps-1, loc=df[f"{s}_mean"], scale=df[f"{s}_stdev"])
    df[f"{s}_t_critical_down"] = t.ppf(q=0.025, df=NbOfReps-1, loc=df[f"{s}_mean"], scale=df[f"{s}_stdev"])
    # margin of error
    df[f"{s}_moe_up"] = df[f"{s}_t_critical_up"] * (df[f"{s}_stdev"]/np.sqrt(NbOfReps))
    df[f"{s}_moe_down"] = df[f"{s}_t_critical_down"] * (df[f"{s}_stdev"]/np.sqrt(NbOfReps))
    ax.plot(x, df[f'{s}_mean'], color=samples.samples[s]["color"], label=samples.samples[s]["label"], marker='o')
    ax.errorbar(x=df["time (h)"], y=df[f"{s}_mean"], yerr=df[f"{s}_moe_up"], ecolor=samples.samples[s]["color"])

plt.xlabel("Time (h)")
plt.ylabel(r"$A_{750 nm}$")
plt.legend()
plt.show()

#y=df.filter(regex=("_mean"))

#print(df)
#print(y)

# Sample DataFrame
#df = pd.DataFrame({'data': np.random.normal(0, 1, size=100)})
# Calculate mean, standard deviation and the size of the dataset
#mean = df['data'].mean()
#std = df['data'].std()
#n = len(df)
# Calculate the critical value
#t_critical = t.ppf(q=0.975, df=n-1)
# Calculate the margin of error
#moe = t_critical * (std/np.sqrt(n))
# Plot error bar
#plt.errorbar(x='Measurement', y=mean, yerr=moe, fmt='o')
#plt.title('95% Confidence Interval Error Bar')
#plt.show()