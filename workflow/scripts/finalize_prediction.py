import pandas as pd
import numpy as np
import os


df = pd.read_csv(snakemake.input.raw, sep=",")
for filename in snakemake.input.predictions:
    classifier = os.path.basename(filename).rsplit(".", 1)[0]
    data = np.genfromtxt(filename, delimiter="\t")
    df[classifier]=data

top_classifier = snakemake.params.top
top_df = df[top_classifier]

df["Top_mean"] = top_df.mean(axis=1)
top_x = top_df >= 0.5

booleanDictionary = {True: 'X', False: ''}
top_x = top_x.replace(booleanDictionary)
top_x = top_x.add_suffix("_class")

df = pd.concat([df, top_x], axis=1)

df.to_csv(snakemake.output.final, sep="\t", index=False)

np.savetxt(snakemake.output.prediction, df["Top_mean"].values, delimiter="\t")