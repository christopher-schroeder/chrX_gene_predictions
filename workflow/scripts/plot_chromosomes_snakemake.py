import pandas as pd
import numpy as np
from plot_chromosomes import plot_sensitivity

df = pd.read_csv(snakemake.input.raw)
prediction = np.genfromtxt(snakemake.input.prediction, delimiter="\t")
chromosomes, precision, sensitivity, text_precision, text_sensitivity = plot_sensitivity(df, prediction, snakemake.output.plot, snakemake.params.classifier)

with open(snakemake.output.table, "w") as o:
    print("chromosomes", "precision", "sensitivity",  "precision text", "sensitivity text", sep="\t", file=o)
    for c, p, s, tp, ts in zip(chromosomes, precision, sensitivity, text_precision, text_sensitivity):
        print(c, p, s, tp, ts, sep="\t", file=o)
