import pandas as pd
import numpy as np
from plot_chromosomes import plot_sensitivity

df = pd.read_csv(snakemake.input.raw)
prediction = np.genfromtxt(snakemake.input.prediction, delimiter="\t")
plot_sensitivity(df, prediction, snakemake.output.plot, snakemake.params.classifier)
