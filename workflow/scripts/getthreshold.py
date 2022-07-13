import plotly.express as px
import argparse
import pandas as pd
import numpy as np
from itertools import chain

df = pd.read_csv(snakemake.input.raw)
valid_class = df["class"].isin(["training_validation_brain_disease", "training_validation_dispensable"])
df = df[valid_class]
true_class = (df["class"] == "training_validation_brain_disease").values

prediction = np.genfromtxt(snakemake.input.prediction, delimiter="\t")
prediction = prediction[valid_class]
prediction = prediction >= 0.5

positive = 0
negative = 0
last_p_value = 0

for i, (index, row) in enumerate(df.iterrows()):
    if row["truth"] == 1:
        positive += 1
    elif row["truth"] == 0:
        negative += 1
    else:
        assert False, "should not happen!"

    fdr = negative / (i + 1)
    print(fdr, positive, negative, (positive + negative),  row["prediction"])
    # if fdr > 0.01 and (positive + negative) > 100:
    #     break
    last_p_value = row["prediction"]

print("threshold", last_p_value)
print("fdr", fdr)
