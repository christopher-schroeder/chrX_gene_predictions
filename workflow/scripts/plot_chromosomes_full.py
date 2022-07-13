import plotly.express as px
import argparse
import pandas as pd
import numpy as np
from itertools import chain

chromosomes = list(chain(map(str, range(1,23)), ["X"]))

df = pd.read_csv("full.tsv", sep="\t")
prediction = df["prediction"].values


is_training_validation_brain_disease = df["class"] == "training_validation_brain_disease"
is_training_validation_dispensable = df["class"] == "training_validation_dispensable"
is_training = is_training_validation_brain_disease | is_training_validation_dispensable
df_training = df[is_training]

class_training_truth = df_training["class"] == "training_validation_brain_disease"
class_training_truth = class_training_truth.values
probability_training = prediction[is_training]
class_training_prediction = probability_training >= 0.5

order = np.flip(np.argsort(probability_training))

class_training_truth = class_training_truth[order]
probability_training = probability_training[order]
class_training_prediction = class_training_prediction[order]

positive = 0
negative = 0
last_p_value = 0

# # try only on chrx

for i, (truth, pred, prob) in enumerate(zip(class_training_truth, class_training_prediction, probability_training)):
    if truth:
        positive += 1
    elif not truth:
        negative += 1
    else:
        assert False, "should not happen!"

    fdr = negative / (i + 1)
    print(fdr, positive, negative, (positive + negative),  pred, prob)
    if fdr > 0.05 and (positive + negative) > 100:
        break
    last_p_value = prob

ratio = []

for chromosome in chromosomes:
    is_chromosome = df["chr"] == chromosome
    is_positive_training = df["class"] == "training_validation_brain_disease"
    is_positive_predicted = prediction > last_p_value

    confirmed_chr = np.sum(is_chromosome & is_positive_training)
    confirmed_predicted_chr = np.sum(is_chromosome & is_positive_training & is_positive_predicted)
    ratio.append(confirmed_predicted_chr / confirmed_chr)

template = "simple_white"

fig = px.bar(x=chromosomes, y=ratio, template=template, labels=chromosomes)
fig.add_hline(y=1.0, line_dash="dot", line_width=1, line_color="black")
fig.update_yaxes(range=[0.0, 1.0])
fig.write_image("test.pdf")