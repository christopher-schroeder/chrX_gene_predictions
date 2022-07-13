import plotly.express as px
import argparse
import pandas as pd

data = pd.read_csv(snakemake.input.mcc, sep="\t", names=["classifier", "score_name", "mean", "std"])

template = "simple_white"
fig = px.bar(data, x='classifier', y='mean', template=template, error_y="std", labels=dict(classifier="Classifier", mean="Mean Matthews Correlation Coefficient (MCC)"))
fig.add_hline(y=0.7, line_dash="dot", line_width=1, line_color="black")

fig.update_yaxes(range=[-0.2, 1.0])
fig.write_image(snakemake.output.pdf)