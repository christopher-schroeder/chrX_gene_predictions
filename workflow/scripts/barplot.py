from cProfile import label
from heapq import merge
from pdb import line_prefix
from turtle import color
import plotly.express as px
import argparse
import pandas as pd
import plotly.graph_objects as go

print(snakemake.input.mcc)

data = pd.read_csv(snakemake.input.mcc, sep="\t", header=None)
data.columns = ["classifier", "score_name", "mean", "std"] + list(map(lambda x: f"v{x}", range(len(data.columns) - 4)))

v_names = [f"v{i}" for i in range(len(data.columns) - 4)]

plus_minus = u"\u00B1"
def sign(x):
    if x > 0:
        return " "
    else:
        return ""

labels = [f'{sign(row["mean"])}{row["mean"]:.3f} {plus_minus} {row["std"]:.3f}    {row["classifier"]}' for _, row in data.iterrows()]
print(labels)

template = "simple_white"
# fig = px.bar(data, x='classifier', y='mean', template=template, error_y="std", labels=dict(classifier="Classifier", mean="Mean Matthews Correlation Coefficient (MCC)"), height=700)
# fig = px.box(merged_data, x="classifier", y='v', points="all", boxmean=True, template=template, height=700, labels=dict(classifier="Classifier", mean="Mean Matthews Correlation Coefficient (MCC)"))
# fig.add_hline(y=0.7, line_dash="dot", line_width=1, line_color="black")


fig = go.Figure()

for i, row in data.iterrows():
    # fig.add_trace(go.Box(
    #     y=row[v_names],
    #     name=f'mean: {sign(row["mean"])}{row["mean"]:.3f}    {row["classifier"]}',
    #     boxmean=True,
    #     jitter=0,
    #     pointpos=0,
    #     boxpoints='all', # represent all points
    #     # line=dict(width=0.75),
    #     marker = dict(size=3, color = "#d62728"),
    #     line = dict(color = 'rgba(0,0,0,0)'),
    #     fillcolor = 'rgba(0,0,0,0)'
    # ))
    fig.add_trace(go.Scatter(
        mode='markers',
        x=[labels[i]] * len(v_names),
        # name=labels[i],
        y=row[v_names],
        marker = dict(
            color = "#e377c2",
        )
    ))

    fig.add_trace(go.Scatter(
        mode='markers',
        x=[labels[i]],
        # name=f'mean: {sign(row["mean"])}{row["mean"]:.3f}    {row["classifier"]}',
        y=[row['mean']],
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[row['std']],
            visible=True,
            thickness=1,
            color = 'black',),
        marker = dict(
            size = 10,
            line = dict(
                width = 1
            ),
            symbol="line-ew"
        ),
    ))




# for idx in range(len(fig.data)):
#     fig.data[idx].x = labels

# fig.data[0].x = labels

fig.add_hline(y=0.0, line_dash="dot", line_width=1, line_color="black")

fig.update_layout(
    template=template,
    showlegend=False,
    height=800,
    width=800)
fig.update_yaxes(range=[-0.2, 1.0])
fig.update_xaxes(tickangle=90)
# # fig.data = fig.data[::-1] 
fig.write_image(snakemake.output.pdf)