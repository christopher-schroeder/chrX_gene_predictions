import sys

import plotly.express as px
import plotly.graph_objects as go

import argparse
import pandas as pd
import numpy as np
from itertools import chain
from plotly.subplots import make_subplots

chromosomes = list(chain(map(str, range(1,23)), ["X"]))

# # try only on chrx

def plot_sensitivity(df, prediction, filename, classifier=""):
    is_training_validation_brain_disease = df["class"] == "training_validation_brain_disease"
    is_training_validation_dispensable = df["class"] == "training_validation_dispensable"
    is_training = is_training_validation_brain_disease | is_training_validation_dispensable
    df_training = df[is_training]

    class_training_truth = df_training["class"] == "training_validation_brain_disease"
    class_training_truth = class_training_truth.values
    probability_training = prediction[is_training]
    class_training_prediction = probability_training >= 0.5

    assert(len(class_training_truth) == len(class_training_prediction))
    assert(len(probability_training) == len(class_training_prediction))

    order = np.flip(np.argsort(probability_training))

    class_training_truth = class_training_truth[order]
    probability_training = probability_training[order]
    class_training_prediction = class_training_prediction[order]

    positive = 0
    negative = 0
    # last_p_value = 0
    prob = 0
    for i, (truth, pred, prob) in enumerate(zip(class_training_truth, class_training_prediction, probability_training)):
        if truth:
            positive += 1
        elif not truth:
            negative += 1
        else:
            assert False, "should not happen!"

        fdr = negative / (i + 1)
        # print(fdr, positive, negative, (positive + negative),  pred, prob)
        if fdr > 0.05 and (positive + negative) > 100:
            break
        # last_p_value = prob

        # print(last_p_value)

    accuracy = []
    sensitivity = []
    precision = []
    text_sensitivity = []
    text_precision = []

    for chromosome in chromosomes:
        is_chromosome = df["chr"] == chromosome
        is_positive_training = df["class"] == "training_validation_brain_disease"
        is_negative_training = df["class"] == "training_validation_dispensable"
        is_positive_predicted = prediction > prob
        # print(prediction)
        # print(np.sum(is_positive_predicted))

        TP = np.sum(is_chromosome & is_positive_training & is_positive_predicted)
        TN = np.sum(is_chromosome & is_negative_training & ~is_positive_predicted)
        FP = np.sum(is_chromosome & is_negative_training & is_positive_predicted)
        FN = np.sum(is_chromosome & is_positive_training & ~is_positive_predicted)

        # confirmed_chr = np.sum(is_chromosome & is_positive_training)
        # confirmed_predicted_chr = np.sum(is_chromosome & is_positive_training & is_positive_predicted)
        sensitivity.append(TP / (TP + FN + sys.float_info.epsilon))
        precision.append(TP / (TP + FP + sys.float_info.epsilon))
        text_sensitivity.append(f"{TP / (TP + FN + sys.float_info.epsilon):.2f} ({TP} / {(TP + FN)})")
        text_precision.append(f"{TP / (TP + FP + sys.float_info.epsilon):.2f} ({TP} / {(TP + FP)})")
    template = "simple_white"

    precision = np.array(precision)

    # df_chromosome = pd.DataFrame()
    # df_chromosome["chromosome"] = chromosomes
    # df_chromosome["sensitivity"] = sensitivity
    # df_chromosome["specificity"] = specificity

    # fig = px.bar(x=chromosomes, y=sensitivity, template=template, labels=chromosomes, text=text)
    fig = make_subplots(rows=1, cols=2, column_widths=[0.5, 0.5], horizontal_spacing = 0.0,
        # specs=[[{"r":-0.015}, {"l":0.05}, {"l":-0.03}]],
        subplot_titles=("precision", "sensitivity"),
        shared_yaxes = True,
        )
    fig.add_trace(go.Bar(y=chromosomes, x=precision, orientation='h', name="precision", text=text_precision, base=0, marker_color='DarkSeaGreen'), row=1, col=1)
    # fig.add_trace(go.Bar(y=chromosomes, x=[0]*len(chromosomes), orientation='h', base=0), row=1, col=2)
    fig.add_trace(go.Bar(y=chromosomes, x=sensitivity, orientation='h', name="sensitivity", text=text_sensitivity, base=0, marker_color='SlateGrey'), row=1, col=2)
    fig.add_vline(x=1.0, line_dash="dot", line_width=1, line_color="black", row=1, col=2)
    fig.add_vline(x=1.0, line_dash="dot", line_width=1, line_color="black", row=1, col=1)
    # fig.update_xaxes(range=[0.0, 0.1], visible=False, row=1, col=2)
    # # fig.update_yaxes(zeroline=False, showline=False, showgrid=False, tickmode = 'array', tickcolor="white", row=1, col=2)
    fig.update_yaxes(showline=False, tickcolor="white", row=1, col=1)
    fig.update_yaxes(visible=False, row=1, col=2)
    fig.update_xaxes(autorange="reversed", row=1, col=1)
    # fig.update_yaxes(visible=False, row=1, col=2)

    fig.update_layout(
        # barmode='stack',
        yaxis_title="chromosome",
        template=template,
        # legend=dict(
        #     x=0,
        #     y=1.0,
        #     bgcolor='rgba(255, 255, 255, 0)',
        #     bordercolor='rgba(255, 255, 255, 0)'
        # ),
        margin=dict(l=20, r=20, t=50, b=20),
        title_text="Stacked Subplots"
    )
    fig.update_layout(showlegend=False, title=f"{classifier}, FDR \u2264 0.05, probability \u2265 {prob:.4f}", title_font_size=12, title_x=0.5)
    fig.update_annotations(font_size=11)
    fig.update_yaxes(
        ticktext=chromosomes,
    )
    # fig.update_xaxes(ticktext=chromosomes)
    fig.write_image(filename)