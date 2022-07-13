import importlib
import yaml

import plotly.express as px
import plotly.graph_objects as go

from predict import importance
with open(snakemake.input.best_params, "r") as stream:
    config = yaml.safe_load(stream)

label = snakemake.wildcards.c
config = config[label]
module = config["module"]
class_name = config["class"]
classifier_class = getattr(importlib.import_module(module), class_name)
params_best = config.get("params_best", {})
clf = classifier_class(**params_best)
importance, feature_names = importance(clf)
with open(snakemake.output.importance, "w") as o:
    print(*feature_names, sep="\t", file=o)
    print(*importance["importances_mean"], sep="\t", file=o)

classifier=snakemake.params.classifier

fig = go.Figure()
template = "simple_white"
fig.add_trace(go.Bar(y=importance["importances_mean"], x=feature_names, name="precision", base=0, marker_color='DarkSeaGreen'))
fig.update_layout(
    width=1200,
    height=800,
    title=f"{classifier} Feature Importance"
)
fig.update_xaxes(tickmode='linear')
fig.write_image(snakemake.output.plot)