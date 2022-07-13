import importlib
import yaml
# save the classifier
import pickle

from train_best import train
with open(snakemake.input.best_params, "r") as stream:
    config = yaml.safe_load(stream)

label = snakemake.wildcards.c
config = config[label]
module = config["module"]
class_name = config["class"]
classifier_class = getattr(importlib.import_module(module), class_name)
params_best = config.get("params_best", {})
clf = classifier_class(**params_best)
clf = train(clf, filename=snakemake.input.raw)
with open(snakemake.output.model, 'wb') as o:
    pickle.dump(clf, o)