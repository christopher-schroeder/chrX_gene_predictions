import importlib
import yaml

from scripts.predict import predict

with open(snakemake.input.best_params, "r") as stream:
    config = yaml.safe_load(stream)

label = snakemake.wildcards.c
config = config[label]
module = config["module"]
class_name = config["class"]
classifier_class = getattr(importlib.import_module(module), class_name)
params_bet = config.get("params_best", {})
clf = classifier_class(**params_bet)

from sklearn.model_selection import cross_val_score
clf = svm.SVC(kernel='linear', C=1, random_state=42)
scores = cross_val_score(clf, X, y, cv=5)