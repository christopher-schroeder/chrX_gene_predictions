# import sys
# 
# from nested_cross import cross_validation
# from contextlib import contextmanager

# @contextmanager
# def custom_redirection(fileobj):
#     old_stdout = sys.stdout
#     old_stderr = sys.stderr
#     sys.stdout = fileobj
#     sys.stderr = fileobj
#     try:
#         yield fileobj
#     finally:
#         sys.stdout = old_stdout
#         sys.stderr = old_stderr

import yaml
import importlib
import sys
from contextlib import contextmanager
from sklearn.ensemble import StackingClassifier
import pandas as pd
from sklearn import preprocessing
from sklearn.calibration import CalibratedClassifierCV
import numpy as np

from predict import predict

clfs = []

for filename in snakemake.input.best_params:
    with open(filename, "r") as stream:
        best_params = yaml.safe_load(stream)
        classifier = list(best_params.keys())[0]
        config = best_params[classifier]
        module = config["module"]
        class_name = config["class"]
        classifier_class = getattr(importlib.import_module(module), class_name)
        params_best = config.get("params_best", {})
        clf = classifier_class(**params_best)
        clfs.append((classifier, clf))

clf = StackingClassifier(clfs)
prediction = predict(clf)

with open(str(snakemake.output.prediction), "w") as o:
    print(*prediction, file=o, sep="\t")