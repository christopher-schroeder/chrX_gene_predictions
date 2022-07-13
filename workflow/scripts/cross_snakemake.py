import sys
import importlib
from nested_cross import cross_validation
from contextlib import contextmanager

@contextmanager
def custom_redirection(fileobj):
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    sys.stdout = fileobj
    sys.stderr = fileobj
    try:
        yield fileobj
    finally:
        sys.stdout = old_stdout
        sys.stderr = old_stderr


config = snakemake.params.config
spacename = snakemake.wildcards.space
label = snakemake.wildcards.c
space = config.get("space", {"space1":{}})[spacename]
module = config["module"]
class_name = config["class"]
classifier_class = getattr(importlib.import_module(module), class_name)
init_parameter = config.get("params", {})

clf = classifier_class(**init_parameter)

# clf = LogisticRegression(random_state=seed, max_iter=10000, solver="saga", dual=False)

with open(str(snakemake.log), "w") as log:
    with custom_redirection(log):
        mean_score, std_score = cross_validation(clf, label=label, fill_na=True, space=space)

with open(snakemake.output.mcc, "w") as o:
    print(label, "Matthews Correlation: %.2f (%.2f)" % (mean_score, std_score), file=o)