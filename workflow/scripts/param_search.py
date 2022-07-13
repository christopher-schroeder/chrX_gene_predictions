import sys
import importlib
from nested_cross import parameter_search
from contextlib import contextmanager
import pandas as pd


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


with open(str(snakemake.log), "w") as log:
    # with custom_redirection(log):
    params = parameter_search(clf, space, fill_na=True)
    df = pd.DataFrame.from_dict(params).sort_values("rank_test_score")
    df.to_csv(snakemake.output.params, sep="\t", index=False)

# with open(snakemake.output.mcc, "w") as o:
#     print(clf)
#     print(label, "Matthews Correlation: %.2f (%.2f)" % (mean_score, std_score), file=o)