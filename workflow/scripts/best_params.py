import pandas as pd
import yaml
import sys
import ast

df = pd.concat(pd.read_csv(filename, sep="\t") for filename in snakemake.input.params)
df = df.sort_values("mean_test_score", ascending=False)

data = df.iloc[0].to_dict()
data["params_best"] = ast.literal_eval(data["params"]) | snakemake.params.classifier.get("params", {})
data = data | snakemake.params.classifier
data = {snakemake.params.name : data}
# yaml.dump(data, sys.stderr, allow_unicode=True)
with open(snakemake.output.yaml, "w") as o:
    yaml.dump(data, o, allow_unicode=True)