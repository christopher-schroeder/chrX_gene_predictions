from termios import FIOCLEX
import pandas as pd
import yaml

df = pd.read_csv(snakemake.input.mcc, sep="\t", names=["classifier", "method", "mean", "std"])
df.set_index("classifier", inplace=True)

classifiers = {}

for filename in snakemake.input.best_params:
    data = yaml.safe_load(open(filename, "r"))
    classifiers |= data


with open(snakemake.output.table, "w") as o:
    for c, config in classifiers.items():
        # print(c, df.loc[c])
        d = df.loc[c]
        space = config.get("space", {"space1":{}})
        params = config.get("params", {})
        if type(params) == type(""):
            params = eval(params)
        params = {k:[v] for k, v in params.items()}
        s_string = []
        for s, s_params in space.items():
            s_params = {**params, **s_params}
            s_string.append("\n".join([f"{k}: {', '.join(map(str, v))}" for k,v in s_params.items()]))

        b_params = config["params_best"]
        b_string = []
        b_params = {**params, **b_params}
        b_string.append("\n".join([f"{k}: {v}" for k,v in b_params.items()]))

        s_string="\n\n".join(s_string)
        b_string="\n\n".join(b_string)
        print(c, d["mean"], d["std"], f'"{s_string}"', f'"{b_string}"', sep="\t", file=o)