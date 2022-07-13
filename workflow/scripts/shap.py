import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt

import pickle
from predict import shap
from common import CombinedClassifier

clfs = []
for filename in snakemake.input.models:
    print(filename)
    with open(filename, 'rb') as f:
        clfs.append(pickle.load(f))

clf = CombinedClassifier(clfs)
shap_values = shap(clf, snakemake.wildcards.gene, filename=snakemake.input.raw)
plt.subplots_adjust(left=0.05, right=0.95, top=0.5, bottom=-0.4)
ax = plt.gca()
plt.text(0, 1.3, snakemake.wildcards.gene, fontsize="xx-large", c="dimgray", transform=ax.transAxes)
plt.savefig(snakemake.output.plot)
with open(snakemake.output.values, "w") as o:
    print(*shap_values, sep="\t", file=o)