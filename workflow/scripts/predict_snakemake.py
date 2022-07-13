import pickle
from predict import predict
from common import CombinedClassifier

with open(snakemake.input.model, 'rb') as f:
    clf = pickle.load(f)
prediction = predict(clf, filename=snakemake.input.raw)

with open(snakemake.output.prediction, "w") as o:
    print(*prediction, sep="\t", file=o)