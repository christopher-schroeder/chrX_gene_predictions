import pickle
from predict import predict
from common import CombinedClassifier

clfs = []
for filename in snakemake.input.models:
    print(filename)
    with open(filename, 'rb') as f:
        clfs.append(pickle.load(f))

clf = CombinedClassifier(clfs)
prediction = predict(clf, filename=snakemake.input.raw)

with open(snakemake.output.prediction, "w") as o:
    print(*prediction, sep="\t", file=o)