import sys

from contextlib import contextmanager
import pandas as pd
from sklearn import preprocessing
import numpy as np
from sklearn.calibration import CalibratedClassifierCV
from sklearn.inspection import permutation_importance

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


def get_data_full(filename="raw/nn_data.csv", fill_na=True):
    print("loading data", file=sys.stderr)
    seed = 1094795585
    np.random.seed(seed)
    df = pd.read_csv(filename, sep=",")

    # df_symbol = df["Approved symbol"]
    df = df.drop("Approved symbol", axis=1)
    df = df.drop("chr", axis=1)
    df = df.drop("class", axis=1)

    x = df.values #returns a numpy array
    min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = min_max_scaler.fit_transform(x)
    df = pd.DataFrame(x_scaled)

    if fill_na:
        print("filling NA", file=sys.stderr)
        df.fillna(-1, inplace=True)

    X = df.values.astype(np.float32)
    return X


def get_data(filename="raw/nn_data.csv", fill_na=True):
    print("loading data", file=sys.stderr)
    seed = 1094795585
    np.random.seed(seed)
    df = pd.read_csv(filename, sep=",")

    # df_symbol = df["Approved symbol"]
    df = df.drop("Approved symbol", axis=1)
    df = df.drop("chr", axis=1)

    valid_class = df["class"].isin(["training_validation_brain_disease", "training_validation_dispensable"])
    df = df[valid_class]
    y = df["class"] == "training_validation_brain_disease"

    df = df.drop("class", axis=1)
    features_names = df.columns
    x = df.values #returns a numpy array
    min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = min_max_scaler.fit_transform(x)
    df = pd.DataFrame(x_scaled)

    if fill_na:
        print("filling NA", file=sys.stderr)
        df.fillna(-1, inplace=True)

    X = df.values.astype(np.float32)
    y = y.values.astype(np.float32)
    

    return X, y, features_names

def importance(clf):
    X, y, feature_names = get_data()
    if "predict_proba" in dir(clf) or "decision_function" in dir(clf):
        clf = CalibratedClassifierCV(clf)
    print("fit")
    clf.fit(X, y)
    print("importance")
    return permutation_importance(clf, X, y, n_repeats=30, random_state=1094795585, n_jobs=80, scoring="neg_log_loss"), feature_names


def train(clf, filename="raw/nn_data.csv"):
    X, y, _ = get_data(filename=filename)

    if "predict_proba" in dir(clf) or "decision_function" in dir(clf):
        clf = CalibratedClassifierCV(clf)
    print("fit")
    clf.fit(X, y)
    return clf


# exit()

# # with open(str(snakemake.log), "w") as log:
# #     # with custom_redirection(log):
# #     params = parameter_search(clf, space, fill_na=True)
# #     df = pd.DataFrame.from_dict(params).sort_values("rank_test_score")
# #     df.to_csv(snakemake.output.params, sep="\t", index=False)
