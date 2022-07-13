import sys

from contextlib import contextmanager
from tkinter.tix import X_REGION
import pandas as pd
from sklearn import preprocessing
import numpy as np
from sklearn.calibration import CalibratedClassifierCV
from sklearn.inspection import permutation_importance
import shap as shp

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

def get_data_full(filename="raw/data.csv", fill_na=True, gene=None):
    print("loading data", file=sys.stderr)
    seed = 1094795585
    np.random.seed(seed)
    df = pd.read_csv(filename, sep=",")

    symbols = df["Approved symbol"]

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

    df = df[symbols == gene]

    X = df.values.astype(np.float32)
    return X, df.columns


def get_data(filename="raw/data.csv", fill_na=True):
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
    print("importance")
    X, y, feature_names = get_data()
    return permutation_importance(clf, X, y, n_repeats=30, random_state=1094795585, n_jobs=80, scoring="neg_log_loss"), feature_names


def shap(clf, gene, filename="raw/data.csv"):
    X, feature_names = get_data_full(filename=filename, gene=gene)

    def _predict(X):
        if "predict_proba" in dir(clf):
            prediction = clf.predict_proba(X)
        else:
            prediction = clf.predict(X)
        if len(prediction.shape) == 2:
            prediction = prediction[:,1]
        return prediction

    prediction = _predict(X)
    print(prediction)
    print("shap")
    X_train, y_train, feature_names = get_data()
    print("shap2")
    explainer = shp.KernelExplainer(_predict, shp.sample(X_train, 100))
    print("shap3")
    shap_values = explainer.shap_values(X, nsamples=500)
    X = np.array([[f"{x:0.3}" for x in X[0]]])
    shp.force_plot(explainer.expected_value, shap_values, X, matplotlib=True, show=False, figsize=(20,3), feature_names=feature_names)
    return shap_values[0]


def predict(clf, filename="raw/data.csv"):
    X, _ = get_data_full(filename=filename)

    def _predict(X):
        if "predict_proba" in dir(clf):
            prediction = clf.predict_proba(X)
        else:
            prediction = clf.predict(X)
        if len(prediction.shape) == 2:
            prediction = prediction[:,1]
        return prediction

    prediction = _predict(X)

    # prediction = np.round(prediction).astype(int)
    return prediction


# exit()

# # with open(str(snakemake.log), "w") as log:
# #     # with custom_redirection(log):
# #     params = parameter_search(clf, space, fill_na=True)
# #     df = pd.DataFrame.from_dict(params).sort_values("rank_test_score")
# #     df.to_csv(snakemake.output.params, sep="\t", index=False)
