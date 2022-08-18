from fileinput import filename
import repackage
repackage.up()

import numpy as np
from numpy import mean, std

from sklearn.model_selection import cross_val_score, cross_validate, KFold, GridSearchCV, RandomizedSearchCV
from sklearn import preprocessing
from sklearn.metrics import make_scorer
from common import matthews_correlation
import pandas as pd
import sys
from math import prod

score = make_scorer(matthews_correlation, greater_is_better=True)

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
    x = df.values #returns a numpy array
    min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = min_max_scaler.fit_transform(x)
    df = pd.DataFrame(x_scaled)

    if fill_na:
        print("filling NA", file=sys.stderr)
        df.fillna(-1, inplace=True)

    X = df.values.astype(np.float32)
    y = y.values.astype(np.float32)

    return X, y


def cross_validation(model, space={}, label="", fill_na=True):
    if len(space) == 0:
        scores = simple_cross_validation(model, fill_na=fill_na)
    else:
        scores = nested_cross_validation(model, space)
    return scores
    # print(label, "Matthews Correlation: %.3f (%.3f)" % (mean(scores), std(scores)))


def simple_cross_validation(model, fill_na=True):
    seed = 1094795585
    X, y = get_data(fill_na=fill_na)
    print("perform cross validation", file=sys.stderr)
    cv_results = cross_validate(model, X, y, cv=10, verbose=1, n_jobs=10, scoring=score)
    scores = cv_results["test_score"]
    return scores


def parameter_search(model, space, fill_na=True, n_jobs=-1, method="grid"):
    X, y = get_data()

    # configure the cross-validation procedure
    cv_inner = KFold(n_splits=5, shuffle=True, random_state=1)

    # define search

    if prod(len(x) for x in space.values()) > 100:
        method = "random"

    print(space)

    if method == "grid":
        print("Set method to GridSearchCV", file=sys.stderr)
        search = GridSearchCV(model, space, scoring=score, n_jobs=8, cv=cv_inner, refit=True, verbose=1)
    elif method == "random":
        print("Set method to RandomizedSearchCV", file=sys.stderr)
        search = RandomizedSearchCV(model, space, scoring=score, n_jobs=8, cv=cv_inner, refit=True, verbose=1, n_iter=100)
    else:
        assert(False)
    
    search.fit(X,y)
    return search.cv_results_


def nested_cross_validation(model, space, n_jobs=-1, method="grid"):
    X, y = get_data()

    # configure the cross-validation procedure
    cv_inner = KFold(n_splits=5, shuffle=True, random_state=1)

    # define search

    if prod(len(x) for x in space.values()) > 100:
        method = "random"

    print(space)

    if method == "grid":
        print("Set method to GridSearchCV", file=sys.stderr)
        search = GridSearchCV(model, space, scoring=score, n_jobs=9, cv=cv_inner, refit=True, verbose=1)
    elif method == "random":
        print("Set method to RandomizedSearchCV", file=sys.stderr)
        search = RandomizedSearchCV(model, space, scoring=score, n_jobs=9, cv=cv_inner, refit=True, verbose=1, n_iter=100)
    else:
        assert(False)
    # configure the cross-validation procedure
    cv_outer = KFold(n_splits=10, shuffle=True, random_state=1)
    # execute the nested cross-validation
    scores = cross_val_score(search, X, y, scoring=score, cv=cv_outer, n_jobs=8, verbose=1)
    # report performance
    # print('Matthews Correlation: %.3f (%.3f)' % (mean(scores), std(scores)))
    return scores