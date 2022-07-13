import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin

def matthews_correlation(y_true, y_pred):
    y_pred_pos = np.round(np.clip(y_pred, 0, 1))
    y_pred_neg = 1 - y_pred_pos

    y_pos = np.round(np.clip(y_true, 0, 1))
    y_neg = 1 - y_pos

    tp = np.sum(y_pos * y_pred_pos)
    tn = np.sum(y_neg * y_pred_neg)

    fp = np.sum(y_neg * y_pred_pos)
    fn = np.sum(y_pos * y_pred_neg)

    numerator = (tp * tn - fp * fn)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

    return numerator / (denominator + np.finfo(np.float32).eps)

class CombinedClassifier(BaseEstimator, ClassifierMixin):
     def __init__(self, classifiers):
        self.classifiers = classifiers
     def fit(self, X, y):
        pass
     def predict(self, X):
        return (self.predict_proba(X) > 0.5).astype(np.float16)
     def predict_proba(self, X):
        probabilities = [clf.predict_proba(X) for clf in self.classifiers]
        probabilities = np.array(probabilities)
        return np.average(probabilities, axis=0)