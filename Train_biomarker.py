import os
import numpy as np
import pandas as pd
import joblib
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve


df = pd.read_excel('Data/DEPs.xlsx','protein biomarker combination')
indp = 'CC1	CC2	CC3	CC4	CC5	CC6	CC7	CC8	CC9	CC10	CC11	CC12	CC13	CC14	CC15	CC16	CC17	CC18'.split('\t')
indn = 'HC1	HC2	HC3	HC4	HC5	HC6	HC7	HC8	HC9	HC10	HC11	HC12'.split('\t')
ind = indp + indn
lis=[0,1,2,3,4]
X = df.loc[lis, ind].T.values
y = np.array([1 if x in indp else 0 for x in ind])


def get_rmse(records_real, records_predict):
    if len(records_real) == len(records_predict):
        return math.sqrt(sum([(x - y) ** 2 for x, y in zip(records_real, records_predict)]) / len(records_real))
    else:
        return None


def Train_fold(i,r, X_train, y_train, X_test, y_test):
    Clist = [1]
    clf = LogisticRegressionCV(Cs=Clist, penalty='l2', fit_intercept=False, cv=5, solver='lbfgs', n_jobs=4,
                               refit=True, class_weight='balanced', multi_class='ovr')
    clf.fit(X_train, y_train)
    result = []
    coef = clf.coef_.ravel()
    result.append(clf.intercept_.tolist() + coef.tolist())
    y_test_scores = clf.predict_proba(X_test)[:, 1]
    AUC = roc_auc_score(y_test, y_test_scores)
    joblib.dump(clf, filename=str(n) + f'/Protein_CC vs HC.model')
    return y_test, y_test_scores, result, AUC


tests = np.array([])
y_tests = np.array([])
y_test_scores = np.array([[]])
n='Model'
for r in range(20):
    Skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=r)
    y_tests, y_test_scores = [], []
    for i, (train, test) in enumerate(Skf.split(X, y)):
        X_train = X[train]
        y_train = y[train]
        X_test = X[test]
        y_test = y[test]
        y_test, y_test_score, result, AUC = Train_fold(i,r, X_train, y_train, X_test, y_test)
        y_tests = y_tests + y_test
        y_test_scores = y_test_scores + X_test
        print(AUC)
    rmse1 = get_rmse(y_tests, y_test_scores)
















