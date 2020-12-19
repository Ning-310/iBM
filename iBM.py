import os
import numpy as np
import pandas as pd
import joblib
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
import random
import math


df = pd.read_excel('Data/DEPs.xlsx','DEPs_CC-sp')
df1 = pd.read_excel('Data/DEPs.xlsx','DEPs_CC vs HC')
df2 = pd.read_excel('Data/DEPs.xlsx','DEPs_CC vs AC')
indp = 'CC1	CC2	CC3	CC4	CC5	CC6	CC7	CC8	CC9	CC10	CC11	CC12	CC13	CC14	CC15	CC16	CC17	CC18'.split('\t')
indn = 'HC1	HC2	HC3	HC4	HC5	HC6	HC7	HC8	HC9	HC10	HC11	HC12'.split('\t')
ind = indp + indn
indp1 = 'CC1	CC2	CC3	CC4	CC5	CC6	CC7	CC8	CC9	CC10	CC11	CC12	CC13	CC14	CC15	CC16	CC17	CC18'.split('\t')
indn1 = 'AC1	AC2	AC3	AC4	AC5	AC6	AC7	AC8	AC9	AC10	AC11	AC12	AC13	AC14	AC15	AC16	AC17	AC18	AC19	AC20	AC21	AC22	AC23	AC24	AC25	AC26	AC27	AC28	AC29	AC30	AC31	AC32	AC33	AC34	AC35	AC36	AC37	AC38	AC39	AC40	AC41	AC42	AC43'.split('\t')
ind1 = indp1 + indn1


#CCG
def random_int_list(start, stop, length):
    start, stop = (int(start), int(stop)) if start <= stop else (int(stop), int(start))
    length = int(abs(length)) if length else 0
    random_list = []
    for i in range(length):
        random_list.append(random.randint(start, stop))
    return random_list
rd1,rd2=[],[]
while rd1.__len__()<10:
    rd01, rd02= [], []
    rd0=random_int_list(0, 43, 5)
    lis = df['Accession'].T.values.tolist()
    lis1=df1['Accession'].T.values.tolist()
    lis2 = df2['Accession'].T.values.tolist()
    for k in rd0:
        rd01.append(lis1.index(lis[k]))
        rd02.append(lis2.index(lis[k]))
    if rd01 not in rd1 and rd02 not in rd2:
        rd1.append(rd01)
        rd2.append(rd02)


#FCP
def get_rmse(records_real, records_predict):
    if len(records_real) == len(records_predict):
        return math.sqrt(sum([(x - y) ** 2 for x, y in zip(records_real, records_predict)]) / len(records_real))
    else:
        return None


for li,lis in enumerate(rd1):
    X = df1.loc[lis, ind].T.values
    y = np.array([1 if x in indp else 0 for x in ind])
    X1 = df2.loc[rd2[li], ind1].T.values
    y1 = np.array([1 if x in indp1 else 0 for x in ind1])

    def Train_fold(X_train, y_train, X_test, y_test, l):
        Clist = [1]
        clf = LogisticRegressionCV(Cs=Clist, penalty=l, fit_intercept=False, cv=5, solver='liblinear', n_jobs=4,
                                   refit=True, class_weight='balanced', multi_class='ovr')
        clf.fit(X_train, y_train)
        result = []
        coef = clf.coef_.ravel()
        result.append(coef.tolist())
        y_test_scores = clf.predict_proba(X_test)[:, 1]
        AUC = roc_auc_score(y_test, y_test_scores)
        joblib.dump(clf, filename=str(n) + '/Protein_CC vs HC.model')
        return y_test, y_test_scores, result, AUC

    tests = np.array([])
    y_tests = np.array([])
    y_test_scores = np.array([[]])
    n='Model'
    for r in range(20):
        Skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=r)
        for i, (train, test) in enumerate(Skf.split(X, y)):
            X_train = X[train]
            y_train = y[train]
            X_test = X[test]
            y_test = y[test]
            y_test, y_test_score, result, AUC = Train_fold(X_train, y_train, X_test, y_test,'l1')
        lis1=[]
        for w in result:
            for wi, w0 in enumerate(w):
                if w0 != 0:
                    lis1.append(lis[wi])
            XX = df1.loc[lis1, ind].T.values
            Skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=r)
            y_tests, y_test_scores = [], []
            for i, (train, test) in enumerate(Skf.split(XX, y)):
                X_train = XX[train]
                y_train = y[train]
                X_test = XX[test]
                y_test = y[test]
                y_test, y_test_score, result, AUC = Train_fold(X_train, y_train, X_test, y_test,'l2')
                y_tests = y_tests + y_test
                y_test_scores = y_test_scores + X_test
                print(lis1, AUC)
            rmse = get_rmse(y_tests, y_test_scores)
        Skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=r)
        for i, (train, test) in enumerate(Skf.split(X1, y1)):
            X_train = X1[train]
            y_train = y1[train]
            X_test = X1[test]
            y_test = y1[test]
            y_test, y_test_score, result, AUC = Train_fold(X_train, y_train, X_test, y_test,'l1')
        lis1=[]
        for w in result:
            for wi, w0 in enumerate(w):
                if w0 != 0:
                    lis1.append(rd2[li][wi])
            XX1 = df2.loc[lis1, ind1].T.values
            Skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=r)
            y_tests, y_test_scores=[],[]
            for i, (train, test) in enumerate(Skf.split(XX1, y1)):
                X_train = XX1[train]
                y_train = y1[train]
                X_test = XX1[test]
                y_test = y1[test]
                y_test, y_test_score, result, AUC = Train_fold(X_train, y_train, X_test, y_test,'l2')
                y_tests=y_tests+y_test
                y_test_scores=y_test_scores+X_test
                print(lis1, AUC)
            rmse1 = get_rmse(y_tests, y_test_scores)


















