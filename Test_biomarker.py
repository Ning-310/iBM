import os
import numpy as np
import pandas as pd
import joblib

model=joblib.load("Model/Protein_CC vs HC.model")
df = pd.read_excel('Data/DEPs.xlsx','protein biomarker combination')
indp = 'CC1	CC2	CC3	CC4	CC5	CC6	CC7	CC8	CC9	CC10	CC11	CC12	CC13	CC14	CC15	CC16	CC17	CC18'.split('\t')
indn = 'HC1	HC2	HC3	HC4	HC5	HC6	HC7	HC8	HC9	HC10	HC11	HC12'.split('\t')
ind = indp + indn
lis=[0,1,2]
X_test = df.loc[lis, ind].T.values
y_test_scores = model.predict_proba(X_test)
print(y_test_scores)










