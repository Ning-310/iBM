import numpy as np
import pandas as pd
import scipy.stats as st
import math

df = pd.read_excel('Data/DEPs.xlsx','575')

dict={"CC":'CC1	CC2	CC3	CC4	CC5	CC6	CC7	CC8	CC9	CC10	CC11	CC12	CC13	CC14	CC15	CC16	CC17	CC18'.split('\t'),
      "HC":'HC1	HC2	HC3	HC4	HC5	HC6	HC7	HC8	HC9	HC10	HC11	HC12'.split('\t')
      }
tit=['CC/HC']
tim=["CC","HC"]
w=open('T-test.txt','w')
for xi in range(df[["CC1"]].T.values.shape[1]):
    X = '\t'.join(np.array(df.loc[xi, :].T.values,dtype=str).tolist())
    indm = []
    for ttsp1 in tim:
        Xm = df[dict[ttsp1]].T.values
        indm.append(str(np.mean(Xm[:, xi])))
    rat=[]
    lg=[]
    pv = []
    for tt in tit:
        ttsp = tt.split('/')
        indp = dict[ttsp[0]]
        indn = dict[ttsp[1]]
        X1 = df[indp].T.values
        X2 = df[indn].T.values
        rat.append(str(np.mean(X1[:,xi])/np.mean(X2[:,xi])))
        lg.append(str(math.log(np.mean(X1[:, xi]) / np.mean(X2[:, xi]),2)))
        p = st.ttest_ind(X1[:, xi], X2[:, xi], equal_var=False)[1]
        pv.append(str(p))
    tst = X + '\t' + '\t' + '\t'.join(indm) + '\t' + '\t'+ '\t'.join(rat)+ '\t' + '\t'+ '\t'.join(lg)+ '\t' + '\t'+ '\t'.join(pv)
    w.write(tst+'\n')

