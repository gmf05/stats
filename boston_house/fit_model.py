# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 08:38:20 2016

@author: gmf
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
#import matplotlib.pyplot as plt

# Load data
headers = open('headers.txt','rb').readline().strip().split('\t')
data = pd.read_csv('housing.csv', header=None)
data.columns = headers

data = data[data.RAD<20] # Limit to "within or near town" homes
x1 = data.INDUS # Proportion of non-retail business acres per town
x2 = data.ZN # Proportion of "big (residential) lots" 
y = data.TAX # Tax rate

# Correlation between predictors (2,3) and (1,2) : 0.518 and -0.128
print 'Correlations: %f, %f ' % (np.corrcoef(x1,y)[1,0], np.corrcoef(x2,y)[1,0])

# ZN & INDUS strongly anti-correlated, as expected: -0.470
#  since these two types of land do not usually coexist
print 'Correlations: %f' % np.corrcoef(x1,x2)[1,0]

# Fit models
X1 = sm.add_constant(x1)
X2 = sm.add_constant(x2)
Xfull = sm.add_constant(data[['INDUS', 'ZN']])

m1 = sm.OLS(y, X1).fit()
m2 = sm.OLS(y, X2).fit()
mfull = sm.OLS(y, Xfull).fit()

# Print results
print m1.rsquared, m2.rsquared # individual R^2
print m1.rsquared+m2.rsquared # sum of individual R^2
print mfull.rsquared # full model R^2 > sum