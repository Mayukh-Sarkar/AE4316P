# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 10:12:33 2021

@author: Mayukh
"""

import numpy as np
import scipy.io
from matplotlib import pyplot as plt

rms_e = scipy.io.loadmat('rms_m_e.mat')

yfit_me = scipy.io.loadmat('yfit_m_e.mat')
rms_m_e = rms_e['RMS_M_e']
y_m_e = yfit_me['yfit_e']

x = np.linspace(1,60,num = 60).reshape(1,60)
plt.figure(1)
#plt.boxplot(rms_m_e)
plt.plot(x.tolist(),y_m_e.tolist())
plt.show()