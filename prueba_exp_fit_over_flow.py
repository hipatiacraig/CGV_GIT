# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import numpy as np
import lmfit


ts = np.array([ 15,  32,  51, 106, 123, 142, 160, 177, 196, 213, 232, 249, 269, 286, 323, 340, 359, 375, 394, 466, 484, 520, 539, 645, 681])
amps = np.array([78, 64, 64, 42, 42, 15, 34, 29, 34, 31, 31, 22,  5,  6,  8,  4, 11, 14, 14,  1,  2, 10,  4,  3,  1])
emodel = lmfit.Model(lambda x,t,A: A*np.exp(-x/t))

plt.errorbar(ts, amps, xerr=2, yerr=np.sqrt(amps), fmt="ko-", capsize = 5)
plt.plot(ts, emodel.fit(amps, x=ts, t=150, A=140).best_fit, 'r-', label='best fit')
plt.plot(ts, emodel.fit(amps, x=ts, weights=1/np.sqrt(amps), t=150, A=140).best_fit, 'r--', label='weighted best fit (1/err)')
plt.plot(ts, emodel.fit(amps, x=ts, weights=1/amps, t=150, A=140).best_fit, 'r:', label='weighted best fit (1/err²)')
plt.legend()