import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

algorithm = "SingleCluster"

file_name = "./Singlecluster_Auto/" + algorithm + "_Corr_L_tauA_E.txt"
data = pd.read_csv(file_name, sep=' ', header=None)
tauA_E = np.asarray(data[1])

file_name = "./Singlecluster_Auto/" + algorithm + "_Corr_L_tauA_M.txt"
data = pd.read_csv(file_name, sep=' ', header=None)
tauA_M = np.asarray(data[1])

file_name = "./Singlecluster_Auto/" + algorithm + "_Corr_L.txt"
data = pd.read_csv(file_name, header=None)
L = np.asarray(data[0])


# fit with y = a * x + b
def fitting_function(x, a, b):
    return a * x + b


# Fit for auto correlation from <m>
popt, pcov = curve_fit(fitting_function, np.log10(L), np.log10(tauA_M))
print(popt)
# Result
#
# Metropolis
#        a           b
# [ 2.05188199 -1.52904038]
#
# Heatbath
#        a           b
# [ 1.53576901 -0.2867087 ]
#
# Singlecluster
#        a           b
# [-1.23845871  1.10069914]

x = np.linspace(np.min(np.log10(L)), np.max(np.log10(L)), 1000, endpoint=True)

plt.figure(num=None, figsize=(10, 5), dpi=80, facecolor='w', edgecolor='k')
plt.subplot(1, 2, 2)
plt.suptitle(algorithm, fontsize=16)
plt.plot(x, fitting_function(x, *popt), 'r--', label="fit y = ax + b")
plt.scatter(np.log10(L), np.log10(tauA_M), label="data from <m>")
plt.title(r'$log_{10}\tau_A$' + "-" + r'$log_{10}L$', fontsize=14)
plt.xlabel("Box Length " + r'$log_{10}L$', fontsize=14)
plt.ylabel(r'$log_{10}\tau_A$', fontsize=14)
plt.legend()
plt.subplot(1, 2, 1)
plt.scatter(np.log10(L), np.log10(tauA_E), label="data from <e>")
plt.title(r'$log_{10}\tau_A$' + "-" + r'$log_{10}L$', fontsize=14)
plt.xlabel("Box Length " + r'$log_{10}L$', fontsize=14)
plt.ylabel(r'$log_{10}\tau_A$', fontsize=14)
plt.legend()
plt.show()
