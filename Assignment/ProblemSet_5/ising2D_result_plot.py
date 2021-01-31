import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Read data e, m, E2, M2
# e   -> energy density
# m   -> magnetization density
# E2  -> total energy square
# M2  -> total magnetization square
data = pd.read_csv('heatbath_thermalized.txt', sep=' ', header=0)
swp = np.asarray(data["sweep"])
e = np.asarray(data["e"])
m = np.asarray(data["m"])
E2 = np.asarray(data["E2"])
M2 = np.asarray(data["M2"])

print(swp)
print(E2)
print(M2)

data = pd.read_csv('heatbath_thermalized_e_errorJK.txt', sep=' ', header=0)
nb = np.asarray(data["nb"])
ave = np.asarray(data["Bjk_ave"])
errorJK = np.asarray(data["errorJK"])

print(nb)
print(ave)
print(errorJK)
print(np.max(errorJK))

# nb = np.log2(nb)

# print(nb)
# print(errorJK)

# fit with y = a * x + b
# def fitting_function(x, a, b):
#     return a * x + b


# plt.scatter(nb, errorJK)
# plt.show()

# popt, pcov = curve_fit(fitting_function, N, fit_d)
# plt.plot(N, fit_d, ".", label=r'$\frac{<d>}{\sqrt{N}}$')
# plt.plot(N, fitting_function(N, *popt), "-", label="fit y = ax + b")
# plt.title("Relationship of Mean Distance and Steps Moved", fontsize=16)
# plt.legend(bbox_to_anchor=(1, 1))
# plt.xlabel("Steps Moved", fontsize=14)
# plt.ylabel(r'$\frac{<d>}{\sqrt{N}}$', fontsize=14)
# plt.show()
# print(popt)
