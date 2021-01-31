import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Read data
data = pd.read_csv('ave_E-T.txt', sep=' ', header=0)
ave_E = np.asarray(data["ave_E"])
T = np.asarray(data["T"])

# fit with y = a * x + b
def fitting_function(x, a, b):
    return a * x + b


# Draw data
plt.scatter(T[0], ave_E[0], color='r', label="abandon data")
plt.scatter(T[1:], ave_E[1:], color='b', label="data")
plt.xlabel('Temperature T', fontsize=14)
plt.ylabel('Energy Density <e>', fontsize=14)
plt.title('<e> - T', fontsize=16)
# Draw fit data
popt, pcov = curve_fit(fitting_function, T[1:], ave_E[1:])
plt.plot(T[1:], fitting_function(T[1:], *popt), 'g--', label="fit y = ax + b")
plt.legend()
plt.show()

print(popt)
# [     a            b    ]
# [ 1.30615079 -4.39145985]
