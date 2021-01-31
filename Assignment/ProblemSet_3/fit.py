import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = pd.read_csv('output.txt', sep='\t', header=0)
N = np.asarray(data["N"])
d = np.asarray(data["d"])

fit_d = d / np.sqrt(N)

# fit with y = a * x + b


def fitting_function(x, a, b):
    return a * x + b


popt, pcov = curve_fit(fitting_function, N, fit_d)
plt.plot(N, fit_d, ".", label=r'$\frac{<d>}{\sqrt{N}}$')
plt.plot(N, fitting_function(N, *popt), "-", label="fit y = ax + b")
plt.title("Relationship of Mean Distance and Steps Moved", fontsize=16)
plt.legend(bbox_to_anchor=(1, 1))
plt.xlabel("Steps Moved", fontsize=14)
plt.ylabel(r'$\frac{<d>}{\sqrt{N}}$', fontsize=14)
plt.show()
print(popt)

#popt
#       a              b
#[1.04810971e-09 8.86492221e-01]