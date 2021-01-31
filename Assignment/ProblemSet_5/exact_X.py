import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

# Read data
data = pd.read_csv('ave_M-B.txt', sep=' ', header=0)
ave_M = np.asarray(data["ave_M"])
error = np.asarray(data["error"])
B = np.asarray(data["B"])

# fit with y = a * x^2 + b * x + c, consider error


def model(x, a, b, c):
    value = 0.
    value = a * np.square(x) + b * x + c
    return value


def fcn(p):
    expt = model(B[:-1], p[0], p[1], p[2])
    delta = (ave_M[:-1] - expt) / error[:-1]
    return (delta ** 2).sum()


p_init = np.array([14, 0.0, 0.0])
r = opt.minimize(fcn, p_init)
fitting_output = r.x
fitting_x = np.linspace(np.min(B[:-1]), np.max(B[:-1]), 10000, endpoint=True)

print(fitting_output)

# Draw data
plt.errorbar(np.log10(B[:-1]), ave_M[:-1], yerr=error[:-1], fmt='.b', label='Data')
plt.errorbar(np.log10(B[-1]), ave_M[-1], yerr=error[-1], fmt='.r', label='Abandob Data')
plt.ylabel('Magnetization Density <m>', fontsize=14)
plt.xlabel('External Magnetic Field ' + r'$log_{10}(B)$', fontsize=14)
plt.title('<m> - B', fontsize=16)

# Draw fit data
plt.plot(np.log10(fitting_x), np.polyval(fitting_output, fitting_x), 'g--', label="fit y = a * x^2 + b * x + c")
plt.legend()
plt.show()


# Fitting with first 4 data
# [     a           b            c    ]
# [14.02700005 31.05066151  0.63445071]

# Fitting with all data.
# Ignore this, since the last data is too far from zero
# [     a           b            c    ]
# [-2.24966442e+03  3.44880025e+01  6.33879727e-01]
