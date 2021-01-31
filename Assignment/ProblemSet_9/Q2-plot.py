import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

folder = "HMC/"
data = pd.read_csv(folder + "HMC_x.txt", header=None)
x = np.asarray(data[0])
data = pd.read_csv(folder + "HMC_Prop.txt", header=None, sep=' ')
P_ave = np.asarray(data[0])
P_err = np.asarray(data[1])
data = pd.read_csv(folder + "HMC_Sub.txt", header=None, sep=' ')
S_ave = np.asarray(data[0])
S_err = np.asarray(data[1])

data = pd.read_csv(folder + "exact.txt", header=None, sep=' ')
exact_P = np.asarray(data[1])
exact_S = np.asarray(data[0])

sortindex = np.argsort(x)
x = np.sort(x)
P_ave = P_ave[sortindex]
P_err = P_err[sortindex]
S_ave = S_ave[sortindex]
S_err = S_err[sortindex]

plt.errorbar(x, P_ave, P_err, label='Propagator')
plt.plot(x, exact_P, label='Exact Solution')
plt.title("Amplitude of the Propagator with HMC", fontsize=16)
plt.xlabel("x", fontsize=14)
plt.ylabel("<P>", fontsize=14)
plt.legend()
plt.show()

plt.errorbar(x, S_ave, S_err, label='Subpropagator')
plt.plot(x, exact_S, label='Exact Solution')
plt.title("Amplitude of the Subpropagator with HMC", fontsize=16)
plt.xlabel("x", fontsize=14)
plt.ylabel("<S>", fontsize=14)
plt.legend()
plt.show()
