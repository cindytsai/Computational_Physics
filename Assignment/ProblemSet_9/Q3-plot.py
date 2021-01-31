import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

folder = ["HMC_lambda0_4/", "HMC_lambda1E-1/", "HMC_lambda1E-2/", "HMC_lambda1E-3/", "HMC_lambda1E-4/", "HMC_lambda1E-5/", "HMC_lambda1E-6/"]
lambda_value = [0.4, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]

# Get exact solution, lambda = 0
data = pd.read_csv("HMC/exact.txt", header=None, sep=' ')
exact_P = np.asarray(data[1])
exact_S = np.asarray(data[0])


for i in range(len(lambda_value)):

    data = pd.read_csv(folder[i] + "HMC_x.txt", header=None)
    x = np.asarray(data[0])
    data = pd.read_csv(folder[i] + "HMC_Prop.txt", header=None, sep=' ')
    P_ave = np.asarray(data[0])
    P_err = np.asarray(data[1])

    sortindex = np.argsort(x)
    x = np.sort(x)
    P_ave = P_ave[sortindex]
    P_err = P_err[sortindex]

    plt.errorbar(x[1:], P_ave[1:], yerr=P_err[1:], label=r'$\lambda$' + " = " + str(lambda_value[i]), fmt='o')


plt.plot(x, exact_P, label=r'$\lambda$' + " = 0")
plt.title("Propagator with HMC", fontsize=16)
plt.xlabel("x", fontsize=14)
plt.ylabel("<P>", fontsize=14)
plt.legend()
plt.show()


for i in range(len(lambda_value)):
    data = pd.read_csv(folder[i] + "HMC_x.txt", header=None)
    x = np.asarray(data[0])
    data = pd.read_csv(folder[i] + "HMC_Sub.txt", header=None, sep=' ')
    S_ave = np.asarray(data[0])
    S_err = np.asarray(data[1])

    sortindex = np.argsort(x)
    x = np.sort(x)
    S_ave = S_ave[sortindex]
    S_err = S_err[sortindex]

    plt.errorbar(x[1:], S_ave[1:], yerr=S_err[1:], label=r'$\lambda$' + " = " + str(lambda_value[i]), fmt='o')

plt.plot(x, exact_S, label=r'$\lambda$' + " = 0")
plt.title("Subpropagator with HMC", fontsize=16)
plt.xlabel("x", fontsize=14)
plt.ylabel("<S>", fontsize=14)
plt.legend()
plt.show()
