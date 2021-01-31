import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

folder = "Q3_result/"
Lsize = "_L200"
Algorithm = ["Jacobi", "GaussSeidel", "SOR", "CGdouble", "CGRestartDouble", "CGRestartSingle"]
# Algorithm = ["Jacobi", "GaussSeidel", "SOR", "CGRestartDouble", "CGRestartSingle"]

r = np.arange(1, 100)
phi = (1. / (4. * np.pi)) * ((1.0 / r) - 1.0)

for alg in Algorithm:
    data = pd.read_csv(folder + alg + Lsize + ".txt", delimiter=r"\s+", header=None)
    U0 = np.asarray(data[0])

    plt.plot(r, U0[1:100], label=alg)

plt.plot(r, phi, label="ExactSolution")
plt.title("Exact and Simulation Solution", fontsize=16)
plt.xlabel("r", fontsize=14)
plt.ylabel(r"$\phi(r)$", fontsize=14)
plt.legend()
plt.show()
