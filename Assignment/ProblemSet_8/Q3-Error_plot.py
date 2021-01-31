import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

folder = "Q3_result/"
Lsize = "_L200"
error_label = "_error.txt"
Algorithm = ["Jacobi", "GaussSeidel", "SOR", "CGdouble", "CGRestartDouble", "CGRestartSingle"]
# Algorithm = ["Jacobi", "GaussSeidel", "SOR", "CGRestartDouble", "CGRestartSingle"]

for alg in [Algorithm[1]]:
    data = pd.read_csv(folder + alg + Lsize + error_label, delimiter=r"\s+", header=None)
    k = np.asarray(data[0])
    error = np.asarray(data[1])
    print(k)
    print(error)

    plt.plot(k, error, 'b.-', label=alg + " Error")

plt.title("Errors in Different Algorithm", fontsize=16)
plt.xlabel("Iterations", fontsize=14)
plt.ylabel("Errors", fontsize=14)
plt.legend()
plt.show()
