import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


file_name = "Q1_result.txt"
data = pd.read_csv(file_name, sep=' ', header=None)
step = np.asarray(data[0])
x = np.asarray(data[1])
exact_sol = np.asarray(data[2])
euler = np.asarray(data[3])
modified_euler = np.asarray(data[4])
improved_euler = np.asarray(data[5])
rungekutta3 = np.asarray(data[6])

plt.plot(x, np.abs(exact_sol - euler), label='Euler')
plt.plot(x, np.abs(exact_sol - modified_euler), label='Modified Euler')
plt.plot(x, np.abs(exact_sol - improved_euler), label='Improved Euler')
plt.plot(x, np.abs(exact_sol - rungekutta3), label='Runge Kutta 3rd')
plt.title("Errors For Each Algorithm", fontsize=14)
plt.xlabel("x (step made)", fontsize=14)
plt.ylabel("Errors", fontsize=14)
plt.legend()
plt.show()
