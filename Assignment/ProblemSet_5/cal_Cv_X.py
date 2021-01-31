import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Data files file name string
methods = ["metropolis_thermalized", "heatbath_thermalized", "singlecluster_thermalized"]
column = ["_e_errorJK", "_E2_errorJK", "_m_errorJK", "_M2_errorJK"]
filetype = ".txt"

T = 2.26
kB = 1.0


for i in range(3):
    # Read data e, m, E2, M2
    # e   -> energy density
    # m   -> magnetization density
    # E2  -> total energy square
    # M2  -> total magnetization square
    data = pd.read_csv(methods[i] + filetype, sep=' ', header=0)
    E = 10000 * np.asarray(data["e"])
    M = 10000 * np.asarray(data["m"])
    E2 = np.asarray(data["E2"])
    M2 = np.asarray(data["M2"])

    # Compute average
    E = np.average(E)
    M = np.average(M)
    E2 = np.average(E2)
    M2 = np.average(M2)

    # Compute Cv, X
    Cv = (1.0 / (kB * (T ** 2))) * (E2 - E ** 2)
    X = (1.0 / (kB * T)) * (M2 - M ** 2)

    # Print out the result
    print(methods[i])
    print("Cv")
    print(Cv)
    print("X")
    print(X)


# Result
# metropolis_thermalized
# Cv
# 17746.77482746587
# X
# 226623.13079665194
# heatbath_thermalized
# Cv
# 31865.75977493635
# X
# 2855191.510615729
# singlecluster_thermalized
# Cv
# 21420.38337121379
# X
# 685925.9055736846
#

# Calculate Error

for i in range(3):
        # Cv
    # Read data
    data = pd.read_csv(methods[i] + column[0] + filetype, sep=' ', header=0)
    ave_E = 10000 * np.asarray(data["Bjk_ave"])[0]
    error_E = 10000 * np.max(np.asarray(data["errorJK"]))

    data = pd.read_csv(methods[i] + column[1] + filetype, sep=' ', header=0)
    error_E2 = np.max(np.asarray(data["errorJK"]))

    # Find error of Cv
    error_Cv = (1.0 / (kB * (T ** 2))) * (error_E2 - 2.0 * ave_E * error_E)
    print(methods[i] + "  Cv error")
    print(error_Cv)

    # X
    # Read data
    data = pd.read_csv(methods[i] + column[2] + filetype, sep=' ', header=0)
    ave_M = 10000 * np.asarray(data["Bjk_ave"])[0]
    error_M = 10000 * np.max(np.asarray(data["errorJK"]))

    data = pd.read_csv(methods[i] + column[3] + filetype, sep=' ', header=0)
    error_M2 = np.max(np.asarray(data["errorJK"]))

    # Find error of Cv
    error_X = (1.0 / (kB * T)) * (error_M2 - 2.0 * ave_M * error_M)
    print(methods[i] + "  X error")
    print(error_X)

# Result
# metropolis_thermalized  Cv error
# 206153.74226251084
# metropolis_thermalized  X error
# -37475.6620566372
# heatbath_thermalized  Cv error
# 1585779.4013391812
# heatbath_thermalized  X error
# 7244563.805114159
# singlecluster_thermalized  Cv error
# 139921.11088573895
# singlecluster_thermalized  X error
# 144573.52458407084
