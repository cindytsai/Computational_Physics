import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

col_names = ["N", "SS", "sign1", "SSerror", "RM", "sign2", "RMerror", "MA", "sign3", "MAerror"]
data = pd.read_csv('integration.txt', sep=' ', names=col_names, header=0)

N = np.asarray(data["N"])
N_log = np.log2(N)
SS = np.asarray(data["SS"])
SSerror = np.asarray(data["SSerror"])
RM = np.asarray(data["RM"])
RMerror = np.asarray(data["RMerror"])
MA = np.asarray(data["MA"])
MAerror = np.asarray(data["MAerror"])


# plt.subplot(3, 1, 1)
# plt.errorbar(N_log, SS, yerr=SSerror, fmt='.b', label='Simple Sampling')
# plt.legend(bbox_to_anchor=(1, 1), fontsize=12)
# plt.title("High Dimensional Integration with Different Method", fontsize=16)

# plt.subplot(3, 1, 2)
# plt.errorbar(N_log, RM, yerr=RMerror, fmt='.r', label='Rejection Method')
# plt.legend(bbox_to_anchor=(1, 1), fontsize=12)

# plt.subplot(3, 1, 3)
# plt.errorbar(N_log, MA, yerr=MAerror, fmt='.g', label='Metropolis Algorithm')
# plt.legend(bbox_to_anchor=(1, 1), fontsize=12)
# plt.xlabel("Sample Points  " + r'$2^n$', fontsize=14)

# plt.show()

print(data[["N", "SS", "SSerror", "RM", "RMerror", "MA", "MAerror"]])
