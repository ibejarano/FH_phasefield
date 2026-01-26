import matplotlib.pyplot as plt
import numpy as np

#data1 = np.loadtxt("output_job_1/output_refactored.csv", delimiter=",", skiprows=1)
data1 = np.loadtxt("output_kres_sym/results.csv", delimiter=",", skiprows=1)

ts1, ps1 = data1[:, 0], data1[:, 1]

data2 = np.loadtxt("output_adapt_sym/results.csv", delimiter=",", skiprows=1)

ts2, ps2 = data2[:, 0], data2[:, 1]

plt.plot(ts1, ps1, label="kres")
plt.plot(ts2, ps2, marker="o", label="dt adapt")

plt.plot(ts1, 0.9e4*ts1**(-1/3), "k--", alpha=0.5)
plt.legend()
plt.xscale("log")
plt.yscale("log")
plt.show()


wp1, wm1 = data1[:, 3], data1[:, 4]
op1 = wp1-wm1

wp2, wm2 = data2[:, 3], data2[:, 4]
op2 = wp2-wm2

plt.plot(ts1, op1, label="full")
plt.plot(ts2, op2, label="simetrico")

#plt.plot(ts2, 0.9e4*ts2**(-1/3), "k--", alpha=0.5)
plt.legend()
plt.xscale("log")
plt.yscale("log")
plt.show()