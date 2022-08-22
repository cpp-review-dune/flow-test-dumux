#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

x = np.array([58.2, 39.43, 36.95, 34.48, 28.30])
y = np.array([0.1, 1, 3, 5, 10])

plt.clf()
plt.plot(x, y, marker="o", color="b", linewidth=2)
plt.savefig("plot.pdf")
