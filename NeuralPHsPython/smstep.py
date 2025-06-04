import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import comb

def cb(n, k):
    return comb(n, k, exact=True)

def sig(r):
    delta = 0.01
    return np.sqrt(r**2 + delta**2) - delta

def step(sigma, b, n):
    temp = 0
    if sigma <= b:
        for k in range(n + 1):
            temp += cb(n + k, k) * cb(2 * n + 1, n - k) * (-sigma / b) ** k
        temp *= (sigma / b) ** (n + 1)
    else:
        temp = 1.0
    return temp

def d_step(sigma, b, n):
    temp = 0
    if sigma <= b:
        for k in range(n + 1):
            temp += cb(n + k, k) * cb(2 * n + 1, n - k) * (-1) ** k * (sigma / b) ** (n + k) * (n + k + 1) / b
    return temp

range_vals = np.arange(-3, 3.01, 0.01)
h = np.zeros(len(range_vals))
dh = np.zeros(len(range_vals))

rc = 2.0
b = sig(rc)
n = 2

for i, r in enumerate(range_vals):
    h[i] = step(sig(r), b, n)
    dh[i] = d_step(sig(r), b, n)

# plt.figure(figsize=(6, 5))
# plt.plot(range_vals, h, linestyle='solid', linewidth=3, label='h')
# plt.plot(range_vals, dh, linestyle='solid', linewidth=3, label="∂h")
# plt.xlim(min(range_vals), max(range_vals))
# plt.ylim(-0.1, 1.2)
# plt.xticks([-rc, 0, rc], ["$-b$", "$0$", "$b$"])
# plt.yticks([0, 1], ["$0$", "$1$"])
# plt.annotate("$h(σ(\Vert\mathbf{r}\Vert_2))$", (-1.1, 1.1), color='red', fontsize=14)
# plt.annotate("${∂h}/{∂σ}$", (3.0, 0.2), color='red', fontsize=14)
# plt.annotate("$\Vert\mathbf{r}\Vert_2$", (2.9, -0.19), color='red', fontsize=14)
# plt.legend()
# plt.savefig(os.path.join("..","figs", "hr.png")
# plt.show()

# X, Y = np.meshgrid(range_vals, range_vals)
# Z = np.array([[step(sig(np.sqrt(r1**2 + r2**2)), b, n) for r1, r2 in zip(row_x, row_y)] for row_x, row_y in zip(X, Y)])

# fig = plt.figure(figsize=(6, 5))
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X, Y, Z, cmap=cm.viridis)
# ax.set_xticks([-rc, 0, rc])
# ax.set_xticklabels(["$-b$", "$0$", "$b$"])
# ax.set_yticks([-rc, 0, rc])
# ax.set_yticklabels(["$-b$", "$0$", "$b$"])
# ax.text(-1, -1, 1.1, "$h(σ(\mathbf{r}))$", color='red', fontsize=14)
# ax.text(1, -1, 0, "$\mathbf{r}$", color='red', fontsize=14)
# plt.savefig(os.path.join("..","figs", "3Dhr.png")
# plt.show()
