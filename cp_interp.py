import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.interpolate import interp1d, UnivariateSpline, griddata
import numpy as np

rc('text', usetex=True)

cp = pd.read_csv('cp.csv')


angles = [15, 20, 25, 30, 35, 40, 45]
N = 100

points = []
values = []

for angle in angles:
    x = cp[str(angle)+'X'].dropna().values
    y = cp[str(angle)+'Y'].dropna().values

    indx = x.argsort()
    x = x[indx]
    y = y[indx]

    xp = np.linspace(min(x), max(x), N)
    f = UnivariateSpline(x, y, k=5)

    points.extend([(xi, yi) for xi, yi in zip(xp, f(xp))])
    values.extend([angle] * N)

grid_x_ = np.linspace(0, 3, 100)
grid_y_ = np.linspace(0, 0.4, 100)

points = np.array(points)
values = np.array(values)

grid_x, grid_y = np.meshgrid(grid_x_, grid_y_)

grid_z = griddata(points, values, (grid_x, grid_y), method='cubic')

J = 1.136
Cp = 0.293

beta = griddata(points, values, (J, Cp), method='cubic')

im = plt.contourf(grid_x, grid_y, grid_z, levels=20)

plt.scatter(J, Cp, s=100, marker='*', color='r')
cbar = plt.colorbar(im)
cbar.set_label(r'$\beta$ (deg)')

plt.xlabel('$J$')
plt.ylabel('$C_P$')
plt.savefig('cp_interp.pdf')
# plt.show()
print(beta)
