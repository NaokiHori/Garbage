import sys
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from matplotlib import pyplot as plt


lx = 1.5
ly = 1.
print("domain size (x, y): ({}, {})".format(lx, ly))
print("NOTE: x is inner row")

itot0 = 51
jtot0 = 21
dx0 = lx/itot0
dy0 = ly/jtot0
print("original domain: {} by {}".format(itot0, jtot0))

itot1 = 101
jtot1 =  51
dx1 = lx/itot1
dy1 = ly/jtot1
print("interpolated domain: {} by {}".format(itot1, jtot1))

xs0 = np.linspace(0., lx, itot0, endpoint=True)
ys0 = np.linspace(0., ly, jtot0, endpoint=True)
mgxs0, mgys0 = np.meshgrid(xs0, ys0)
data0 = np.zeros((jtot0, itot0))
for j in range(jtot0):
    for i in range(itot0):
        data0[j, i] = np.cos(2.*np.pi*xs0[i])*np.sin(2.*np.pi*ys0[j])

xs1 = np.linspace(0., lx, itot1, endpoint=True)
ys1 = np.linspace(0., ly, jtot1, endpoint=True)
mgxs1, mgys1 = np.meshgrid(xs1, ys1)
interp_func = RegularGridInterpolator((ys0, xs0), data0, method="linear", bounds_error=False, fill_value=None)
xs1d = np.ravel(mgxs1)
ys1d = np.ravel(mgys1)
grid = np.vstack([ys1d, xs1d]).T
data1 = interp_func(grid)
data1 = np.reshape(data1, (jtot1, itot1))

fig = plt.figure()
ax0 = fig.add_subplot(121)
ax1 = fig.add_subplot(122)
ax0.set_aspect(1.)
ax1.set_aspect(1.)
ax0.pcolor(mgxs0, mgys0, data0, shading="auto")
ax1.pcolor(mgxs1, mgys1, data1, shading="auto")
plt.show()
plt.close()

