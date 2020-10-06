import sys
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import h5py as h5
from matplotlib import pyplot as plt
from tqdm import trange


def cordin(which):
    with h5.File("{}/continua_master.h5".format(which), "r") as f:
        ylen = f["ylen"][()][0]
        zlen = f["zlen"][()][0]
    with h5.File("{}/cordin_info.h5".format(which), "r") as f:
        xm = f["xm"][()]
        xc = f["xc"][()]
        ym = f["ym"][()]
        zm = f["zm"][()]
    jtot = ym.shape[0]
    ktot = zm.shape[0]
    yc = list()
    yc.append(0.)
    for j in range(jtot-1):
        yc.append(0.5*(ym[j]+ym[j+1]))
    yc.append(ylen)
    yc = np.array(yc)
    zc = list()
    zc.append(0.)
    for k in range(ktot-1):
        zc.append(0.5*(zm[k]+zm[k+1]))
    zc.append(zlen)
    zc = np.array(zc)
    return xm, xc, ym, yc, zm, zc

def master():
    xm,  xc,  ym,  yc,  zm,  zc  = cordin("original")
    rxm, rxc, rym, ryc, rzm, rzc = cordin("refined")
    return xm, xc, ym, yc, zm, zc, rxm, rxc, rym, ryc, rzm, rzc

def ndmesh(*args):
    # https://stackoverflow.com/a/16409517
    args = map(np.asarray,args)
    return np.broadcast_arrays(*[x[(slice(None),)+(None,)*i] for i, x in enumerate(args)])

def interp(xs0, ys0, zs0, xs1, ys1, zs1, varname):
    itot0 = xs0.shape[0]
    jtot0 = ys0.shape[0]
    ktot0 = zs0.shape[0]
    itot1 = xs1.shape[0]
    jtot1 = ys1.shape[0]
    ktot1 = zs1.shape[0]
    with h5.File("original/continua_{}.h5".format(varname), "r") as f:
        data0 = f["var"][()]
    print("argument sizes of RegularGridInterpolator: {} {} {} {}".format(zs0.shape, ys0.shape, xs0.shape, data0.shape))
    interp_func = RegularGridInterpolator((zs0, ys0, xs0), data0, method="linear", bounds_error=False, fill_value=None)
    mgxs1, mgys1, mgzs1 = ndmesh(xs1, ys1, zs1)
    mgxs1 = np.ravel(mgxs1)
    mgys1 = np.ravel(mgys1)
    mgzs1 = np.ravel(mgzs1)
    grid1d = np.vstack([mgzs1, mgys1, mgxs1]).T
    # interpolation main
    print("interpolating...")
    data1 = interp_func(grid1d)
    print("done")
    data1 = np.reshape(data1, (ktot1, jtot1, itot1))
    with h5.File("refined/continua_{}.h5".format(varname), "w") as f:
        f.create_dataset("var", data=data1)
    ## output
    # mgxs0, mgys0 = np.meshgrid(xs0, ys0)
    # mgxs1, mgys1 = np.meshgrid(xs1, ys1)
    # fig = plt.figure()
    # ax0 = fig.add_subplot(121)
    # ax1 = fig.add_subplot(122)
    # ax0.set_aspect(1.)
    # ax1.set_aspect(1.)
    # ax0.pcolor(mgxs0, mgys0, data0[10,:,:], shading="auto")
    # ax1.pcolor(mgxs1, mgys1, data1[10,:,:], shading="auto")
    # plt.show()
    # plt.close()


if __name__ == "__main__":
    xm, xc, ym, yc, zm, zc, rxm, rxc, rym, ryc, rzm, rzc = master()
    interp(xc, ym, zm, rxc, rym, rzm, "vx")
    interp(xm, yc, zm, rxm, ryc, rzm, "vy")
    interp(xm, ym, zc, rxm, rym, rzc, "vz")

