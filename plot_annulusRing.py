# %% codecell
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import sasmodels.compare
import sasmodels.core
import sasmodels.direct_model
import sasmodels.data
import rheoSANS_fitOpt_Functions as rsf
import annular_sector_extraction_loop as ansect
import itertools
# %% codecell
sans = rsf.sans2d()
sans.nq = 500
sans.q = np.linspace(-sans.qmax, sans.qmax, sans.nq)
sans.xxq, sans.yyq = np.meshgrid(sans.q, sans.q)
sans.getSim('5')
simInt, a = sans.interpData_noMask(None, sans.simImport)

x, y, I, I_ann = ansect.annular(simInt, 0.07, 0.01)

x2, y2, I2 = ansect.sector(simInt, 0, np.pi/20)


# %% codecell
an = sasmodels.data.Data2D(x=x, y=y, z=I)
sans.sasPlot(an)


# %% codecell

simInt.data[~I_ann] = np.nan


a = np.reshape(simInt.qx_data, (sans.nq, sans.nq))
b = np.reshape(simInt.qy_data, (sans.nq, sans.nq))
c = np.reshape(simInt.data, (sans.nq, sans.nq))

# d = np.reshape(I_ann, (50,50))
# ax = plt.subplot()
plt.figure(figsize=[5, 5], dpi=200)
plt.plot(np.array([0, 0.06]), [0, 0.06])
plt.pcolormesh(a, b, c, cmap='jet')
plt.xlim(-0.08, 0.08)
plt.ylim(-0.08, 0.08)
plt.xticks([-0.06, -0.03, 0, 0.03, 0.06], [-60, -30, 0, 30, 60])
plt.yticks([-0.06, -0.03, 0, 0.03, 0.06], [-60, -30, 0, 30, 60])
plt.xlabel('q/10^3 A')
plt.ylabel('q/10^3 A')
plt.rc('figure', frameon=False)
plt.rc('axes.spines', top=False)
plt.rc('axes.spines', right=False)

# plt.scatter([-30,60],[30,-60])

# %% codecell

# %% codecell
simSet = set(itertools.product(self.qx_unique, self.qy_unique))
expSet = {(vals[0], vals[1]) for vals in self.expData_sort}

setDif = simSet - expSet  # (x,y) pairs in simulation not in experiment
# print(setDif)
arrayDif = np.array(list(setDif)
# print(arrayDif)
ranI=1000
# (x,y) pairs in simulation not in experiment with constant z
self.arrayDif_z=np.insert(arrayDif, 2, ranI, axis=1)
# print(arrayDif_z)
exp_data=self.expData_sort[:, 0:3]
# print(exp_data)
self.expData_fill=np.vstack((exp_data, self.arrayDif_z))
self.expData_fill_sort=np.array(
    sorted(self.expData_fill, key=lambda col: (col[1], col[0])))


zz=np.full_like(xx, np.nan)
zz[~(xx < 0)]=10
yy
np.shape(I)
plt.pcolormesh(xx, yy, zz, cmap='jet', vmin=2,
               vmax=70, norm=mcolors.LogNorm())
# %% codecell
