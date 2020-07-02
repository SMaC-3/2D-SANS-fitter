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
# %% codecell
sans = rsf.sans2d()
sans.nq = 50
sans.q = np.linspace(-sans.qmax, sans.qmax, sans.nq)
sans.xxq, sans.yyq = np.meshgrid(sans.q, sans.q)
sans.getSim('5')
simInt, a = sans.interpData(None, sans.simImport)
len(simInt.data)
x, y, I = ansect.annular(simInt, 0.07, 0.01)
len(x)
# %% codecell
plt.errorbar(b, av, err, marker='', lineStyle='')
plt.loglog(b, av, marker='')
# %% codecell
an = sasmodels.data.Data2D(x=x, y=y, z=I)
sans.sasPlot(an)
# dir(an)
len(an.data)
# %% codecell
plt.pcolor(an.qx_data, an.qy_data, an.data, cmap='jet', vmin=2, vmax=70, norm=mcolors.LogNorm())
