
import scipy
import numpy as np
import matplotlib.pyplot as pl
import scipy.stats as st
#number of observations could matter here
data = np.random.multivariate_normal((1., 1), [[1, .478], [.478, 1]], 100000)

x = data[:, 0]
y = data[:, 1]
xmin, xmax = -4, 3
ymin, ymax = -3, 3.2

# Peform the kernel density estimate
xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([x, y])
kernel = st.gaussian_kde(values)
f = np.reshape(kernel(positions).T, xx.shape)

fig = pl.figure()
ax = fig.gca()
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
# Contourf plot
cfset = ax.contourf(xx, yy, f,cmap='Blues')#colors='white')
## Or kernel density estimate plot instead of the contourf plot
ax.imshow(np.rot90(f), cmap='Blues', extent=[xmin, xmax, ymin, ymax])

# Contour plot
cset = ax.contour(xx, yy, f, colors='k', linewidths=0.5)
# Label plot
gamma=-0.75
ax.clabel(cset, inline=1, fontsize=10)
ax.set_xlabel('Bankruptcy Score')
ax.set_ylabel('Going Concern Score')
ax.hlines(0, 0, 6, lw=1, color='black')
ax.hlines(0, 0, -6, lw=1, color='black')
ax.vlines(0, 0, 6, lw=1, color='black')
ax.vlines(0, 0, -6, lw=1, color='black')
ax.vlines(gamma, 0, 4, lw=1, color='black')
ax.text(-2, -2, '$\pi_\mathrm{00}$', fontsize=12)
ax.text(2, -2, '$\pi_\mathrm{10}$', fontsize=12)
ax.text(-2, 2, '$\pi_\mathrm{01}$', fontsize=12)
ax.text(2.132, 2.018, '$\pi_\mathrm{11}$',  fontsize=12)

ax.text(-.4, 2.5, 'A', fontsize=12)#

ax.text(gamma-0.04, -0.335, '$-\gamma$', fontsize=12)
import matplotlib.transforms as mtransforms
trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
y1positive=y>0
d = scipy.zeros(len(y))

ax.fill_betweenx(y, gamma, 0, where= y>=d, 
                facecolor='grey', alpha=0.5)
#
pl.ylabel('Exposure score: $\Phi^{-1}(P(G\mid \mathbf{x}))$')
pl.xlabel('Outcome score: $\Phi^{-1}(P(B\mid \mathbf{x}))$')
pl.savefig('C:/Users/demetri/Documents/research/comp/bivar_fit_2_comp_bank.jpg')
#pl.savefig('C:/Users/demetri/Documents/research/comp/bivar_fit_2.pdf')
pl.show()