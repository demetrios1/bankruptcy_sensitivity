import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import matplotlib as mpl
import pandas as pd
mpl.style.use('seaborn-white')
# Define dimension. 
d = 2
# Set mean vector. 
m = np.array([1, 2]).reshape(2, 1)

# Set covariance function. 
K_0 = np.array([[2, 1],
                [1, 2]])

np.linalg.eigvals(K_0)

# Define epsilon.
epsilon = 0.0001

# Add small pertturbation. 
K = K_0 + epsilon*np.identity(d)

L = np.linalg.cholesky(K)
L

np.dot(L, np.transpose(L))

# Number of samples. 
n = 10000

u = np.random.normal(loc=0, scale=1, size=d*n).reshape(d, n)

x = m + np.dot(L, u)
B=x[0]
G=x[1]


sns.jointplot(x=B,
              y=G, 
             kind="kde", space=0);
#plt.xlabel('$P(B\mid \mathbf{x})$')
#plt.ylabel('$P(G\mid \mathbf{x})$')

#plt.savefig('C:/Users/demetri/Documents/ASU/dissertation/going_vs_bank.pdf')
plt.show()
