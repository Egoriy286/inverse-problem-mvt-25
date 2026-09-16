import numpy as np;  from scipy import sparse;from scipy.sparse import linalg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import sparse; from scipy.sparse import linalg
def u0(x,y,t): return .25*np.exp(-.25*(x**2+y**2)/t)/t
def poisson(n,a,b):
    V = []; I = []; J = []
    for j in range(n):
        for i in range(n):
            k=i+n*j;  V.append(b); I.append(k); J.append(k)
            if i>0:   V.append(a); I.append(k); J.append(k-1);
            if i<n-1: V.append(a); I.append(k); J.append(k+1)
            if j>0:   V.append(a); I.append(k); J.append(k-n)
            if j<n-1: V.append(a); I.append(k); J.append(k+n)
    return sparse.coo_matrix((V,(I,J)))
def gr(N,l,a,b,m,z):
    x=np.linspace(0.,l,N+1); y=np.linspace(0.,l,N+1); X,Y = np.meshgrid(x,y);
    p=np.linspace(0.,1.,m); q=np.linspace(a,b,m); plt.contourf(x,y,z,q); plt.colorbar()
    plt.xlabel('$x$'); plt.ylabel('$y$'); plt.show()
def gr1(N,l,z,i,j):
    x=np.linspace(0.,l,N+1); y=np.linspace(0.,l,N+1); X,Y = np.meshgrid(x,y);
    X,Y=np.meshgrid(x,y); fig=plt.figure(); plt.xlabel('$X$'); plt.ylabel('$Y$')
    plt.contourf(X,Y,z,levels=10,cmap='Spectral_r'); plt.colorbar(); plt.show()