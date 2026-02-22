import numpy as np; import matplotlib.pyplot as plt;
L=1; T=1; J=100; N=100; mu1=0; mu2=0
tau=T/J; h=L/N; r=2*h**2/tau
y=np.zeros((J+1, N+1)); x=np.linspace(0, L, N+1); t=np.linspace(0, T, J+1)
y[0,:]=np.sin(np.pi * x); v_e = np.sin(t)
for j in range(1,J):
    A=-1.; B=r+2.; C=-1.
    ksi=np.zeros(N+1); eta=np.zeros(N+1); ksi[0]=mu1; eta[0]=0
    for i in range(1, N-1):
        ksi[i+1]= -C / (A*ksi[i] + B)
        eta[i+1]= (r*y[j-1,i]+y[j-1,i+1]-2*y[j-1,i]+y[j-1,i-1]-v_e[j]*2*h*(y[j-1][i] - y[j-1,i-1]) - A*eta[i]) / (A*ksi[i] + B)
        y[j][N]=mu2
    for i in range(N - 1, 0, -1): y[j][i]=ksi[i+1]*y[j][i+1]+eta[i+1]

p = plt.imshow(y)
plt.colorbar(p)
plt.plot()
plt.show()



