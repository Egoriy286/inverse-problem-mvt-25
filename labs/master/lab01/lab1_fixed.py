import numpy as np; import matplotlib.pyplot as plt;
L=1; T=0.1; J=100; N=100; mu1=0; mu2=0
tau=T/J; h=L/N; r=2*h**2/tau; c=2+r
u=np.zeros((J+1, N+1)); x=np.linspace(0, L, N+1); t=np.linspace(0, T, J+1)
u[0,:]=np.sin(np.pi * x); p_e = np.exp(2*t)*np.sin(np.pi/4*t); phi = np.zeros(J+1)
for j in range(1, J+1):
    A=-1.; B=c-h**2*p_e[j-1]; C=-1.
    xi=np.zeros(N+1); eta=np.zeros(N+1); xi[1]=0; eta[1]=mu1
    for i in range(1, N):
        xi[i+1]= -C / (A*xi[i] + B)
        eta[i+1]= (r*u[j-1][i]+(u[j-1][i+1]-2.*u[j-1][i]+u[j-1][i-1] + h**2*p_e[j]*u[j-1][i]) - A*eta[i]) / (A*xi[i] + B)
    u[j][N]=mu2
    for i in range(N - 1, -1, -1): u[j][i]=xi[i+1]*u[j][i+1]+eta[i+1]
    u[j][0]=mu1

x_ = 0.5
n = int(round(x_ / h))
phi = u[:, n]

y=np.zeros((J+1, N+1)); z=np.zeros(N+1); v=np.zeros(N+1); p=np.zeros(J+1)
y[0,:]=u[0,:]
p[0]=p_e[0] #p_e[0] # p[0]=p_e[0] # точность 10^(-12) # p[0]=p0_est # точность 10^(-4) 
for j in range(1,J+1):
    A=-1.; B=c-h**2*p[j-1]; C=-1.; z[N]=mu2; v[N]=0
    xi_z=np.zeros(N+1); eta_z=np.zeros(N+1); xi_z[1]=0; eta_z[1]=mu1
    xi_v=np.zeros(N+1); eta_v=np.zeros(N+1); xi_v[1]=0;   eta_v[1]=0
    for i in range(1, N):
        xi_z[i+1]= -C / (A*xi_z[i] + B)
        eta_z[i+1]= (r*y[j-1][i]+(y[j-1][i+1]-2.*y[j-1][i]+y[j-1][i-1]) - A*eta_z[i]) / (A*xi_z[i] + B)
        xi_v[i+1]= -C / (A*xi_v[i] + B)
        eta_v[i+1]= ( h**2*y[j-1][i] - A*eta_v[i]) / (A*xi_v[i] + B)
    
    for i in range(N - 1, -1, -1): 
        z[i]=xi_z[i+1]*z[i+1]+eta_z[i+1]; v[i]=xi_v[i+1]*v[i+1]+eta_v[i+1]
    if abs(v[n]) > 1e-12:
        p[j] = (phi[j] - z[n]) / v[n]
    else:
        p[j] = p[j-1]
    y[j,:] = z + p[j]*v; y[j,0]=mu1; y[j, N]=mu2

# Визуализация решения
im = plt.imshow(y)
plt.colorbar(im)
plt.title("Приближенное решение")

plt.figure()
im = plt.imshow(u)
plt.colorbar(im)
plt.title("Точное решение")

plt.figure()
plt.plot(p_e, label='exact')
plt.plot(p, label='aprox')
plt.legend()
plt.title("p(t)=")
plt.grid()

plt.figure()
plt.plot(abs(p-p_e))
plt.grid()
plt.title("Погрешность восстановления p-p_e")

plt.show()
