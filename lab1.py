import numpy as np; import matplotlib.pyplot as plt;
L=1; T=1; J=100; N=100; mu1=0; mu2=0
tau=T/J; h=L/N; r=2*h**2/tau; c=2+r; e=2-r
u=np.zeros((J+1, N+1)); x=np.linspace(0, L, N+1); t=np.linspace(0, T, J+1)
u[0,:]=np.sin(np.pi * x); p_e = np.sin(t); phi = np.zeros(J+1)
for j in range(1,J):
    A=-1.; B=2+r-h**2*p_e[j-1]; C=-1.
    ksi=np.zeros(N+1); eta=np.zeros(N+1); ksi[0]=mu1; eta[0]=0
    for i in range(1, N-1):
        ksi[i+1]= -C / (A*ksi[i] + B)
        eta[i+1]= (r*u[j-1][i]+(u[j-1][i+1]-2.*u[j-1][i]+u[j-1][i-1] + h**2*p_e[j]*u[j-1][i]) - A*eta[i]) / (A*ksi[i] + B)
        u[j][N]=mu2
    for i in range(N - 1, 0, -1): u[j][i]=ksi[i+1]*u[j][i+1]+eta[i+1]
    


# p(0) from operdetermination:
# phi'(0) = u_x(L,0) - u_x(0,0) - p(0) * (u(L,0) - u(0,0))
# phi'(0) = u_xx(x_bar, 0) + p(0) * u0(x_bar)
x_bar = 0.2
n = int(round(x_bar / h))
phi = u[n, :]

u0_xx = (u[0, n + 1] - 2.0 * u[0,n] + u[0,n - 1]) / h**2
phi_t0 = (-3.0 * phi[0] + 4.0 * phi[1] - phi[2]) / (2.0 * tau)
p0_est = (phi_t0 - u0_xx) / u[0,n]

y=np.zeros((J+1, N+1)); z=np.zeros(N+1); v=np.zeros(N+1); p=np.zeros(J+1)
y[0,:]=u[0,:]
p[0]=p0_est
for j in range(1,J):
    A=-1.; B=2+r-h**2*p[j-1]; C=-1.; z[N]=mu2; v[N]=0
    ksi_z=np.zeros(N+1); eta_z=np.zeros(N+1); ksi_z[0]=mu1; eta_z[0]=0
    ksi_v=np.zeros(N+1); eta_v=np.zeros(N+1); ksi_v[0]=0;   eta_v[0]=0
    for i in range(1, N-1):
        ksi_z[i+1]= -C / (A*ksi_z[i] + B)
        eta_z[i+1]= (r*u[j-1][i]+(y[j-1][i+1]-2.*y[j-1][i]+y[j-1][i-1]) - A*eta_z[i]) / (A*ksi_z[i] + B)
        ksi_v[i+1]= -C / (A*ksi_v[i] + B)
        eta_v[i+1]= ( h**2*y[j-1][i] - A*eta_v[i]) / (A*ksi_v[i] + B)
    
    for i in range(N - 1, 0, -1): 
        z[i]=ksi_z[i+1]*z[i+1]+eta_z[i+1]; v[i]=ksi_v[i+1]*v[i+1]+eta_v[i+1]
    p[j] = (phi[j] - z[n]) / v[n]
    for i in range(N+1): y[j,i] = z[i] + p[j]*v[i] 


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
plt.semilogy(abs(p-p_e))
plt.title("Погрешность восстановления p-p_e")

plt.show()