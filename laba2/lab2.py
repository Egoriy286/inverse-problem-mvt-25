import numpy as np; import matplotlib.pyplot as plt;
L=10; T=0.1; J=100; N=100; mu1=0; mu2=0
tau=T/J; h=L/N; r=2*h**2/tau; c=2+r; e=2-r
u=np.zeros((J+1, N+1)); x=np.linspace(0, L, N+1); t=np.linspace(0, T, J+1)
u[0,:]=np.sin(np.pi * x); v_e = np.sin(t); phi = np.zeros(J+1)
for j in range(J):
    A=-1.; B=r+2.; C=-1.
    ksi=np.zeros(N+1); eta=np.zeros(N+1); ksi[1]=mu1; eta[1]=0
    for i in range(1, N):
        ksi[i+1]= -C / (A*ksi[i] + B)
        eta[i+1]= (r*u[j,i]+u[j,i+1]-2*u[j,i]+u[j,i-1]-v_e[j]*2*h*(u[j][i] - u[j,i-1]) - A*eta[i]) / (A*ksi[i] + B)
        u[j][N]=mu2
    for i in range(N - 1, -1, -1): u[j+1][i]=ksi[i+1]*u[j+1][i+1]+eta[i+1]
    
for j in range(J+1): 
    phi[j] = h*sum(u[j])/2#sum(h/2 * (u[j, i] + u[j, i-1]) for i in range(1, N+1))
    
# v(0) from overdetermination:
# phi'(0) = u_x(L,0) - u_x(0,0) - v(0) * (u(L,0) - u(0,0))
phi_t0 = (-3.0 * phi[0] + 4.0 * phi[1] - phi[2]) / (2.0 * tau)
u0_ = u[0, :]
u0_x0 = (-3.0 * u0_[0] + 4.0 * u0_[1] - u0_[2]) / (2.0 * h)
u0_xL = (3.0 * u0_[N] - 4.0 * u0_[N - 1] + u0_[N - 2]) / (2.0 * h)
den_v0 = u0_[N] - u0_[0]
v0_est = (u0_xL - u0_x0 - phi_t0) / den_v0

y=np.zeros((J+1, N+1)); v=np.zeros(J+1)
y[0,:]=u[0,:]
v[0]=v0_est
for j in range(J):
    z=np.zeros(N+1); w=np.zeros(N+1)
    A=-1.; B=c; C=-1.;z[N]=mu2; w[N]=mu2
    ksi_z=np.zeros(N+1); eta_z=np.zeros(N+1); ksi_z[1]=mu1; eta_z[1]=0
    ksi_w=np.zeros(N+1); eta_w=np.zeros(N+1); ksi_w[1]=0;   eta_w[1]=0
    for i in range(1, N):
        ksi_z[i+1]= -C / (A*ksi_z[i] + B)
        eta_z[i+1]= (y[j,i+1]-e*y[j,i]+y[j,i-1] - A*eta_z[i]) / (A*ksi_z[i] + B)
        ksi_w[i+1]= -C / (A*ksi_w[i] + B)
        eta_w[i+1]= (-2*h*(y[j,i] - y[j,i-1]) - A*eta_w[i]) / (A*ksi_w[i] + B)
    
    for i in range(N - 1, -1, -1): 
        z[i]=ksi_z[i+1]*z[i+1]+eta_z[i+1]
        w[i]=ksi_w[i+1]*w[i+1]+eta_w[i+1]
    # v[j] = ((phi[j] - sum(h/2 * (z[i] + z[i-1]) for i in range(1, N+1)))
    #         / (sum(h/2 * (w[i] + w[i-1]) for i in range(1, N+1))))
    # v[j + 1] = (phi[j + 1] - h*sum(z)/2) / (h*sum(w)/2)
    
    den = h*sum(w)
    if abs(den) > 1e-12:
        v[j + 1] = (phi[j + 1] - h*sum(z)) / den
    else:
        v[j + 1] = v[j]
        
    y[j+1,:] = z + v[j+1]*w 
    y[j+1,0]=mu1; y[j+1,N]=mu2


# Визуализация решения
p = plt.imshow(y)
plt.colorbar(p)
plt.title("Приближенное решение")

plt.figure()
p = plt.imshow(u)
plt.colorbar(p)
plt.title("Точное решение")

plt.figure()
plt.plot(v_e, label='exact')
plt.plot(v, label='aprox')
plt.legend()
plt.title("v(t)=")
plt.grid()

plt.figure()
plt.semilogy(abs(v-v_e))
plt.title("Погрешность восстановления v-v_e")

plt.show()