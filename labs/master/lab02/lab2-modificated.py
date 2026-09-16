import numpy as np; import matplotlib.pyplot as plt;
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parents[3] / "results" / "lab02"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
L=2*np.pi; T=1.0; J=100; N=100; mu1=0; mu2=0
tau=T/J; h=L/N; r=2*h**2/tau; c=2+r; e=2-r
u=np.zeros((J+1, N+1)); x=np.linspace(0, L, N+1); t=np.linspace(0, T, J+1)
u[0,:]=np.sin(2*np.pi * x/L); v_e = np.sin(2*np.pi*t)*10; phi = np.zeros(J+1)
for j in range(J):
    A=-1.; B=r+2.; C=-1.
    ksi=np.zeros(N+1); eta=np.zeros(N+1); ksi[1]=0; eta[1]=mu1
    for i in range(1, N):
        ksi[i+1]= -C / (A*ksi[i] + B)
        eta[i+1]= (r*u[j,i]+u[j,i+1]-2*u[j,i]+u[j,i-1]-v_e[j+1]*2*h*(u[j][i] - u[j,i-1]) - A*eta[i]) / (A*ksi[i] + B)
    u[j+1][N]=mu2
    for i in range(N - 1, -1, -1): u[j+1][i]=ksi[i+1]*u[j+1][i+1]+eta[i+1]
    u[j+1,0]=mu1
    
for j in range(J+1): 
    phi[j] = h*(sum(u[j]) - 0.5*(u[j,0] + u[j,N])) #sum(h/2 * (u[j, i] + u[j, i-1]) for i in range(1, N+1))
    
# v(0) from overdetermination:
# phi'(0) = u_x(L,0) - u_x(0,0) - v(0) * (u(L,0) - u(0,0))
phi_t0 = (-3.0 * phi[0] + 4.0 * phi[1] - phi[2]) / (2.0 * tau)
u0_ = u[0, :]
u0_x0 = (-3.0 * u0_[0] + 4.0 * u0_[1] - u0_[2]) / (2.0 * h)
u0_xL = (3.0 * u0_[N] - 4.0 * u0_[N - 1] + u0_[N - 2]) / (2.0 * h)
den_v0 = u0_[N] - u0_[0]
if abs(den_v0) > 1e-12:
    v0_est = (u0_xL - u0_x0 - phi_t0) / den_v0
else:
    v0_est = 0

y=np.zeros((J+1, N+1)); v=np.zeros(J+1)
y[0,:]=u[0,:]
v[0]=v0_est
for j in range(J):
    z=np.zeros(N+1); w=np.zeros(N+1)
    A=-1.; B=c; C=-1.;z[N]=mu2; w[N]=0
    ksi_z=np.zeros(N+1); eta_z=np.zeros(N+1); ksi_z[1]=0; eta_z[1]=mu1
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
    
    den = h*(sum(w) - 0.5*(w[0] + w[N]))
    if abs(den) > 1e-12:
        v[j + 1] = (phi[j + 1] - h*(sum(z) - 0.5*(z[0] + z[N]))) / den
    else:
        v[j + 1] = v[j]
        
    y[j+1,:] = z + v[j+1]*w 
    y[j+1,0]=mu1; y[j+1,N]=mu2

# Фильтр Совицкий Галея Статья и шум добавить.
# Визуализация решения
p = plt.contourf(y); plt.colorbar(p); plt.title("Приближенное решение")
plt.figure(); p = plt.contourf(u); plt.colorbar(p); plt.title("Точное решение")
plt.figure(); plt.plot(v_e, label='exact'); plt.plot(v, label='aprox'); plt.legend(); plt.title("v(t)="); plt.grid()
plt.figure(); plt.plot(abs(v-v_e)); plt.title("Погрешность восстановления v-v_e")


from matplotlib.animation import FuncAnimation, PillowWriter

fig, ax = plt.subplots(); line, = ax.plot([], [])
ax.set(xlim=(0, L), ylim=(1.1*np.min(u), 1.1*np.max(u)))
txt = ax.text(0.02, 0.9, '', transform=ax.transAxes)

skip = max(1, J//50)

FuncAnimation(fig, lambda f: (line.set_data(x, u[f*skip]) or txt.set_text(f"t={t[f*skip]:.2f}") or (line, txt)),
              frames=J//skip, blit=True)\
.save(OUTPUT_DIR / "laba2-generated.gif", writer=PillowWriter(fps=20))

plt.show()
