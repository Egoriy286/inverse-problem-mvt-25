import numpy as np; import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# Параметры
L=1; T=0.3; N=100; J=200; mu1=0; mu2=0
h=L/N; tau=T/J; r=h**2/tau; c=2+r; e=2-r
x=np.linspace(0,L,N+1); t=np.linspace(0,T,J+1)

# Точное решение: u(x,t)=sin(pi*x)*exp(-pi^2*t), v_exact(t)=sin(t)
v_e=np.sin(t)
u=np.zeros((J+1,N+1)); u[0]=np.sin(np.pi*x)

# --- Шаг 0: прямая задача с v_exact (для phi) ---
for j in range(J):
    ksi=np.zeros(N+1); eta=np.zeros(N+1); eta[1]=mu1
    for i in range(1,N):
        d = r*u[j,i] + u[j,i+1]-2*u[j,i]+u[j,i-1] - v_e[j]*2*h*(u[j,i]-u[j,i-1])
        ksi[i+1]= 1/(c-ksi[i]); eta[i+1]=(d+eta[i])/(c-ksi[i])
    u[j+1,N]=mu2
    for i in range(N-1,-1,-1): u[j+1,i]=ksi[i+1]*u[j+1,i+1]+eta[i+1]
    u[j+1,0]=mu1

# Интегральное наблюдение phi(t) = int_0^L u(x,t) dx
phi = h*(u.sum(axis=1) - 0.5*(u[:,0]+u[:,N]))

# v[0] из переопределённости: phi'(0) = u_x(L,0)-u_x(0,0) - v(0)*(u(L,0)-u(0,0))
phi_t0=(-3*phi[0]+4*phi[1]-phi[2])/(2*tau)
ux0=(-3*u[0,0]+4*u[0,1]-u[0,2])/(2*h); uxL=(3*u[0,N]-4*u[0,N-1]+u[0,N-2])/(2*h)
den0=u[0,N]-u[0,0]
v0=(uxL-ux0-phi_t0)/den0 if abs(den0)>1e-12 else 0

# --- Обратная задача: восстановить v(t) ---
y=np.zeros((J+1,N+1)); v=np.zeros(J+1)
y[0]=u[0]; v[0]=v_e[0]

for j in range(J):
    # Два правых вектора: z (без v) и w (коэф. при v)
    # y^{j+1} = z + v^{j+1}*w, v^{j+1} из phi[j+1]
    kz=np.zeros(N+1); ez=np.zeros(N+1); ez[1]=mu1
    kw=np.zeros(N+1); ew=np.zeros(N+1)  # ew[1]=0
    for i in range(1,N):
        dz = y[j,i+1]-e*y[j,i]+y[j,i-1]
        dw = -2*h*(y[j,i]-y[j,i-1])
        kz[i+1]=1/(c-kz[i]); ez[i+1]=(dz+ez[i])/(c-kz[i])
        kw[i+1]=1/(c-kw[i]); ew[i+1]=(dw+ew[i])/(c-kw[i])
    z=np.zeros(N+1); w=np.zeros(N+1); z[N]=mu2
    for i in range(N-1,-1,-1):
        z[i]=kz[i+1]*z[i+1]+ez[i+1]; w[i]=kw[i+1]*w[i+1]+ew[i+1]

    # v^{j+1} из phi[j+1] = int(z+v*w)dx
    phi_z=h*(z.sum()-0.5*(z[0]+z[N])); phi_w=h*(w.sum()-0.5*(w[0]+w[N]))
    v[j+1]=(phi[j+1]-phi_z)/phi_w if abs(phi_w)>1e-12 else v[j]

    y[j+1]=z+v[j+1]*w; y[j+1,0]=mu1; y[j+1,N]=mu2

# --- Графики ---
fig,axes=plt.subplots(2,2,figsize=(10,7))
axes[0,0].imshow(y,aspect='auto',origin='lower'); axes[0,0].set_title('y — обратная задача')
axes[0,1].imshow(u,aspect='auto',origin='lower'); axes[0,1].set_title('u — точное решение')
axes[1,0].plot(t,v_e,label='v точн.'); axes[1,0].plot(t,v,'--',label='v восст.')
axes[1,0].legend(); axes[1,0].set_title('v(t)'); axes[1,0].grid()
axes[1,1].semilogy(t,np.abs(v-v_e)+1e-16)
axes[1,1].set_title('Погрешность |v-v_e|'); axes[1,1].grid()
plt.tight_layout(); plt.show()



fig,ax=plt.subplots(); line,=ax.plot([],[])
ax.set(xlim=(0,L),ylim=(1.1*u.min(),1.1*u.max()))
txt=ax.text(0.02,0.9,'',transform=ax.transAxes)
skip=max(1,J//50)

FuncAnimation(fig,lambda f:(line.set_data(x,u[f*skip]) or txt.set_text(f"t={t[f*skip]:.2f}") or (line,txt)),
              frames=J//skip,blit=True)\
.save("laba5.gif",writer=PillowWriter(fps=20))

plt.show()