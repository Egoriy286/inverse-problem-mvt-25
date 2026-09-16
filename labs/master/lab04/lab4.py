import numpy as np, matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parents[3] / "results" / "lab04"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

l,T,n,J=1.,0.8,100,100
h,tau=l/n,T/J
x=np.linspace(0,l,n+1); t=np.linspace(0,T,J+1)

v_e=0.5+0.2*np.sin(2*np.pi*t)
u0=lambda x: np.exp(-40*(x-0.5)**2)
mu=lambda t: 0.

up=lambda y,v,t: np.r_[mu(t), y[1:]-v*tau/h*(y[1:]-y[:-1])]
psi=lambda y: h/2*(y[:-1].sum()+y[1:].sum())

# прямая
y=u0(x); y[0]=0
U=np.zeros((J+1,n+1)); U[0]=y
Psi=np.zeros(J+1); Psi[0]=psi(y)
for j in range(1,J+1): y=up(y,v_e[j],t[j]); U[j]=y; Psi[j]=psi(y)

# обратная (бисекция)
y=u0(x); y[0]=0
Ur=np.zeros_like(U); Ur[0]=y
v=np.zeros(J+1); vmax=0.95*h/tau
v[0] =v_e[0]
for j in range(1,J+1):
    tj,Ps=t[j],Psi[j]
    F=lambda vv: psi(up(y,vv,tj))-Ps
    a,b=0.,vmax; Fa,Fb=F(a),F(b)
    if Fa*Fb>0: vj=a if abs(Fa)<abs(Fb) else b
    else:
        for _ in range(25):
            m=(a+b)/2; Fm=F(m)
            if Fa*Fm<=0: b,Fb=m,Fm
            else: a,Fa=m,Fm
        vj=(a+b)/2
    v[j]=vj; y=up(y,vj,tj); Ur[j]=y

# графики
fig,ax=plt.subplots(2,2,figsize=(10,8))
c=ax[0,0].contourf(x,t,U); fig.colorbar(c,ax=ax[0,0]); ax[0,0].set_title("точное решение"); 
ax[0,0].set_ylabel('t'); ax[0,0].set_xlabel('x'); 
c=ax[0,1].contourf(x,t,Ur); fig.colorbar(c,ax=ax[0,1]); ax[0,1].set_title("приближенное решение"); 
ax[0,1].set_ylabel('t'); ax[0,1].set_xlabel('x')
ax[1,0].plot(t,v_e,label='точное решение'); ax[1,0].plot(t,v,label='приближенное решение')
ax[1,0].legend(); ax[1,0].grid(); ax[1,0].set_title("v(t)")
ax[1,1].plot(t,abs(v-v_e)); ax[1,1].grid(); ax[1,1].set_title("error")
plt.tight_layout()

# анимация (быстрая)
fig,ax=plt.subplots(); line,=ax.plot([],[])
ax.set(xlim=(0,l),ylim=(1.1*U.min(),1.1*U.max()))
txt=ax.text(0.02,0.9,'',transform=ax.transAxes)
skip=max(1,J//50)

FuncAnimation(fig,lambda f:(line.set_data(x,U[f*skip]) or txt.set_text(f"t={t[f*skip]:.2f}") or (line,txt)),
              frames=J//skip,blit=True)\
.save(OUTPUT_DIR / "laba4-generated.gif",writer=PillowWriter(fps=20))

plt.show()
