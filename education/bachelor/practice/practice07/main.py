import numpy as np
import matplotlib.pyplot as plt

def maps(value, titles):
    fig, ax1 = plt.subplots(figsize=(6, 4))
    im1 = ax1.imshow(value, extent=[ 0, Len, 0, T], aspect='auto', cmap='hot', origin='lower', interpolation='nearest')
    ax1.set_ylabel("Время (t)"); ax1.set_xlabel("Координата (x)"); ax1.set_title(f"{titles}(x,t) "); ax1.grid(True)
    fig.colorbar(im1, ax=ax1, orientation='vertical')
    



Len = 1.0
N = 100
T = 1.0
M = 100
h = Len / N
tau = T / M
x = np.linspace(0,Len,N+1)
t = np.linspace(0,T,  M+1)

def psi(t):
    return (2 * t / T) if (t <= T/2) else (2 * (T - t) / T) 

A = -( 1 / h**2 )
B =  ( 2 / h**2 + 1 / tau)
C = -( 1 / h**2)

u = np.zeros((M+1,N+1))
u[0][:] = 0 

psi = np.array([psi(i) for i in t])

for j in range(M):
    alpha = np.zeros((N+1))
    beta = np.zeros((N+1))
    alpha[0] = 1
    beta[0] = 0
    # Прямой ход прогонки
    for i in range(0, N):
        alpha[i+1] = -C / (A * alpha[i] + B)
        beta[i+1] = (u[j,i]/tau - A * beta[i]) / (A * alpha[i] + B)
    u[j+1,N] = psi[j+1]
    for i in range(N-1, -1, -1):
        u[j+1][i] = alpha[i+1]*u[j+1][i+1] + beta[i+1]

phi = np.zeros(N+1)
for j in range(M+1):
    phi[j] = u[j][0]
maps(u, "u")
delta = 0.005
sigma = np.random.normal(0,1,size=(N+1))
phi_delta = phi + delta * sigma
plt.figure(figsize=(6,4))
plt.plot(t, phi_delta)
plt.title(f"Возмущенный $\phi_\delta(t)$"); plt.ylabel("$u(0,t)$"); plt.xlabel("t"); plt.grid()


q = np.zeros(N+1)
qN = 1.0

alpha = np.zeros((N+1))
beta = np.zeros((N+1))
alpha[0] = 1
beta[0] = 0
# Прямой ход прогонки
for i in range(0, N):
    alpha[i+1] = -C / (A * alpha[i] + B)
    beta[i+1] = (- A * beta[i]) / (A * alpha[i] + B)
q[N] = qN
for i in range(N-1, -1, -1):
    q[i] = alpha[i+1]*q[i+1] + beta[i+1]


u_alpha = np.zeros((M+1, N+1))
alpha_const = 0.05
z = np.zeros((M+1, N+1))
v = np.zeros((M+1))
for j in range(M):
    alpha = np.zeros((N+1))
    beta = np.zeros((N+1))
    alpha[0] = 1
    beta[0] = 0
    # Прямой ход прогонки
    for i in range(0, N):
        alpha[i+1] = -C / (A * alpha[i] + B)
        beta[i+1] = ((u_alpha[j][i] - z[j+1][i]) / tau - A * beta[i]) / (A * alpha[i] + B)
    z[j+1,N] = 0.0
    for i in range(N-1, -1, -1):
        z[j+1][i] = alpha[i+1]*z[j+1][i+1] + beta[i+1]

    v[j+1] = ((phi_delta[j+1] - z[j+1][0]) / (alpha_const + q[0]))
    
    u_alpha[j+1] = z[j+1] + q*v[j+1]


plt.figure(figsize=(6,4))
plt.plot(t, v, label="$v(t)$")
plt.plot(t, psi,'--', label="$\psi(t)$")
plt.title(f"$u(l,t)$"); plt.ylabel("$u(0,t)$"); plt.xlabel("t"); plt.grid(); plt.legend()


plt.figure(figsize=(6,4))
plt.plot(t, phi, label="$v(t)$")
plt.plot(t, phi_delta, label="$\psi(t)$")
plt.title(f"$u(0,t)$"); plt.ylabel("$u(0,t)$"); plt.xlabel("t"); plt.grid(); plt.legend()


plt.figure(figsize=(6,4))
plt.plot(x, q, label="$q(x)$")
plt.title(f"$q(x)$"); plt.ylabel("$q(x)$"); plt.xlabel("x"); plt.grid(); plt.legend()
plt.show()
