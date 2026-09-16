import numpy as np; import matplotlib.pyplot as plt

T = 0.8
L = 2.0
N = 100
M = 100

def v(t):
    return np.sin(np.pi*t/T)
def mu(t):
    return np.cos(np.pi * t)
def u0(x):
    return np.sin(2*np.pi*x/L)

tau = T/M
h = L/N
print("tau:", tau, "h:", h, "\nУсловие устойчивости:", tau<=h)

u = np.zeros((M+1,N+1))
x = np.linspace(0, L, N+1)
t = np.linspace(0, T, M+1)
u[0,:] = u0(x)
psi = np.zeros(M+1)
psi[0] = h/2 * (u[0,:-1].sum() + u[0,1:].sum())

for j in range(1, M+1):           # ФИКС 1: range(1,M) → range(1,M+1), иначе последний шаг пропускается
    u[j, 0] = mu(t[j])            # ФИКС 2: u[j-1,0]=0 → u[j,0]=mu(t[j]), граничное условие на текущий слой
    for i in range(1, N+1):
        u[j,i] = u[j-1,i] - tau/h * v(t[j-1]) * (u[j-1,i]-u[j-1,i-1])
    s = 0
    for i in range(1, N+1):
        s += (u[j, i-1] + u[j, i])
    psi[j] = h/2 * s

plt.figure()
p = plt.contourf(u)
plt.colorbar(p)
plt.xlabel("x")
plt.ylabel("t")

y = np.zeros((M+1,N+1))
y[0,:] = u0(x)
# ФИКС 3: убрана строка y[:,0]=0 — она обнуляла граничное условие для всех t сразу
v_r = np.zeros(M+1)

for j in range(1, M+1):           # ФИКС 4: range(1,M) → range(1,M+1)
    yn   = y[j-1, N]
    yn1  = y[j-1, N-1]
    denom = yn - yn1
    if abs(denom) > 1e-10:        # ФИКС 5: формула (8) с доски вместо psi[j-1]-y[j-1,N]
        v_ = (yn - psi[j]) * h / tau / denom
    else:
        v_ = 0.0
    v_r[j] = v_
    y[j, 0] = mu(t[j])            # граничное условие
    for i in range(1, N+1):
        y[j,i] = y[j-1,i] - v_ * tau/h * (y[j-1,i]-y[j-1,i-1])  # ФИКС 6: была просто копия y[j-1,i]

plt.figure()
p = plt.contourf(y)
plt.colorbar(p)
plt.xlabel("x")
plt.ylabel("t")

plt.figure()
plt.plot(v(t), label="v(t) истинное")
plt.plot(v_r,  label="v(t) восстановленное")
plt.legend()

plt.show()                         # ФИКС 7: убрана сломанная markdown-ссылка