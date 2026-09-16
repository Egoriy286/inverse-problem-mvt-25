import numpy as np; import matplotlib.pyplot as plt
T = 0.8
L = 2.0
N = 100
M = 200

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
for j in range(1, M):
    u[j-1, 0] = 0
    for i in range(1, N+1):
        u[j,i] = u[j-1,i] - tau/h * v(t[j-1]) * (u[j-1,i]-u[j-1,i-1])
    s = 0
    for i in range(1, N):
        s += (u[j, i-1] + u[j, i])
    psi[j] = h/2 * s
    
plt.figure()
p = plt.contourf(u)
plt.colorbar(p)
plt.xlabel("x")
plt.ylabel("t")


y = np.zeros((M+1,N+1))
y[0,:] = u0(x)
y[:, 0] = 0
v_r = np.zeros(M+1)
for j in range(1, M):
    s = 0
    for i in range(1, N):
        s += (y[j-1, i-1] + y[j-1, i])
    v_ = psi[j] - s
    v_r[j] = v_
    for i in range(1, N+1): 
        y[j,i] = y[j-1,i] - tau/h * v_ * (y[j-1,i]-y[j-1,i-1])
        
    


plt.figure()
p = plt.contourf(y)
plt.colorbar(p)
plt.xlabel("x")
plt.ylabel("t")

plt.figure()

plt.plot(v(t))
plt.plot(v_r)
plt.show()