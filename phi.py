import numpy as np
import matplotlib.pyplot as plt

def phi(t):
    return np.exp(-2*t) * np.sin(np.pi/4) 

def p(t):
    return np.exp(-1*t) * np.sin(np.pi/8)

def u0(x):
    return np.sin(np.pi * x)/4 + 0.7




def mu1(t):
    return 0.4

def mu2(t):
    return 0.3



J = 100
T = 1.0
tau = T / (J + 1)
N = 100
L = 1.0
h = L / (N + 1)

t = np.linspace(0, 1, J+1)
x = np.linspace(0, 1, N+1)

A = np.zeros((N+1,N+1))
rhs = np.zeros((N+1))
r = 2*h**2/tau
y0 = u0(x)
for i in range(N+1):
    A[i,i] = -2 + r * p(t[0]) / 2 + r
    if (i<N):
        A[i+1, i] = 1
    if (i>0):
        A[i-1, i] = 1
    if (i < N):
        rhs[i] = r - y0[i-1] + 2 * y0[i] - y0[i+1] - r * p(t[0]) / 2

rhs[0] = mu1(0)
rhs[N] = mu2(0)
print(rhs)
print(A)
y = np.zeros((J+1,N+1))
for j in range(J):
    for i in range(N+1):
        A[i,i] = -2 + r * p(t[j+1]) / 2 + r
        if (i < N):
            rhs[i] = r - y0[i-1] + 2 * y0[i] - y0[i+1] - r * p(t[j+1]) / 2
        rhs[0] = mu1(t[j+1])
        rhs[N] = mu2(t[j+1])
    
    y[j+1][:] = np.linalg.solve(A, rhs)
    