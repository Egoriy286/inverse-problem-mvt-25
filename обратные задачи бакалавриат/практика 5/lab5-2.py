import numpy as np
import matplotlib.pyplot as plt
from tqdm.notebook import tqdm

# Plotting function
def maps(value, titles):
    fig, ax1 = plt.subplots(figsize=(6, 4))
    im1 = ax1.imshow(value, extent=[0, L, 0, T], aspect='auto', cmap='hot', origin='lower', interpolation='nearest')
    ax1.set_ylabel("Время (t)")
    ax1.set_xlabel("Координата (x)")
    ax1.set_title(f"Решение {titles}(x,t)")
    ax1.grid(True)
    fig.colorbar(im1, ax=ax1, orientation='vertical')
    plt.show()

# L2 norm over space-time grid
def Norm(x):
    return np.sqrt(np.sum(x**2))

# Exact source term
def f(x, t):
    return 10 * t * (1 - t) * x * (1 - x)

# Thomas algorithm for tridiagonal system: solves Au[i-1] + Bu[i] + Cu[i+1] = rhs[i] for i=1 to N-1
def solve_tridiagonal(A, B, C, rhs, N):
    alpha = np.zeros(N-1)  # For i=1 to N-1
    beta = np.zeros(N-1)
    u = np.zeros(N+1)      # Includes boundaries

    # Forward elimination
    denom = B
    alpha[0] = -C / denom
    beta[0] = rhs[1] / denom  # rhs[1] corresponds to i=1
    for k in range(1, N-1):
        denom = A * alpha[k-1] + B
        alpha[k] = -C / denom
        beta[k] = (rhs[k+1] - A * beta[k-1]) / denom

    # Backward substitution
    u[N-1] = beta[N-2]  # i=N-1
    for k in range(N-3, -1, -1):
        u[k+1] = alpha[k] * u[k+2] + beta[k]
    # Boundaries u[0] and u[N] remain 0
    return u

# Parameters
L = 1.0
T = 1.0
N = 100
M = 100
h = L / N
tau = T / M
x = np.linspace(0, L, N+1)
t = np.linspace(0, T, M+1)

# Tridiagonal coefficients
A = -1 / h**2
B = 2 / h**2 + 1 / tau
C = -1 / h**2

# Compute exact source
f_values = np.zeros((M+1, N+1))
for j in range(M+1):  # Include all time steps for consistency
    for i in range(N+1):
        f_values[j, i] = f(x[i], t[j])

# Solve direct problem: u_t = u_xx + f
u = np.zeros((M+1, N+1))
for j in range(M):
    rhs = u[j] / tau + f_values[j]
    u[j+1] = solve_tridiagonal(A, B, C, rhs, N)
maps(u, "u")

# Add noise to simulate measured data
delta = 0.001
sigma = np.random.normal(0, 1, size=(M+1, N+1))
u_delta = u #+ delta * sigma

# Solve adjoint problem: y_t = -y_xx + u_delta (backward in time)
y_delta = np.zeros((M+1, N+1))  # y_delta[M] = 0 as final condition
for j in reversed(range(M)):
    rhs = -y_delta[j+1] / tau + u_delta[j+1]
    y_delta[j] = solve_tridiagonal(A, B, C, rhs, N)
maps(y_delta, "y_delta")

# Iterative restoration of f
fk = np.zeros((M+1, N+1))  # Initial guess
v = np.zeros((M+1, N+1))
w = np.zeros((M+1, N+1))
rk = np.zeros((M+1, N+1))
Grk = np.zeros((M+1, N+1))

pbar = tqdm(range(2000))
for k in pbar:
    # Direct problem: v = G f_k
    v.fill(0)  # Reset v
    for j in range(M):
        rhs = v[j] / tau + fk[j]
        v[j+1] = solve_tridiagonal(A, B, C, rhs, N)

    # Adjoint problem: w = G^* v
    w.fill(0)
    for j in reversed(range(M)):
        rhs = -w[j+1] / tau + v[j+1]
        w[j] = solve_tridiagonal(A, B, C, rhs, N)

    # Residual
    rk = w - y_delta

    # Direct problem: Grk = G rk
    Grk.fill(0)
    for j in reversed(range(M)):
        rhs = -Grk[j+1] / tau + rk[j+1]
        Grk[j] = solve_tridiagonal(A, B, C, rhs, N)

    # Iteration step
    tau_k = Norm(rk) / Norm(Grk)
    fk_new = fk - tau_k * rk

    # Stopping criterion
    error = Norm(v - u_delta)
    pbar.set_description(f"tau: {tau_k:.6f}, error: {error:.6f}")
    if error <= delta * np.sqrt(N):
        fk = fk_new
        break

    fk = fk_new

maps(fk, "Восстановленная правая часть f")
maps(f_values, "Точная правая часть f_exact")

# Plot results at t = 0.5
plt.figure(figsize=(10, 6))
plt.plot(x, fk[M//2], label='Восстановленный f')
plt.plot(x, f_values[M//2], '--', label='Точный f')
plt.xlabel('x')
plt.ylabel('f(x, t)')
plt.title('Восстановленный источник при t = 0.5')
plt.legend()
plt.grid(True)
plt.show()

# Compute and print restoration error
error = np.linalg.norm(fk[M//2] - f_values[M//2])
print(f"Ошибка восстановления при t = 0.5: {error:.6f}")