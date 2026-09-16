import numpy as np
import matplotlib.pyplot as plt


# PDE:
# u_t = u_xx - p(t) * u,   0 < x < L, 0 < t <= T
# u(0,t) = mu1(t), u(L,t) = mu2(t), u(x,0) = u0(x)
# extra data for inverse problem: u(x_obs, t) = phi(t)

L = 1.0
T = 1.0
M = 100   # space steps
N = 200   # time steps

h = L / M
tau = T / N

x = np.linspace(0.0, L, M + 1)
t = np.linspace(0.0, T, N + 1)
r = 2.0 * h * h / tau

x_obs = 0.5
n_obs = int(round(x_obs / h))
x_obs = x[n_obs]


def mu1(t_val):
    return 0.0


def mu2(t_val):
    return 0.0


def p_true(t_val):
    return np.exp(t_val)


def u0(x_val):
    return np.sin(np.pi * x_val)


def thomas(lower, diag, upper, rhs):
    n = len(rhs)
    c_mod = np.zeros(n)
    d_mod = np.zeros(n)

    c_mod[0] = upper[0] / diag[0]
    d_mod[0] = rhs[0] / diag[0]

    for i in range(1, n):
        denom = diag[i] - lower[i] * c_mod[i - 1]
        if i < n - 1:
            c_mod[i] = upper[i] / denom
        d_mod[i] = (rhs[i] - lower[i] * d_mod[i - 1]) / denom

    sol = np.zeros(n)
    sol[-1] = d_mod[-1]
    for i in range(n - 2, -1, -1):
        sol[i] = d_mod[i] - c_mod[i] * sol[i + 1]

    return sol


def build_z_v(u_prev, t_next):
    size = M - 1
    lower = -np.ones(size)
    diag = (2.0 + r) * np.ones(size)
    upper = -np.ones(size)

    rhs_z = np.zeros(size)
    rhs_v = np.zeros(size)

    for i in range(1, M):
        idx = i - 1
        lap = u_prev[i + 1] - 2.0 * u_prev[i] + u_prev[i - 1]
        rhs_z[idx] = r * u_prev[i] + lap
        rhs_v[idx] = -h * h * u_prev[i]

    # non-homogeneous boundaries are included only in z
    rhs_z[0] += mu1(t_next)
    rhs_z[-1] += mu2(t_next)

    z_inner = thomas(lower, diag, upper, rhs_z)
    v_inner = thomas(lower, diag, upper, rhs_v)

    z = np.zeros(M + 1)
    v = np.zeros(M + 1)

    z[0] = mu1(t_next)
    z[M] = mu2(t_next)
    z[1:M] = z_inner

    v[1:M] = v_inner

    return z, v


def solve_direct():
    u = np.zeros((N + 1, M + 1))
    u[0, :] = u0(x)

    for j in range(N):
        z, v = build_z_v(u[j, :], t[j + 1])
        p_next = p_true(t[j + 1])
        u[j + 1, :] = z + p_next * v
        u[j + 1, 0] = mu1(t[j + 1])
        u[j + 1, M] = mu2(t[j + 1])

    return u


def solve_inverse(phi_data):
    u_inv = np.zeros((N + 1, M + 1))
    p_rec = np.zeros(N + 1)
    p_rec[0] = p_true(t[0])
    u_inv[0, :] = u0(x)

    for j in range(N):
        z, v = build_z_v(u_inv[j, :], t[j + 1])

        denom = v[n_obs]
        if abs(denom) < 1e-14:
            p_rec[j + 1] = p_rec[j]
        else:
            p_rec[j + 1] = (phi_data[j + 1] - z[n_obs]) / denom

        u_inv[j + 1, :] = z + p_rec[j + 1] * v
        u_inv[j + 1, 0] = mu1(t[j + 1])
        u_inv[j + 1, M] = mu2(t[j + 1])

    return u_inv, p_rec


# 1) direct problem
u_direct = solve_direct()

# 2) synthetic extra data for inverse problem from direct solution
phi_data = u_direct[:, n_obs].copy()

# 3) inverse problem
u_inv, p_rec = solve_inverse(phi_data)

p_vals = p_true(t)
p_err = np.max(np.abs(p_vals - p_rec))
print(f"x_obs = {x_obs:.3f}")
print(f"max |p_true - p_rec| = {p_err:.3e}")

plt.figure(figsize=(8, 5))
plt.plot(t, p_vals, label="p true")
plt.plot(t, p_rec, "--", label="p recovered")
plt.title("Inverse reconstruction of p(t)")
plt.xlabel("t")
plt.ylabel("p(t)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(x, u_direct[-1, :], label="u_direct(x, T)")
plt.plot(x, u_inv[-1, :], "--", label="u_inverse(x, T)")
plt.title("Direct vs inverse solution at final time")
plt.xlabel("x")
plt.ylabel("u(x, T)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
