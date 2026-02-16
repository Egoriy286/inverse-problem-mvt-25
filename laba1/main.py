import numpy as np
import matplotlib.pyplot as plt

# ==========================================================
# ПАРАМЕТРЫ
# ==========================================================
l = 1.0
t_max = 1.0
N = 60
J = 200
x_bar = 0.5

h = l / N
tau = t_max / J
r = 2.0 * h**2 / tau

x = np.linspace(0.0, l, N + 1)
t = np.linspace(0.0, t_max, J + 1)

n = int(round(x_bar / h))
if n <= 0 or n >= N:
    raise ValueError("x_bar должен быть внутренним узлом: 0 < n < N")

eps = 1e-12

print(f"h = {h:.6f}, tau = {tau:.6f}, r = {r:.6f}")
print(f"x_bar = {x[n]:.6f}, индекс n = {n}")


# ==========================================================
# ИСХОДНЫЕ ДАННЫЕ
# ==========================================================
def u0(x_val):
    return np.sin(np.pi * x_val)


def mu1(t_val):
    return 0.0


def mu2(t_val):
    return 0.0


def p_true(t_val):
    return np.sin(t_val) #1.5 + 0.4 * np.sin(2.0 * np.pi * t_val)


def thomas(lower, diag, upper, rhs):
    """Решение трехдиагональной системы."""
    n_sys = rhs.size
    c_star = np.zeros(n_sys)
    d_star = np.zeros(n_sys)

    c_star[0] = upper[0] / diag[0] if n_sys > 1 else 0.0
    d_star[0] = rhs[0] / diag[0]

    for i in range(1, n_sys):
        denom = diag[i] - lower[i] * c_star[i - 1]
        c_star[i] = upper[i] / denom if i < n_sys - 1 else 0.0
        d_star[i] = (rhs[i] - lower[i] * d_star[i - 1]) / denom

    sol = np.zeros(n_sys)
    sol[-1] = d_star[-1]
    for i in range(n_sys - 2, -1, -1):
        sol[i] = d_star[i] - c_star[i] * sol[i + 1]
    return sol


# ==========================================================
# 1) ПРЯМАЯ ЗАДАЧА: ГЕНЕРАЦИЯ phi(t)
# ==========================================================
p_ref = p_true(t)
y_ref = np.zeros((N + 1, J + 1))
y_ref[:, 0] = u0(x)
y_ref[0, 0] = mu1(t[0])
y_ref[N, 0] = mu2(t[0])

for j in range(1, J + 1):
    y_prev = y_ref[:, j - 1]
    p_prev = p_ref[j - 1]
    p_cur = p_ref[j]

    size = N - 1
    lower = np.ones(size)
    diag = np.full(size, -(2.0 + r - h**2 * p_prev))
    upper = np.ones(size)
    lower[0] = 0.0
    upper[-1] = 0.0

    rhs = np.zeros(size)
    for i in range(1, N):
        idx = i - 1
        d_bar = r * y_prev[i] + (y_prev[i + 1] - 2.0 * y_prev[i] + y_prev[i - 1])
        rhs[idx] = -(d_bar + h**2 * p_cur * y_prev[i])

    rhs[0] -= mu1(t[j])
    rhs[-1] -= mu2(t[j])

    y_inner = thomas(lower, diag, upper, rhs)

    y_ref[0, j] = mu1(t[j])
    y_ref[N, j] = mu2(t[j])
    y_ref[1:N, j] = y_inner

phi = y_ref[n, :].copy()

# p(0) из условия переопределения:
# phi'(0) = u_xx(x_bar, 0) + p(0) * u0(x_bar)
u0_grid = u0(x)
u0_xx_bar = (u0_grid[n + 1] - 2.0 * u0_grid[n] + u0_grid[n - 1]) / h**2
if J >= 2:
    phi_t0 = (-3.0 * phi[0] + 4.0 * phi[1] - phi[2]) / (2.0 * tau)
else:
    phi_t0 = (phi[1] - phi[0]) / tau

if abs(u0_grid[n]) < eps:
    raise ValueError("Невозможно определить p(0): u0(x_bar) близко к нулю")

p0_est = (phi_t0 - u0_xx_bar) / u0_grid[n]
print(f"estimated p(0) from phi(t): {p0_est:.6f}")


# ==========================================================
# 2) ОБРАТНАЯ ЗАДАЧА: ВОССТАНОВЛЕНИЕ p(t), u(x,t)
# ==========================================================
y = np.zeros((N + 1, J + 1))
p = np.zeros(J + 1)

y[:, 0] = u0_grid
y[0, 0] = mu1(t[0])
y[N, 0] = mu2(t[0])
p[0] = p0_est
for j in range(1, J + 1):
    y_prev = y[:, j - 1]
    p_prev = p[j - 1]

    size = N - 1
    lower = np.ones(size)
    diag = np.full(size, -(2.0 + r - h**2 * p_prev))
    upper = np.ones(size)
    lower[0] = 0.0
    upper[-1] = 0.0

    # Система для z^j
    rhs_z = np.zeros(size)
    for i in range(1, N):
        idx = i - 1
        d_bar = r * y_prev[i] + (y_prev[i + 1] - 2.0 * y_prev[i] + y_prev[i - 1])
        rhs_z[idx] = -d_bar

    rhs_z[0] -= mu1(t[j])
    rhs_z[-1] -= mu2(t[j])

    z = np.zeros(N + 1)
    z[0] = mu1(t[j])
    z[N] = mu2(t[j])
    z[1:N] = thomas(lower, diag, upper, rhs_z)

    # Система для v^j
    rhs_v = -h**2 * y_prev[1:N]
    v = np.zeros(N + 1)
    v[1:N] = thomas(lower, diag, upper, rhs_v)

    # Восстановление p^j
    if abs(v[n]) < eps:
        p[j] = p_prev
    else:
        p[j] = (phi[j] - z[n]) / v[n]

    # Восстановление y^j
    y[:, j] = z + p[j] * v
    y[0, j] = mu1(t[j])
    y[N, j] = mu2(t[j])


# ==========================================================
# КОНТРОЛЬ ТОЧНОСТИ
# ==========================================================
misfit_phi = np.max(np.abs(y[n, :] - phi))
err_p_max = np.max(np.abs(p - p_ref))
err_u_max = np.max(np.abs(y - y_ref))

print(f"max |y(x_bar,t) - phi(t)| = {misfit_phi:.3e}")
print(f"max |p_rec(t) - p_true(t)| = {err_p_max:.3e}")
print(f"max |u_rec(x,t) - u_ref(x,t)| = {err_u_max:.3e}")


# ==========================================================
# ГРАФИКИ
# ==========================================================
fig1 = plt.figure()
plt.plot(t, p_ref, "k--", lw=2, label="p_true(t)")
plt.plot(t, p, "r", lw=1.8, label="p_rec(t)")
plt.title("Восстановление коэффициента p(t)")
plt.xlabel("t")
plt.ylabel("p(t)")
plt.grid(True, alpha=0.3)
plt.legend()

fig2 = plt.figure()
plt.plot(x, y_ref[:, -1], "k--", lw=2, label="u_ref(x,T)")
plt.plot(x, y[:, -1], "b", lw=1.8, label="u_rec(x,T)")
plt.title("Восстановление решения u(x,T)")
plt.xlabel("x")
plt.ylabel("u(x,T)")
plt.grid(True, alpha=0.3)
plt.legend()

fig3 = plt.figure()
plt.semilogy(t, np.abs(p_ref - p), "k--", lw=2, label="abs(p_exact - p)")
plt.title("Погрешность")
plt.xlabel("t")
plt.ylabel("error")
plt.grid(True, alpha=0.3)
plt.legend()

if "agg" in plt.get_backend().lower():
    fig1.savefig("p_recovery.png", dpi=150, bbox_inches="tight")
    fig2.savefig("u_recovery.png", dpi=150, bbox_inches="tight")
    fig3.savefig("error.png", dpi=150, bbox_inches="tight")
    print("Графики сохранены: p_recovery.png, u_recovery.png")
else:
    plt.show()
    
