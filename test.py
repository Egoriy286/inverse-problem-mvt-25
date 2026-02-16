import numpy as np
import matplotlib.pyplot as plt

# ==========================================================
# ПАРАМЕТРЫ
# ==========================================================
l = 1.0
t_max = 1.0
N = 50
J = 100
x_bar = 0.5

h = l / N
tau = t_max / J

x = np.linspace(0, l, N+1)
t = np.linspace(0, t_max, J+1)

n = int(round(x_bar / h))
r = 2 * h**2 / tau

print(f"h = {h:.4f}, tau = {tau:.4f}, r = {r:.4f}")
print(f"x_bar = {x[n]:.4f}, индекс n = {n}")

# ==========================================================
# НАЧАЛЬНЫЕ И ГРАНИЧНЫЕ УСЛОВИЯ
# ==========================================================
def u0(x):
    return np.sin(np.pi * x)

def mu1(t):
    return 0.4

def mu2(t):
    return 0.3

def phi(t):
    return np.exp(-2*t) * np.sin(np.pi/4)

def p_exact(t):
    return np.exp(-1*t) * np.sin(np.pi/8)

y = np.zeros((N+1, J+1))
p = np.zeros(J+1)

# начальное условие
for i in range(N+1):
    y[i,0] = u0(x[i])

# ==========================================================
# ОСНОВНОЙ ЦИКЛ ПО ВРЕМЕНИ
# ==========================================================
for j in range(1, J+1):

    y_prev = y[:, j-1].copy()
    size = N-1

    # --------------------------
    # 1. Система для z
    # --------------------------
    a = np.ones(size)
    c = np.full(size, -(2 + r))
    b = np.ones(size)

    d = np.zeros(size)

    for i in range(1, N):
        idx = i-1
        d_bar = r*y_prev[i] + (y_prev[i+1] - 2*y_prev[i] + y_prev[i-1])
        d[idx] = -d_bar

    d[0] -= mu1(t[j])
    d[-1] -= mu2(t[j])

    # --- метод прогонки ---
    alpha = np.zeros(size)
    beta = np.zeros(size)

    alpha[0] = -b[0] / c[0]
    beta[0] = d[0] / c[0]

    for i in range(1, size):
        denom = c[i] + a[i] * alpha[i-1]
        alpha[i] = -b[i] / denom if i < size-1 else 0
        beta[i] = (d[i] - a[i]*beta[i-1]) / denom

    z_inner = np.zeros(size)
    z_inner[-1] = beta[-1]

    for i in reversed(range(size-1)):
        z_inner[i] = alpha[i]*z_inner[i+1] + beta[i]

    z = np.zeros(N+1)
    z[0] = mu1(t[j])
    z[N] = mu2(t[j])
    z[1:N] = z_inner

    # --------------------------
    # 2. Система для v
    # --------------------------
    d_v = np.zeros(size)
    for i in range(1, N):
        idx = i-1
        d_v[idx] = -h**2 * y_prev[i]

    alpha[:] = 0
    beta[:] = 0

    alpha[0] = -b[0] / c[0]
    beta[0] = d_v[0] / c[0]

    for i in range(1, size):
        denom = c[i] + a[i] * alpha[i-1]
        alpha[i] = -b[i] / denom if i < size-1 else 0
        beta[i] = (d_v[i] - a[i]*beta[i-1]) / denom

    v_inner = np.zeros(size)
    v_inner[-1] = beta[-1]

    for i in reversed(range(size-1)):
        v_inner[i] = alpha[i]*v_inner[i+1] + beta[i]

    v = np.zeros(N+1)
    v[1:N] = v_inner

    # --------------------------
    # 3. Вычисляем p^j
    # --------------------------
    if abs(v[n]) < 1e-12:
        p[j] = p[j-1]
    else:
        p[j] = (phi(t[j]) - z[n]) / v[n]

    # --------------------------
    # 4. Обновляем решение
    # --------------------------
    y[:, j] = z + p[j]*v


# ==========================================================
# ГРАФИКИ
# ==========================================================

# 1) Решение в момент T
plt.figure()
plt.plot(x, y[:, -1])
plt.title("Решение прямой задачи u(x,T)")
plt.xlabel("x")
plt.ylabel("u")
plt.grid()
plt.show()

# 2) Восстановленный коэффициент
plt.figure()
plt.plot(t, p)
plt.title("Восстановленный коэффициент p(t)")
plt.xlabel("t")
plt.ylabel("p")
plt.grid()
plt.show()
