"""
Лабораторная работа: Численное решение обратной задачи
для параболического уравнения методом конечных элементов.

Задача:
    du/dt - Δu - p(t)u = f(x,t),  x ∈ Ω, t ∈ (0,T]
    u|∂Ω = 0
    u(x,0) = u0(x)
    ∫_Ω η(x) u(x,t) dx = φ(t)   — интегральное переопределение

Метод:
    Неявная схема + разложение y = z + p^j * ω на каждом шаге.

Тестовый пример (точное решение):
    u_exact(x,y,t) = exp(t) * sin(πx) * sin(πy)
    p_exact(t)     = sin(t)   (неизвестный коэффициент)
    η(x,y)         = sin(πx) * sin(πy)

Из уравнения выражается f(x,t).
"""

from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parents[4] / "results" / "lab06"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Параметры ─────────────────────────────────────────────────────────────────
N     = 32          # число ячеек по каждому направлению
T     = 1.0         # конечное время
M     = 100         # число временных шагов
tau   = T / M       # шаг по времени
pi    = np.pi

# ── Сетка и пространство ──────────────────────────────────────────────────────
mesh = UnitSquareMesh(N, N)
V    = FunctionSpace(mesh, 'Lagrange', 1)
bc   = DirichletBC(V, Constant(0.0), 'on_boundary')

# ── Точное решение и вспомогательные функции ──────────────────────────────────
# u_exact = exp(t) * sin(πx)*sin(πy)
# Δu = -2π² * exp(t) * sin(πx)*sin(πy)
# du/dt = exp(t) * sin(πx)*sin(πy)
# => f = du/dt - Δu - p(t)*u
#      = exp(t)*sin(πx)*sin(πy) * [1 + 2π² - p(t)]

def p_exact(t):
    return np.sin(t)

def phi_exact(t):
    """φ(t) = ∫_Ω η·u_exact dx = exp(t) * ∫_Ω sin²(πx)sin²(πy)dx = exp(t)/4"""
    return np.exp(t) / 4.0

# Начальное условие u0 = u_exact(x,y,0) = sin(πx)*sin(πy)
u0_expr = Expression('sin(pi*x[0])*sin(pi*x[1])', pi=pi, degree=4)

# Весовая функция η(x,y) = sin(πx)*sin(πy)
eta_expr = Expression('sin(pi*x[0])*sin(pi*x[1])', pi=pi, degree=4)
eta_h    = interpolate(eta_expr, V)

# ── Инициализация ─────────────────────────────────────────────────────────────
y_old = interpolate(u0_expr, V)  # ŷ = u0

# Проверка согласованности: φ(0) = ∫_Ω η·u0 dx
phi0_check = assemble(eta_h * y_old * dx)
print(f"Согласованность: φ(0) = {phi_exact(0):.6f}, ∫η·u0 dx = {phi0_check:.6f}")

# ── Вариационная форма (общая матрица A = M + τK) ─────────────────────────────
u_tr = TrialFunction(V)
v    = TestFunction(V)

a_form = (u_tr * v + tau * inner(grad(u_tr), grad(v))) * dx

# ── Хранение результатов ──────────────────────────────────────────────────────
coords      = mesh.coordinates()
cells_arr   = mesh.cells()
n_nodes     = coords.shape[0]

u_snapshots = np.zeros((M + 1, n_nodes))
p_computed  = np.zeros(M + 1)
p_reference = np.zeros(M + 1)
time_arr    = np.zeros(M + 1)

# Начальный момент
u_snapshots[0] = y_old.compute_vertex_values(mesh)
p_computed[0]  = p_exact(0.0)     # согласованное начальное значение
p_reference[0] = p_exact(0.0)
time_arr[0]    = 0.0

# ── Основной цикл ─────────────────────────────────────────────────────────────
z     = Function(V)
omega = Function(V)
A     = assemble(a_form)
bc.apply(A)

for j in range(1, M + 1):
    t_j = j * tau

    # Правая часть f^j: f = exp(t)*(1+2π²-p(t))*sin(πx)*sin(πy)
    coeff_f = np.exp(t_j) * (1.0 + 2.0 * pi**2 - p_exact(t_j))
    f_expr  = Expression('c*sin(pi*x[0])*sin(pi*x[1])',
                         c=coeff_f, pi=pi, degree=4)
    f_h = interpolate(f_expr, V)

    # Правые части двух СЛАУ (матрица A одинакова):
    # Az = Mŷ + τ·f^j
    # Aω = τ·Mŷ
    Myold = assemble(y_old * v * dx)
    tauf  = assemble(tau * f_h * v * dx)
    tauMy = assemble(tau * y_old * v * dx)

    rhs_z = Myold.copy(); rhs_z += tauf
    rhs_w = tauMy.copy()

    bc.apply(rhs_z)
    bc.apply(rhs_w)

    solve(A, z.vector(), rhs_z)
    solve(A, omega.vector(), rhs_w)

    # Вычисление p^j = (φ^j - ∫η·z dx) / ∫η·ω dx
    Iz  = assemble(eta_h * z     * dx)
    Iw  = assemble(eta_h * omega * dx)
    phi_j = phi_exact(t_j)

    if abs(Iw) < 1e-14:
        raise RuntimeError(f"Шаг {j}: знаменатель ≈ 0, задача вырождена!")

    p_j = (phi_j - Iz) / Iw

    # Обновление: y = z + p^j · ω
    y_new = Function(V)
    y_new.vector()[:] = z.vector()[:] + p_j * omega.vector()[:]

    # Сохранение
    u_snapshots[j] = y_new.compute_vertex_values(mesh)
    p_computed[j]  = p_j
    p_reference[j] = p_exact(t_j)
    time_arr[j]    = t_j

    y_old.assign(y_new)

    if j % 10 == 0:
        err_p = abs(p_j - p_exact(t_j))
        print(f"  t = {t_j:.3f}  p_computed = {p_j:+.6f}  "
              f"p_exact = {p_exact(t_j):+.6f}  |err| = {err_p:.2e}")

# ── Ошибки ────────────────────────────────────────────────────────────────────
# Ошибка по p(t)
err_p_max  = np.max(np.abs(p_computed - p_reference))
err_p_l2   = np.sqrt(tau * np.sum((p_computed - p_reference)**2))

# Ошибка по u в момент T (L2-норма через errornorm)
u_exact_T = Expression('exp(T)*sin(pi*x[0])*sin(pi*x[1])',
                        T=T, pi=pi, degree=4)
u_ex_h  = interpolate(u_exact_T, V)
y_final = Function(V)
y_final.assign(y_old)          # y_old уже содержит y^M после цикла
err_u_l2 = errornorm(u_ex_h, y_final, norm_type='L2', degree_rise=3)

print(f"\n{'='*55}")
print(f"  N = {N},  M = {M},  τ = {tau:.4f}")
print(f"  max|p_h - p_ex|  = {err_p_max:.4e}")
print(f"  L2-норма по p    = {err_p_l2:.4e}")
print(f"  L2-норма по u(T) = {err_u_l2:.4e}")
print(f"{'='*55}\n")

# ── Визуализация ──────────────────────────────────────────────────────────────
triangulation = tri.Triangulation(coords[:, 0], coords[:, 1], cells_arr)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.suptitle(r'Обратная задача: восстановление $p(t)$ для параболического уравнения',
             fontsize=13)

# 1) p(t): вычисленный vs точный
ax = axes[0]
ax.plot(time_arr, p_reference, 'k-',  lw=2,   label=r'$p_{exact}(t) = \sin t$')
ax.plot(time_arr, p_computed,  'r--', lw=1.5, label=r'$p_h(t)$')
ax.set_xlabel('t');  ax.set_ylabel('p(t)')
ax.set_title('Восстановленный коэффициент $p(t)$')
ax.legend();  ax.grid(True, alpha=0.4)

# 2) Численное решение u_h в момент T
ax = axes[1]
tcf = ax.tricontourf(triangulation, u_snapshots[M], levels=20, cmap='plasma')
fig.colorbar(tcf, ax=ax, label='$u_h$')
ax.set_title(f'$u_h(x,y,T)$,  T={T}')
ax.set_xlabel('x');  ax.set_ylabel('y')

# 3) Ошибка |u_h - u_exact| в момент T
u_ex_vals = u_ex_h.compute_vertex_values(mesh)
u_num_vals = y_final.compute_vertex_values(mesh)
err_vals   = np.abs(u_num_vals - u_ex_vals)
ax = axes[2]
tcf2 = ax.tricontourf(triangulation, err_vals, levels=20, cmap='Reds')
fig.colorbar(tcf2, ax=ax, label='$|u_h - u_{ex}|$')
ax.set_title(r'Ошибка $|u_h - u_{ex}|$ при $t=T$')
ax.set_xlabel('x');  ax.set_ylabel('y')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'inverse-result-generated.png', dpi=150)
plt.show()
print("График сохранён: inverse_result.png")
