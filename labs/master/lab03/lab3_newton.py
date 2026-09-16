import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parents[3] / "results" / "lab03"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Параметры сетки ────────────────────────────────────────────────────────
L = 1.0;  T = 0.5
J = 500;  N = 50

tau = T / J;  h = L / N
r   = h**2 / tau     # r = h²/τ

x = np.linspace(0, L, N + 1)
t = np.linspace(0, T, J + 1)

# ── Тестовый пример ────────────────────────────────────────────────────────
# p(t) = 1 + 0.5·cos(2πt),  u(x,0) = sin(πx),  μ₁=μ₂=0
p_exact = 1. + 0.5 * np.cos(2. * np.pi * t)
u0 = np.sin(np.pi * x)

# ═══════════════════════════════════════════════════════════════════════════
#  Прогонка (метод Томаса)
# ═══════════════════════════════════════════════════════════════════════════
def thomas(a, b, c, f):
    """Трёхдиагональная система: a·y_{i-1} + b·y_i + c·y_{i+1} = f."""
    n  = len(f)
    b_ = b.copy(); f_ = f.copy()
    for i in range(1, n):
        w     = a[i] / b_[i-1]
        b_[i] -= w * c[i-1]
        f_[i] -= w * f_[i-1]
    sol = np.zeros(n)
    sol[-1] = f_[-1] / b_[-1]
    for i in range(n-2, -1, -1):
        sol[i] = (f_[i] - c[i]*sol[i+1]) / b_[i]
    return sol

# ═══════════════════════════════════════════════════════════════════════════
#  ПРЯМАЯ ЗАДАЧА
#  Неявная схема:  p^j·(y^j − ŷ)/τ = Λ_h y^j
#  => −y_{i-1} + (2 + r·p)·y_i − y_{i+1} = r·p·ŷ_i
# ═══════════════════════════════════════════════════════════════════════════
u = np.zeros((J + 1, N + 1))
u[0, :] = u0.copy()

for j in range(J):
    pj = p_exact[j + 1];  rp = r * pj;  ni = N - 1
    a_ = np.full(ni, -1.); b_ = np.full(ni, 2.+rp); c_ = np.full(ni, -1.)
    u[j+1, 1:N] = thomas(a_, b_, c_, rp * u[j, 1:N])

def trapz(f):
    return h * (np.sum(f) - 0.5*(f[0] + f[N]))

phi = np.array([trapz(u[j]) for j in range(J+1)])

# ═══════════════════════════════════════════════════════════════════════════
#  ОБРАТНАЯ ЗАДАЧА — квазирешение, итерации Ньютона по p
#
#  На каждом шаге j→j+1 находим p^{j+1} из нелинейного уравнения:
#
#    F(p) = trapz( y(p) ) − φ^{j+1} = 0,
#
#  где y(p) — решение прямой задачи при данном p:
#    A(p)·y = r·p·ŷ,   A(p)_{ii} = 2+rp,  A_{i,i±1} = −1
#
#  Производная (дифференцируем A(p)·y = r·p·ŷ по p):
#    A(p)·(dy/dp) = r·(ŷ − y)
#
#  Шаг Ньютона:  p ← p − F(p)/F'(p),  2 прогонки на итерацию.
#  Сходимость квадратичная, ~3–5 итераций до машинной точности.
# ═══════════════════════════════════════════════════════════════════════════
y     = np.zeros((J + 1, N + 1))
p_rec = np.zeros(J + 1)
y[0, :]  = u0.copy()
p_rec[0] = p_exact[0]

def y_of_p(p_val, y_hat_inner):
    """Решить (2+rp)·y_i − y_{i-1} − y_{i+1} = rp·ŷ_i → вернуть y_inner."""
    rp = r * p_val
    return thomas(np.full(N-1,-1.), np.full(N-1, 2.+rp), np.full(N-1,-1.),
                  rp * y_hat_inner)

def F_and_Fp(p_val, y_hat_inner, phi_target):
    """
    F(p)  = trapz(y(p)) − φ
    F'(p) = trapz(dy/dp)

    dy/dp: дифференцируем A(p)·y = rp·ŷ по p:
      A'(p)·y + A(p)·(dy/dp) = r·ŷ
      A(p)·(dy/dp) = r·ŷ − r·I·y  = r·(ŷ − y)
    """
    rp = r * p_val
    a_ = np.full(N-1,-1.); b_ = np.full(N-1, 2.+rp); c_ = np.full(N-1,-1.)
    y_in   = thomas(a_, b_, c_, rp * y_hat_inner)
    dydp   = thomas(a_, b_, c_,  r * (y_hat_inner - y_in))

    y_full = np.zeros(N+1); y_full[1:N] = y_in
    d_full = np.zeros(N+1); d_full[1:N] = dydp
    return trapz(y_full) - phi_target, trapz(d_full)

for j in range(J):
    y_hat_inner = y[j, 1:N]
    phi_j1      = phi[j + 1]
    p_k         = p_rec[j]           # начальное приближение

    for _ in range(60):
        Fv, Fp = F_and_Fp(p_k, y_hat_inner, phi_j1)
        if abs(Fp) < 1e-30:
            break
        delta = Fv / Fp
        p_k  -= delta
        if abs(delta) < 1e-14 * (abs(p_k) + 1e-30):
            break

    p_rec[j+1]  = p_k
    y[j+1, 1:N] = y_of_p(p_k, y_hat_inner)

# ═══════════════════════════════════════════════════════════════════════════
#  ВИЗУАЛИЗАЦИЯ
# ═══════════════════════════════════════════════════════════════════════════
TEAL = "#00d4c8"; AMBER = "#ffb347"; RED = "#ff5555"; WHITE = "#e8e8e8"

fig = plt.figure(figsize=(16, 10), facecolor="#ffffff")
fig.suptitle(
    r"Обратная задача: $\;p(t)\,\partial_t u = \partial_{xx} u$"
    r"$\quad\longrightarrow\quad$ квазирешение  $y = p\cdot w$",
    color=WHITE, fontsize=13, fontfamily="monospace", y=0.97)

gs = gridspec.GridSpec(2, 3, figure=fig,
                       hspace=0.48, wspace=0.38,
                       left=0.06, right=0.97, top=0.91, bottom=0.07)

ax1 = fig.add_subplot(gs[0, 0])
im1 = ax1.imshow(u.T, aspect="auto", origin="lower", cmap="viridis",
                 extent=[0, T, 0, L])
fig.colorbar(im1, ax=ax1, pad=0.02)
ax1.set_title("u(x,t) — прямая задача (точная p)")
ax1.set_xlabel("t"); ax1.set_ylabel("x")

ax2 = fig.add_subplot(gs[0, 1])
im2 = ax2.imshow(y.T, aspect="auto", origin="lower", cmap="viridis",
                 extent=[0, T, 0, L])
fig.colorbar(im2, ax=ax2, pad=0.02)
ax2.set_title("y(x,t) — квазирешение (восст. p)")
ax2.set_xlabel("t"); ax2.set_ylabel("x")

ax3 = fig.add_subplot(gs[0, 2])
im3 = ax3.imshow(np.abs(u-y).T, aspect="auto", origin="lower", cmap="viridis",
                 extent=[0, T, 0, L])
fig.colorbar(im3, ax=ax3, pad=0.02)
ax3.set_title("|u − y| — погрешность поля")
ax3.set_xlabel("t"); ax3.set_ylabel("x")

ax4 = fig.add_subplot(gs[1, 0])
ax4.plot(t, p_exact, color=TEAL,  lw=2,   label="p(t) точная")
ax4.plot(t, p_rec,   color=AMBER, lw=1.5, ls="--", label="p(t) восст.")
ax4.set_title("p(t): точная vs восстановленная")
ax4.set_xlabel("t"); ax4.set_ylabel("p(t)")
ax4.legend(fontsize=8, facecolor="#1a1a1a", labelcolor=WHITE, edgecolor="#333")

ax5 = fig.add_subplot(gs[1, 1])
ax5.semilogy(t, np.abs(p_exact - p_rec) + 1e-16, color=RED, lw=1.5)
ax5.set_title("|p − p_rec|  (лог. масштаб)")
ax5.set_xlabel("t"); ax5.set_ylabel("погрешность")

phi_y = np.array([trapz(y[j]) for j in range(J+1)])
ax6 = fig.add_subplot(gs[1, 2])
ax6.plot(t, phi,   color=TEAL,  lw=2,   label=r"$\phi$ (прямая)")
ax6.plot(t, phi_y, color=AMBER, lw=1.5, ls="--", label=r"$\phi$ (квазирешение)")
ax6.set_title(r"$\phi(t) = \int_0^l u\,dx$")
ax6.set_xlabel("t"); ax6.set_ylabel(r"$\phi$")
ax6.legend(fontsize=8, facecolor="#1a1a1a", labelcolor=WHITE, edgecolor="#333")

plt.savefig(OUTPUT_DIR / "inverse_problem.png",
            dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.show()

print(f"‖p − p_rec‖_∞  = {np.max(np.abs(p_exact - p_rec)):.3e}")
print(f"‖u − y‖_∞      = {np.max(np.abs(u - y)):.3e}")
print(f"‖φ_u − φ_y‖_∞  = {np.max(np.abs(phi - phi_y)):.3e}")
