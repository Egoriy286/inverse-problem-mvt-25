import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
from inverse_problems.numerics import solve_tridiagonal

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
    return solve_tridiagonal(a, b, c, f)

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
#  ОБРАТНАЯ ЗАДАЧА — квазирешение  y = z + p·w
#
#  Уравнение на слое j+1 при неизвестном p = p^{j+1}:
#    −y_{i-1} + (2 + r·p)·y_i − y_{i+1} = r·p·ŷ_i
#
#  Разложение  y = z + p·w  (линейно по p, т.к. отбрасываем p²·r·w):
#    z ≡ 0  (нулевые границы, нулевая правая часть)
#    w: −w_{i-1} + (2+r)·w_i − w_{i+1} = r·ŷ_i,   w₀=w_N=0
#
#  Из φ^{j+1} = trapz(p·w) = p·trapz(w):
#    p^{j+1} = φ^{j+1} / trapz(w)
# ═══════════════════════════════════════════════════════════════════════════
y     = np.zeros((J + 1, N + 1))
p_rec = np.zeros(J + 1)
y[0, :]  = u0.copy()
p_rec[0] = p_exact[0]

ni  = N - 1
a_w = np.full(ni, -1.); b_w = np.full(ni, 2.+r); c_w = np.full(ni, -1.)

for j in range(J):
    w_inner = thomas(a_w, b_w, c_w, r * y[j, 1:N])
    w = np.zeros(N + 1); w[1:N] = w_inner

    denom    = trapz(w)
    p_rec[j+1] = phi[j+1] / denom if abs(denom) > 1e-14 else p_rec[j]
    y[j+1, 1:N] = p_rec[j+1] * w_inner

# ═══════════════════════════════════════════════════════════════════════════
#  ВИЗУАЛИЗАЦИЯ
# ═══════════════════════════════════════════════════════════════════════════
plt.rcParams.update({
    "figure.facecolor": "#0d0d0d", "axes.facecolor": "#111111",
    "text.color": "#e8e8e8", "axes.labelcolor": "#e8e8e8",
    "xtick.color": "#888", "ytick.color": "#888",
    "axes.edgecolor": "#333", "axes.grid": True,
    "grid.color": "#222", "grid.linewidth": 0.6,
    "font.family": "monospace", "axes.titlesize": 11, "axes.labelsize": 9,
})
TEAL = "#00d4c8"; AMBER = "#ffb347"; RED = "#ff5555"; WHITE = "#e8e8e8"

fig = plt.figure(figsize=(16, 10), facecolor="#0d0d0d")
fig.suptitle(
    r"Обратная задача: $\;p(t)\,\partial_t u = \partial_{xx} u$"
    r"$\quad\longrightarrow\quad$ квазирешение  $y = p\cdot w$",
    color=WHITE, fontsize=13, fontfamily="monospace", y=0.97)

gs = gridspec.GridSpec(2, 3, figure=fig,
                       hspace=0.48, wspace=0.38,
                       left=0.06, right=0.97, top=0.91, bottom=0.07)

ax1 = fig.add_subplot(gs[0, 0])
im1 = ax1.imshow(u.T, aspect="auto", origin="lower", cmap="plasma",
                 extent=[0, T, 0, L])
fig.colorbar(im1, ax=ax1, pad=0.02)
ax1.set_title("u(x,t) — прямая задача (точная p)")
ax1.set_xlabel("t"); ax1.set_ylabel("x")

ax2 = fig.add_subplot(gs[0, 1])
im2 = ax2.imshow(y.T, aspect="auto", origin="lower", cmap="plasma",
                 extent=[0, T, 0, L])
fig.colorbar(im2, ax=ax2, pad=0.02)
ax2.set_title("y(x,t) — квазирешение (восст. p)")
ax2.set_xlabel("t"); ax2.set_ylabel("x")

ax3 = fig.add_subplot(gs[0, 2])
im3 = ax3.imshow(np.abs(u-y).T, aspect="auto", origin="lower", cmap="inferno",
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

plt.savefig(OUTPUT_DIR / "inverse_problem_laba3.png",
            dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.show()

print(f"max|p - p_rec|  = {np.max(np.abs(p_exact - p_rec)):.3e}")
print(f"max|u - y|      = {np.max(np.abs(u - y)):.3e}")
print(f"max|phi_u - phi_y| = {np.max(np.abs(phi - phi_y)):.3e}")
