"""
Обратная задача для уравнения конвекции-диффузии
=================================================
    ∂u/∂t = ∂²u/∂x² − v(t)·∂u/∂x,  x∈(0,1), t∈(0,T]
    u(0,t) = u(1,t) = 0,   u(x,0) = sin(πx)
    φ(t) = ∫₀¹ u(x,t) dx              ← наблюдение

Разностная схема (Кранк–Николсон по диффузии, явный upwind по конвекции):
    (yᵢʲ − yᵢʲ⁻¹)/τ = (yᵢ₋₁ʲ − 2yᵢʲ + yᵢ₊₁ʲ)/h² − v·(yᵢʲ⁻¹ − yᵢ₋₁ʲ⁻¹)/h

Матрица СЛАУ: −yᵢ₋₁ʲ + (2+r)·yᵢʲ − yᵢ₊₁ʲ = rhs_i,  r = h²/τ

Ключ обратной задачи — линейность yʲ⁺¹ = z + v·w по v:
    φ(z + v·w) = φ_obs  ⟹  v = (φ_obs − ∫z) / ∫w
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parents[4] / "results" / "lab05"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ── параметры ──────────────────────────────────────────────────────────────
L, T = 2.0*np.pi, 0.5
N, J = 100, 200       # сетка: N — по x, J — по t
mu1, mu2 = 0.0, 0.0  # граничные условия
h, tau = L/N, T/J
r = h**2 / tau        # параметр схемы r = h²/τ
c = 2.0 + r           # диагональный элемент трёхдиагональной матрицы
x = np.linspace(0, L, N+1)
t = np.linspace(0, T, J+1)
v_exact = np.sin(2*np.pi*t)*10   # истинная скорость (для теста)

# ── прогонка (Thomas algorithm) ────────────────────────────────────────────
def tdma(rhs):
    """Решение −y[i-1] + c·y[i] − y[i+1] = rhs[i],  i=1..N-1."""
    n = N - 1
    k, e = np.zeros(n), np.zeros(n)
    k[0] = 1.0 / c
    e[0] = rhs[0] / c
    for i in range(1, n):
        d = c - k[i-1]
        k[i] = 1.0 / d
        e[i] = (rhs[i] + e[i-1]) / d
    y = np.zeros(N+1)
    y[0] = mu1; y[N] = mu2; y[N-1] = e[n-1]
    for i in range(n-2, -1, -1):
        y[i+1] = k[i+1]*y[i+2] + e[i+1]
    return y

def build_rhs(yj):
    """
    Правая часть СЛАУ = z-часть + v·w-часть.
    rhs_z_i = r·yᵢ + ГУ,   rhs_w_i = −h·(yᵢ − yᵢ₋₁)
    """
    rz = r * yj[1:N].copy()
    rz[0] += mu1; rz[-1] += mu2          # граничные условия
    rw = -h * (yj[1:N] - yj[:N-1])      # конвективный вклад
    return rz, rw

intg = lambda y: h * (y.sum() - 0.5*(y[0] + y[-1]))   # трапеции

# ── ПРЯМАЯ задача ──────────────────────────────────────────────────────────
u = np.zeros((J+1, N+1))
u[0] = np.sin(2*np.pi * x/L) #np.sin(np.pi * x)
for j in range(J):
    rz, rw = build_rhs(u[j])
    u[j+1] = tdma(rz + v_exact[j]*rw)   # СЛАУ с известным v

phi = np.array([intg(u[j]) for j in range(J+1)])   # наблюдение φ(t)

# ── ОБРАТНАЯ задача ────────────────────────────────────────────────────────
y = np.zeros((J+1, N+1))
v = np.zeros(J+1)
y[0] = u[0].copy()
v[0] = v_exact[0]

for j in range(J):
    rz, rw = build_rhs(y[j])
    z = tdma(rz)           # решение без v
    w = tdma(rw)           # вклад от v=1
    phi_w = intg(w)
    if abs(phi_w) > 1e-12:
        v[j+1] = (phi[j+1] - intg(z)) / phi_w
    else:
        v[j+1] = v[j]
    y[j+1] = z + v[j+1]*w
    y[j+1, 0] = mu1; y[j+1, -1] = mu2

err_v = np.abs(v - v_exact)
err_u = np.abs(y - u)
print(f"Погрешность v: max={err_v.max():.4e},  O(τ)={tau:.4f}")
print(f"Погрешность u: max={err_u.max():.2e}")

# ── визуализация ────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
fig.suptitle(r"Обратная задача: $\partial_t u = \partial_{xx}u - v(t)\,\partial_x u$",
             fontsize=13, fontweight='bold')

for ax, data, title in zip(
        [axes[0,0], axes[0,1]], [y, u],
        ["y(x,t) — восстановленное", "u(x,t) — прямая задача"]):
    im = ax.imshow(data, aspect='auto', origin='lower',
                   extent=[0, L, 0, T], cmap='viridis')
    plt.colorbar(im, ax=ax)
    ax.set_title(title); ax.set_xlabel("x"); ax.set_ylabel("t")

axes[1,0].plot(t, v_exact, 'b-',  lw=2.5, label=r"$v_{exact}=\sin t$")
axes[1,0].plot(t, v,       'r--', lw=2,   label=r"$v_{rec}$")
axes[1,0].set_title("Скорость переноса v(t)")
axes[1,0].set_xlabel("t"); axes[1,0].legend(fontsize=10); axes[1,0].grid(alpha=0.4)

axes[1,1].semilogy(t, err_v + 1e-16, 'r-', lw=1.5)
axes[1,1].set_title(f"$|v_{{rec}} - v_{{exact}}|$,  max = {err_v.max():.2e}")
axes[1,1].set_xlabel("t"); axes[1,1].grid(alpha=0.4)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "result-claud.png", dpi=150, bbox_inches='tight')
print("График сохранён: result.png")

# ── анимация ────────────────────────────────────────────────────────────────
fig2, ax2 = plt.subplots(figsize=(7, 4))
ax2.set(xlim=(0,L), ylim=(u.min()*1.1, u.max()*1.1))
ax2.set_xlabel("x"); ax2.set_ylabel("u"); ax2.grid(alpha=0.4)
ax2.set_title(r"Эволюция $u(x,t)$")
lu, = ax2.plot([], [], 'b-',  lw=2, label="u прямая")
ly, = ax2.plot([], [], 'r--', lw=2, label="y обратная")
txt = ax2.text(0.03, 0.90, '', transform=ax2.transAxes, fontsize=11)
ax2.legend()
skip = max(1, J // 60)

def frame(f):
    j = f * skip
    lu.set_data(x, u[j]); ly.set_data(x, y[j])
    txt.set_text(f"t = {t[j]:.3f},  v = {v[j]:.3f}")
    return lu, ly, txt

anim = FuncAnimation(fig2, frame, frames=J//skip, blit=True)
anim.save(OUTPUT_DIR / "laba5-claud.gif", writer=PillowWriter(fps=15))
print("Анимация сохранена: laba5.gif")
plt.show()
