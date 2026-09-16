import numpy as np; import matplotlib.pyplot as plt

# ── Параметры ──────────────────────────────────────────────────────────────
L=1; T=1; N=50; mu1=0.; mu2=0.
h=L/N; x=np.linspace(0,L,N+1)

# Точное v(t) и начальное u₀(x)
v_exact_fn = lambda t: 1.0 + 0.5*np.sin(np.pi*t)   # любая >0 функция
u0 = np.sin(np.pi*x)                                 # начальные данные

# Датчик: фиксированная точка x̄ (условие (4): u(x̄,t)=ψ(t))
n_s = N//3                                           # индекс узла датчика
# ψ(t) будем брать из точного решения на ходу

# ── Хранилища (динамические списки, т.к. τⱼ адаптивный) ──────────────────
Y   = [u0.copy()]          # сеточное решение по слоям
V   = []                   # восстановленные vʲ
TAU = []                   # шаги τⱼ
T_  = [0.0]                # моменты времени

# ── Прогонка: решить -y_{i-1} + c*y_i - y_{i+1} = d_i ───────────────────
# возвращает y[0..N], y[0]=bc0, y[N]=bcN
def thomas(d, c, bc0, bcN):
    n=len(d); k=np.zeros(n); e=np.zeros(n); y=np.zeros(n)
    k[1]=0; e[1]=bc0
    for i in range(1,n-1):
        denom=c-k[i]; k[i+1]=1/denom; e[i+1]=(d[i]+e[i])/denom
    y[n-1]=bcN
    for i in range(n-2,-1,-1): y[i]=k[i+1]*y[i+1]+e[i+1]
    y[0]=bc0; y[n-1]=bcN; return y

# ── Начальный τ₀ из v_exact(0) по формуле (8): τ = h/(2v) ────────────────
v0_init = v_exact_fn(0.)
tau0 = h/(2*v0_init)
V.append(v0_init); TAU.append(tau0)

# ── Главный цикл по времени ───────────────────────────────────────────────
while T_[-1] < T - 1e-10:
    yj  = Y[-1]; tj = T_[-1]
    tau = TAU[-1]; vj = V[-1]

    # ── Этап 1: уравнение переноса (метод характеристик, доска: ȳᵢ=y^{j-1}_{i-1}) ──
    # τ_{j-1}/h = 1/(2v) → сдвиг на «полузла»; аппроксимация вверх по потоку
    # ȳᵢ = yᵢ - (τ/h)*vʲ*(yᵢ - yᵢ₋₁)  — upwind, учитывает дробный сдвиг
    ybar = np.zeros(N+1)
    ybar[0] = mu1
    for i in range(1,N+1):
        ybar[i] = yj[i] - (tau/h)*vj*(yj[i]-yj[i-1])
    ybar[N] = mu2

    # ── Этап 2: уравнение диффузии (неявная схема, прогонка) ─────────────
    # h²/(2τ)*(yⁿ⁺¹ᵢ - ȳᵢ) = y_{i+1} - 2y_i + y_{i-1}
    # → -y_{i-1} + c*y_i - y_{i+1} = ȳᵢ,   c = 2 + h²/(2τ) · ... 
    # (из доски: h²*(yᵢ-ȳᵢ)/(2τ) = y_{i+1}-2yᵢ+y_{i-1})
    r2 = h**2/(2*tau); cv = 2 + r2
    d  = np.zeros(N+1)
    for i in range(1,N): d[i] = r2*ybar[i]   # правая часть
    y_new = thomas(d, cv, mu1, mu2)

    # ── Восстановление vʲ⁺¹ из условия (7) на доске: ─────────────────────
    # vʲ = (h/τ_{j-1}) * (ȳ_{n+j+1} - y_{n+j}) / (y_{n+j+1} - ψ_{n+j})
    # Используем датчик: ψ = u(x_s, t_{j+1}) — точное значение в узле n_s
    t_new = tj + tau
    psi   = v_exact_fn(t_new) * 0 + np.sin(np.pi*x[n_s])*np.exp(-(np.pi**2)*t_new)
    # (в реальной задаче ψ — измерение; здесь берём из точного решения)

    num = ybar[n_s] - yj[n_s]           # ȳ_{n_s} - y_{n_s}^j
    den = y_new[n_s] - psi              # y_{n_s}^{j+1} - ψ^{j+1}
    if abs(den) > 1e-13:
        v_new = (h/tau) * (num/den)
    else:
        v_new = vj                       # если знаменатель мал — держим прежнее

    v_new = max(v_new, 1e-3)            # v должна быть >0 (физика)

    # ── Новый адаптивный шаг (8): τⱼ = h/(2vʲ) ──────────────────────────
    tau_new = h/(2*v_new)
    # Ограничиваем, чтобы не выйти за T
    if t_new + tau_new > T: tau_new = T - t_new

    Y.append(y_new); V.append(v_new)
    TAU.append(tau_new); T_.append(t_new)

    if t_new >= T: break

# ── Сравнение с точным v(t) ───────────────────────────────────────────────
T_ = np.array(T_); V = np.array(V)
v_e = np.array([v_exact_fn(t) for t in T_])

# ── Графики ──────────────────────────────────────────────────────────────
fig,axes=plt.subplots(2,2,figsize=(12,8))

# Решение u(x,t) как тепловая карта
Y_arr = np.array(Y)
im0=axes[0,0].imshow(Y_arr,aspect='auto',origin='lower',
    extent=[0,L,0,T_[-1]]); plt.colorbar(im0,ax=axes[0,0])
axes[0,0].set_title('u(x,t) — численное решение'); axes[0,0].set_xlabel('x'); axes[0,0].set_ylabel('t')

# Адаптивные шаги τⱼ
axes[0,1].plot(T_[:-1],TAU,'.-',markersize=3)
axes[0,1].set_title('Адаптивный шаг τⱼ = h/(2vʲ)'); axes[0,1].set_xlabel('t'); axes[0,1].grid()

# v(t): точная vs восстановленная
axes[1,0].plot(T_,v_e,label='v точн.',linewidth=2)
axes[1,0].plot(T_,V,'--',label='v восст.',linewidth=1.5)
axes[1,0].legend(); axes[1,0].set_title('v(t)'); axes[1,0].grid()

# Погрешность
err=np.abs(V-v_e); err[err<1e-16]=1e-16
axes[1,1].semilogy(T_,err,'r'); axes[1,1].set_title('|v-v_e| погрешность')
axes[1,1].grid(); axes[1,1].set_xlabel('t')

plt.tight_layout(); plt.show()
print(f"Шагов: {len(T_)-1}, max|v-v_e|={err.max():.2e}, T_final={T_[-1]:.4f}")