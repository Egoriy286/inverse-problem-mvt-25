from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

N = 51
mesh = UnitSquareMesh(N, N)
V = FunctionSpace(mesh, 'Lagrange', 1)
f = Expression("beta*x[0]*x[1]", beta=1, degree=1)
gamma = 100
T = 0.1

def p(t):
    return -1000*t/(1 + np.exp(gamma*(t - 0.5*T)))

p0_val = Constant(p(0))
u0_val = Constant(1.0)
u0 = interpolate(u0_val, V)
p0 = interpolate(p0_val, V)
bcs = []

M = 100
tau = T / M

u_trial = TrialFunction(V)
v = TestFunction(V)

a = (1 / tau) * u_trial * v * dx + inner(grad(u_trial), grad(v)) * dx - p0 * u_trial * v * dx
L = (1 / tau) * u0 * v * dx + f * v * dx

u = Function(V)

# Координаты узлов сетки
coords = mesh.coordinates()          # (n_nodes, 2)
cells  = mesh.cells()                # (n_cells, 3)  — треугольники

# Хранилище: список массивов shape (n_nodes,), по одному на каждый временной шаг
u_snapshots = []
time_steps  = []

t = 0
while t < T - 1e-14:
    t += tau

    solve(a == L, u, bcs)

    # Сохраняем значения в узлах как numpy-массив
    u_np = u.compute_vertex_values(mesh)   # shape: (n_nodes,)
    u_snapshots.append(u_np.copy())
    time_steps.append(t)

    u0.assign(u)
    pt = interpolate(Constant(p(t)), V)
    p0.assign(pt)

# Собираем в единый массив: (M, n_nodes)
u_array = np.array(u_snapshots)
time_steps = np.array(time_steps)

print(f"u_array.shape = {u_array.shape}")   # (100, n_nodes)

# ── Визуализация ──────────────────────────────────────────────────────────────
x = coords[:, 0]
y = coords[:, 1]
triangulation = tri.Triangulation(x, y, cells)
data = {
    'x':x,
    'y':y,
    'cells':cells,
    'triang': triangulation,
    'u': u_array
}

np.save('laba6.npy', data)


# Выбираем нужный шаг (последний, или любой другой индекс)
step_idx = -1   # последний шаг
u_plot   = u_array[step_idx]

fig, ax = plt.subplots(figsize=(7, 6))
tcf = ax.tricontourf(triangulation, u_plot, levels=20, cmap='viridis')
fig.colorbar(tcf, ax=ax, label='u')
ax.set_title(f'u(x,y) at t = {time_steps[step_idx]:.4f}')
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.tight_layout()
plt.savefig('contourf_result.png', dpi=150)