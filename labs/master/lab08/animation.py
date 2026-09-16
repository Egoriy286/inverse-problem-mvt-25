import numpy as np
import matplotlib.pyplot as plt
from dolfin import *
from matplotlib.animation import FuncAnimation, PillowWriter


def gr1_gif(V, sol_list, levels=20, fs=(6.5, 4), filename="solution.gif", fps=20):

    fig, ax = plt.subplots()
    u = Function(V)

    u.vector()[:] = sol_list[0]
    p = plot(u)
    cbar = fig.colorbar(p)

    def update(i):
        ax.collections.clear()
        u.vector()[:] = sol_list[i]
        plot(u)

    ani = FuncAnimation(fig, update, frames=len(sol_list), blit=False)

    ani.save(filename, writer=PillowWriter(fps=fps))

    plt.show()