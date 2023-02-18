import numpy as np
import matplotlib.pyplot as plt
from numba import njit
from simulation import GameOfLife


sim = GameOfLife()
fig, ax = plt.subplots()
im=ax.imshow(sim.grid, animated=True)

while True:
    grid = sim.update_grid()
    fig.canvas.flush_events()
    im=ax.imshow(grid, animated=True)
    ax.draw_artist(im)
    plt.pause(0.011)
