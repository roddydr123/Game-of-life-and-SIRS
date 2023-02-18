import numpy as np
import matplotlib.pyplot as plt
from numba import njit



class GameOfLife:

    def __init__(self) -> None:
        self.grid_size = 50
        self.grid = np.zeros((self.grid_size, self.grid_size))

        # Initialize the self.grid with some live cells
        num_live_cells = 1000
        live_cells = np.random.choice(self.grid_size**2, num_live_cells, replace=False)
        self.grid[np.unravel_index(live_cells, (self.grid_size, self.grid_size))] = 1

    # Update the self.grid based on the Game of Life rules
    @njit
    def update_grid(self):
        new_grid = np.zeros_like(self.grid)
        for i in range(self.grid_size):
            for j in range(self.grid_size):
                # Count the number of live neighbors
                num_live_neighbors = np.sum(self.grid[max(i-1, 0):min(i+2, self.grid_size), max(j-1, 0):min(j+2, self.grid_size)]) - self.grid[i, j]
                # Apply the Game of Life rules
                if self.grid[i, j] == 1 and num_live_neighbors in [2, 3]:
                    new_grid[i, j] = 1
                elif self.grid[i, j] == 0 and num_live_neighbors == 3:
                    new_grid[i, j] = 1
        return new_grid