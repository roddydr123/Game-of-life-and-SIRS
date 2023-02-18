import numpy as np
import matplotlib.pyplot as plt

# Set up the grid
grid_size = 50
grid = np.zeros((grid_size, grid_size))

# Initialize the grid with some live cells
num_live_cells = 1000
live_cells = np.random.choice(grid_size**2, num_live_cells, replace=False)
grid[np.unravel_index(live_cells, (grid_size, grid_size))] = 1

# Update the grid based on the Game of Life rules
def update_grid(grid):
    new_grid = np.zeros_like(grid)
    for i in range(grid_size):
        for j in range(grid_size):
            # Count the number of live neighbors
            num_live_neighbors = np.sum(grid[max(i-1, 0):min(i+2, grid_size), max(j-1, 0):min(j+2, grid_size)]) - grid[i, j]
            # Apply the Game of Life rules
            if grid[i, j] == 1 and num_live_neighbors in [2, 3]:
                new_grid[i, j] = 1
            elif grid[i, j] == 0 and num_live_neighbors == 3:
                new_grid[i, j] = 1
    return new_grid

# Visualize the grid
def plot_grid(grid):
    plt.imshow(grid, cmap='binary')

# Run the simulation
num_generations = 100

fig, ax = plt.subplots()
im=ax.imshow(grid, animated=True)

for generation in range(num_generations):
    grid = update_grid(grid)
    fig.canvas.flush_events()
    im=ax.imshow(grid, animated=True)
    ax.draw_artist(im)
    plt.pause(0.011)
