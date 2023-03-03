import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys


def linear(x, m, c):
    return (m * x) + c


def update_grid(frame, im, grid, grid_size, mode):

    active_sites = np.sum(grid)

    # Copy the current grid to avoid overwriting values
    new_grid = grid.copy()
    # Loop over each cell in the grid
    for i in range(grid_size):
        for j in range(grid_size):
            # Count the number of live neighbours for the current cell
            num_neighbours = (
                grid[(i-1)%grid_size,(j-1)%grid_size] +
                grid[(i-1)%grid_size,j] +
                grid[(i-1)%grid_size,(j+1)%grid_size] +
                grid[i,(j-1)%grid_size] +
                grid[i,(j+1)%grid_size] +
                grid[(i+1)%grid_size,(j-1)%grid_size] +
                grid[(i+1)%grid_size,j] +
                grid[(i+1)%grid_size,(j+1)%grid_size]
            )
            # Apply the rules of the game of life
            if grid[i,j] == 1 and (num_neighbours < 2 or num_neighbours > 3):
                new_grid[i,j] = 0
            elif grid[i,j] == 0 and num_neighbours == 3:
                new_grid[i,j] = 1

    grid[:] = new_grid[:]

    try:
        im.set_data(new_grid)
    except:
        pass
    return im, active_sites, grid


def random_grid(grid_size):
    grid = np.random.randint(2, size=(grid_size, grid_size))
    return grid


def visualisation(grid, grid_size, mode):

    fig, ax = plt.subplots()

    # Create an image object to display the grid
    im = ax.imshow(grid)

    # Create the animation object
    ani = animation.FuncAnimation(fig, update_grid, frames=1000, fargs=(im, grid, grid_size), interval=10)

    # Display the animation
    plt.show()


def main():
    cmd_args = sys.argv
    if len(cmd_args) != 5:
        print("Usage SIRS.py <grid_size> <p1> <p2> <p3>")
        sys.exit()

    grid_size = cmd_args[1]

    grid = random_grid(grid_size)

    visualisation(grid, grid_size)

main()
