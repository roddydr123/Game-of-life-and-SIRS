import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys



def check_active(grid):
    active = np.sum(grid)
    with open("active_sites.dat", "a") as f:
        f.write(str(active) + "\n")


def update_grid(i, im, grid, grid_size):
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
    # Update the grid
    im.set_data(new_grid)
    grid[:] = new_grid[:]
    return im


def random_grid(grid_size):
    grid = np.random.randint(2, size=(grid_size, grid_size))
    return grid


def oscillator_grid(grid_size):
    grid = np.zeros((grid_size, grid_size))
    mid = int(grid_size / 2)
    grid[mid, mid] = 1
    grid[mid + 1, mid] = 1
    grid[mid - 1, mid] = 1
    return grid


def glider_grid(grid_size):
    grid = np.zeros((grid_size, grid_size))
    mid = int(grid_size / 2)
    grid[mid, mid] = 1
    grid[mid + 1, mid] = 1
    grid[mid - 1, mid] = 1
    grid[mid + 1, mid + 1] = 1
    grid[mid, mid + 2] = 1
    return grid


def visualisation(grid, grid_size, monitor_sites):

    fig, ax = plt.subplots()

    # Create an image object to display the grid
    im = ax.imshow(grid)

    active_list = []

    # Create the animation object
    ani = animation.FuncAnimation(fig, update_grid, frames=1000, fargs=(im, grid, grid_size), interval=10)

    # Display the animation
    plt.show()


def main():
    cmd_args = sys.argv
    if len(cmd_args) != 2:
        print("Usage game_of_life.py <mode> (Random, Oscillator, Glider)")
        sys.exit()
    mode = cmd_args[1]

    grid_size = 50
    monitor_sites = True

    if mode == "R":
        grid = random_grid(grid_size)
    elif mode == "O":
        grid = oscillator_grid(grid_size)
    elif mode == "G":
        grid = glider_grid(grid_size)
    else:
        print("invalid mode")
        sys.exit()
    visualisation(grid, grid_size, monitor_sites)

main()
