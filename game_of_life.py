import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys



def check_active(grid, i, stop=False):
    active = np.sum(grid)
    with open("active_sites.dat", "a") as f:
        f.write(str(active) + "\n")
    if i % 10 == 0 and i != 0:
        with open("active_sites.dat", "r") as f:
            lines = [line.strip() for line in f]
            if len(set(lines[-10:])) <= 1:
                with open("equlilibration_time.dat", "a") as f:
                    f.write(str(i) + "\n")
                stop = True
    return stop



def update_grid(i, im, grid, grid_size, mode):

    # count active cells
    if mode == "R":
        stop = check_active(grid, i)
    else:
        stop = False

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

    if type(im) != str:
        im.set_data(new_grid)
        if stop:
            print(f"stopped after {i} iterations")
            sys.exit(0)
        return im
    else:
        return stop


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


def visualisation(grid, grid_size, mode):

    fig, ax = plt.subplots()

    # Create an image object to display the grid
    im = ax.imshow(grid)

    # Create the animation object
    ani = animation.FuncAnimation(fig, update_grid, frames=1000, fargs=(im, grid, grid_size, mode), interval=10)

    # Display the animation
    plt.show()


def equilibrations():

    im = "equil"
    grid_size = 50
    grid = random_grid(grid_size)

    for j in range(10):

        f = open("equlilibration_time.dat", "w")
        f.close()
        i = 0
        stop = False
        while stop is False:
            stop = update_grid(i, im, grid, grid_size, "R")
            i += 1



def main():
    cmd_args = sys.argv
    if len(cmd_args) != 2:
        print("Usage game_of_life.py <mode> (Random, Oscillator, Glider)")
        sys.exit()
    mode = cmd_args[1]

    grid_size = 50

    if mode == "R":
        grid = random_grid(grid_size)
    elif mode == "O":
        grid = oscillator_grid(grid_size)
    elif mode == "G":
        grid = glider_grid(grid_size)
    else:
        print("invalid mode")
        sys.exit()
    visualisation(grid, grid_size, mode)

# main()
equilibrations()