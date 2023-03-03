import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
from scipy.optimize import curve_fit


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


def oscillator_grid(grid_size):
    grid = np.zeros((grid_size, grid_size))
    mid = int(grid_size / 2)
    grid[mid, mid] = 1
    grid[mid + 1, mid] = 1
    grid[mid - 1, mid] = 1
    return grid


def glider_grid(grid_size):
    grid = np.zeros((grid_size, grid_size))
    mid = int(grid_size / 2.3)
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
    grid_size = 50

    # number of runs
    for j in range(190):

        grid = random_grid(grid_size)

        stop = False
        sites_list = []
        i = 0
        # run the simulation
        while stop is False:
            if i > 10000:
                stop = True
            im, active_sites, grid = update_grid(i, None, grid, grid_size, "R")
            sites_list.append(active_sites)

            # if unchanged in last 10 iterations
            if len(set(sites_list[-10:])) == 1 and i > 10:
                stop = True
                print(str(j) + " done")
                with open("equlibrations.dat", "a") as f:
                    f.write(str(i) + "\n")
            i += 1


def plot_equilibrations():
    with open("equlibrations.dat", "r") as f:
        lines = [line for line in f]

    x = list(map(int, lines))
    plt.hist(x, bins=30)
    plt.xlabel("No. sweeps for equilibration")
    plt.ylabel("Count")
    print(len(x))
    plt.show()


def c_o_m(grid, grid_size):

    x_list = []
    y_list = []
    banned = [0, grid_size]
    for i in range(grid_size):
        for j in range(grid_size):
            if grid[i,j] == 1:
                if i in banned or j in banned:
                    return None
                x_list.append(i)
                y_list.append(j)
                
    com_x = np.average(x_list)
    com_y = np.average(y_list)
    return com_x, com_y


def plot_glider_com():
    grid_size = 50
    grid = glider_grid(grid_size)

    com_x_list = []
    com_y_list = []
    times_list = []

    for i in range(200):
        grid = update_grid(None, None, grid, grid_size, "G")[2]
        if i % 10 == 0:
            try:
                x,y = c_o_m(grid, grid_size)
                com_x_list.append(x)
                com_y_list.append(y)
                times_list.append(i)
            except:
                pass
    times_listv = np.array(times_list)[10:15]
    com_x_listv = np.array(com_x_list)[10:15]
    com_y_listv = np.array(com_y_list)[10:15]

    xpopt, xpcov = curve_fit(linear, times_listv, com_x_listv)
    x_velocity = xpopt[0]

    ypopt, ypcov = curve_fit(linear, times_listv, com_y_listv)
    y_velocity = ypopt[0]

    total_v = (x_velocity**2 + y_velocity**2)**(1/2)

    print(f"{total_v} pixels per sweep")

    # plt.plot(times_list, linear(times_list, *xpopt), label="Linear fit")
    # plt.scatter(times_list, com_x_list, label="COM position")
    # plt.xlabel("Time (sweeps)")
    # plt.ylabel("X position")

    # plt.plot(times_list, linear(times_list, *ypopt), label="Linear fit")
    # plt.scatter(times_list, com_y_list, label="COM position")
    # plt.xlabel("Time (sweeps)")
    # plt.ylabel("Y position")

    plt.xlabel("X coordinate")
    plt.ylabel("Y coordinate")
    plt.title(f"Velocity: {round(total_v, 2)} pixels per sweep")
    plt.scatter(com_x_list, com_y_list)
    plt.legend()
    plt.show()


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
# equilibrations()
# plot_equilibrations()
plot_glider_com()
