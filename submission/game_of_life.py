import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
from scipy.optimize import curve_fit


def linear(x, m, c):
    # function for fitting a straight line.
    return (m * x) + c


def update_grid(frame, im, grid, grid_size, mode):

    # count the number of active elements.
    active_sites = np.sum(grid)

    # avoid overwriting.
    new_grid = grid.copy()

    # loop over each cell in the grid.
    for i in range(grid_size):
        for j in range(grid_size):
            # live neighbour count
            num_neighbours = (
                grid[(i - 1) % grid_size, (j - 1) % grid_size]
                + grid[(i - 1) % grid_size, j]
                + grid[(i - 1) % grid_size, (j + 1) % grid_size]
                + grid[i, (j - 1) % grid_size]
                + grid[i, (j + 1) % grid_size]
                + grid[(i + 1) % grid_size, (j - 1) % grid_size]
                + grid[(i + 1) % grid_size, j]
                + grid[(i + 1) % grid_size, (j + 1) % grid_size]
            )
            # if the cell is alive...
            if grid[i, j] == 1 and (num_neighbours < 2 or num_neighbours > 3):
                new_grid[i, j] = 0

            # if cell is dead
            elif grid[i, j] == 0 and num_neighbours == 3:
                new_grid[i, j] = 1

    # set the grid to new values simultaneously
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
    # put a small ocillator in the middle of the grid
    grid = np.zeros((grid_size, grid_size))
    mid = int(grid_size / 2)
    grid[mid, mid] = 1
    grid[mid + 1, mid] = 1
    grid[mid - 1, mid] = 1
    return grid


def glider_grid(grid_size):
    # put a glider near the centre of the grid
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
    ani = animation.FuncAnimation(
        fig, update_grid, frames=1000, fargs=(im, grid, grid_size, mode), interval=10
    )

    # Display the animation
    plt.show()


def equilibrations():
    # run without visualisation
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


def c_o_m(grid:np.ndarray, grid_size:int) -> tuple[float,float]:
    """Find the location of the centre of mass.

    Args:
        grid (np.ndarray): Any grid.
        grid_size (int): Grid dimensions.

    Returns:
        tuple: x and y positions of the centre of mass.
    """

    x_list = []
    y_list = []
    banned = [0, grid_size]
    for i in range(grid_size):
        for j in range(grid_size):
            if grid[i, j] == 1:
                if i in banned or j in banned:
                    return None
                x_list.append(i)
                y_list.append(j)

    com_x = np.average(x_list)
    com_y = np.average(y_list)
    return com_x, com_y


def plot_glider_com():
    """Make a line plot of the glider's centre of mass. This can be a plot of
    x or y position vs number of sweeps, or CoM plotted in x and y once per sweep.
    """

    grid_size = 50
    grid = glider_grid(grid_size)

    com_x_list = []
    com_y_list = []
    times_list = []

    for i in range(200):
        grid = update_grid(None, None, grid, grid_size, "G")[2]
        if i % 10 == 0:
            try:
                x, y = c_o_m(grid, grid_size)
                com_x_list.append(x)
                com_y_list.append(y)
                times_list.append(i)
            except:
                pass

    times_list = np.array(times_list)
    com_x_list = np.array(com_x_list)
    com_y_list = np.array(com_y_list)
    times_listv = np.array(times_list)[10:15]
    com_x_listv = np.array(com_x_list)[10:15]
    com_y_listv = np.array(com_y_list)[10:15]

    xpopt, xpcov = curve_fit(linear, times_listv, com_x_listv)
    x_velocity = xpopt[0]

    ypopt, ypcov = curve_fit(linear, times_listv, com_y_listv)
    y_velocity = ypopt[0]

    total_v = (x_velocity**2 + y_velocity**2) ** (1 / 2)

    print(f"{total_v} pixels per sweep")

    # plt.plot(times_listv, linear(times_listv, *xpopt), label=f"y = {round(xpopt[0], 2)} x + {round(xpopt[1], 1)}")
    # plt.scatter(times_listv, com_x_listv, label="COM position")
    # plt.xlabel("Time (sweeps)")
    # plt.ylabel("X position")
    # np.savetxt("GoL_data/x_velocity_fit.dat", np.stack((times_listv, com_x_listv), axis=1))

    # plt.plot(times_listv, linear(times_listv, *ypopt), label=f"y = {round(ypopt[0], 2)} x + {round(ypopt[1], 1)}")
    # plt.scatter(times_listv, com_y_listv, label="COM position")
    # plt.xlabel("Time (sweeps)")
    # plt.ylabel("Y position")
    # np.savetxt("GoL_data/y_velocity_fit.dat", np.stack((times_listv, com_y_listv), axis=1))

    plt.xlabel("X coordinate")
    plt.ylabel("Y coordinate")
    plt.title(f"Velocity: {round(total_v, 2)} pixels per sweep")
    plt.scatter(com_x_list, com_y_list)
    np.savetxt(
        "GoL_data/glider_position.dat",
        np.stack((com_x_list, com_y_list), axis=1),
        header="glider x position, glider y position",
    )
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


main()
# equilibrations()
# plot_equilibrations()
# plot_glider_com()
