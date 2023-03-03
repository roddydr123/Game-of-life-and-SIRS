import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
from tqdm import tqdm



def update_grid(frame, im, grid, grid_size, p_vals):

    for step in range(grid_size**2):

        i, j = np.random.randint(0, grid_size, size=(2,))

        # if susceptible
        if grid[i,j] == 0:
            neighbours = (
                grid[(i-1)%grid_size,j],
                grid[i,(j-1)%grid_size],
                grid[i,(j+1)%grid_size],
                grid[(i+1)%grid_size,j]
            )
            # if there's an infected next door
            if 1 in neighbours:
                if p_vals[0] > np.random.rand():
                    grid[i,j] = 1

        # if infected
        elif grid[i,j] == 1:
            if p_vals[1] > np.random.rand():
                grid[i,j] = 2

        # if recovered
        elif grid[i,j] == 2:
            if p_vals[2] > np.random.rand():
                grid[i,j] = 0

    if im:
        im.set_data(grid)
    inf_sites = np.sum(grid == 1)

    return im, inf_sites


def visualisation(grid, grid_size, p_vals):

    fig, ax = plt.subplots()

    # Create an image object to display the grid
    im = ax.imshow(grid, cmap="gray")

    # Create the animation object
    ani = animation.FuncAnimation(fig, update_grid, frames=1000, fargs=(im, grid, grid_size, p_vals), interval=10)

    # Display the animation
    plt.show()


def phase(grid_size):

    p_space = np.arange(0, 1, 0.05)
    p2 = 0.5

    with open("colour_plot.dat", "w") as f:
        f.write("p1,p2,p3,average I\n")


    for p1 in p_space:
        for p3 in p_space:
            # for each set of probabilities
            p_vals = [p1, p2, p3]
            inf_sites_list = np.zeros(1000)

            grid = np.random.randint(3, size=(grid_size, grid_size))

            # run 1100 times to get measurements
            for k in tqdm(range(1100)):
                inf_sites = update_grid(None, None, grid, grid_size, p_vals)[1]
                if inf_sites == 0:
                    break

                # only take measurements past 100 sweeps
                if k >= 100:
                    inf_sites_list[k-100] = inf_sites

            with open("colour_plot.dat", "a") as f:
                f.write(f"{p1},{p2},{p3},{np.average(inf_sites_list)}\n")



def plot_phase():
    colour_grid = np.zeros((len(p_space), len(p_space)))
    colour_grid[i,j] = np.average(inf_sites_list)


def main():
    cmd_args = sys.argv
    if len(cmd_args) == 1:
        print("Usage SIRS.py <grid_size> <p1> <p2> <p3>")
        sys.exit()

    grid_size = int(cmd_args[1])

    grid = np.random.randint(3, size=(grid_size, grid_size))

    p_vals = [float(x) for x in cmd_args[2:]]

    if len(cmd_args) == 5:
        visualisation(grid, grid_size, p_vals)
    elif len(cmd_args) == 2:
        phase(grid_size)

main()
