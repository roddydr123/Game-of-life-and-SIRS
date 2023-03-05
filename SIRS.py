import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
from tqdm import tqdm
import timeit


def update_grid(frame, im, grid, grid_size, p_vals):

    # array of random indices for one sweep.
    ijs = np.random.randint(0, grid_size, size=(grid_size**2, 2))

    # array of probabilities for S -> I - > R -> S.
    random_probs = np.random.rand(grid_size**2)

    for step in range(grid_size**2):

        i, j = ijs[step]

        # if susceptible
        if grid[i, j] == 0:
            neighbours = (
                grid[(i - 1) % grid_size, j],
                grid[i, (j - 1) % grid_size],
                grid[i, (j + 1) % grid_size],
                grid[(i + 1) % grid_size, j],
            )
            # if there's an infected next door
            if 1 in neighbours:
                if p_vals[0] > random_probs[step]:
                    grid[i, j] = 1

        # if infected
        elif grid[i, j] == 1:
            if p_vals[1] > random_probs[step]:
                grid[i, j] = 2

        # if recovered
        elif grid[i, j] == 2:
            if p_vals[2] > random_probs[step]:
                grid[i, j] = 0

    if im:
        im.set_data(grid)
    inf_sites = np.sum(grid == 1)

    return im, inf_sites


def visualisation(grid, grid_size, p_vals):

    fig, ax = plt.subplots()

    # Create an image object to display the grid
    im = ax.imshow(grid, cmap="gray")

    # Create the animation object
    ani = animation.FuncAnimation(
        fig, update_grid, frames=1000, fargs=(im, grid, grid_size, p_vals), interval=10
    )

    # Display the animation
    plt.show()


def variance(grid_size):

    p_space = np.arange(0.2, 0.5, 0.005)
    p2 = 0.5
    p3 = 0.5

    f = open("SIRS_data/variance_plot.dat", "w")
    f.write("p1, p2, p3, var(I)\n")
    f.close()

    no_iterations = 10100

    random_grids = np.random.randint(3, size=(len(p_space), grid_size, grid_size))

    for p1, grid in tqdm(zip(p_space, random_grids), total=len(p_space)):
        # for each set of probabilities
        p_vals = [p1, p2, p3]

        inf_sites_list = np.zeros(no_iterations)

        # run 10100 times to get measurements
        for k in range(no_iterations + 100):
            inf_sites = update_grid(None, None, grid, grid_size, p_vals)[1]
            if inf_sites == 0:
                break

            # only take measurements past 100 sweeps
            if k >= 100:
                inf_sites_list[k - 100] = inf_sites

        norm_variance = np.var(inf_sites_list) / grid_size**2

        np.savetxt(f"SIRS_data/var.{p1}.{p2}.{p3}.dat", inf_sites_list)

        f = open("SIRS_data/variance_plot.dat", "a")
        f.write(f"{p1},{p2},{p3},{norm_variance}\n")
        f.close()


def get_bootstrap_error(filename):

    inf_list = np.loadtxt(filename)

    k = 100
    n = int(0.8 * len(inf_list))

    variances = np.zeros(k)

    for i in range(k):
        infecteds = np.random.choice(inf_list, n)
        variances[i] = np.var(infecteds)

    error = np.std(variances / 50**2)
    return error


def plot_variance():
    # load in the variances.
    data = np.genfromtxt(
        "SIRS_data/variance_plot.dat", delimiter=",", skip_header=1, dtype=float
    )
    norm_var_I = np.array(data[:, 3])

    errors = []

    # get the error for each variance using the bootstrap method.
    for p1, p2, p3 in zip(data[:, 0], data[:, 1], data[:, 2]):
        errors.append(get_bootstrap_error(f"SIRS_data/var.{p1}.{p2}.{p3}.dat"))

    # plot variance as a function of p1.
    p1s = np.array(data[:, 0])
    fig, ax = plt.subplots()
    ax.errorbar(p1s, norm_var_I, yerr=errors)
    ax.set_xlabel("P1")
    ax.set_ylabel("normalised variance of infected")
    plt.show()


def phase(grid_size):

    p_space = np.arange(0, 1, 0.05)
    p2 = 0.5

    with open("SIRS_data/infected_plot.dat", "w") as f:
        f.write("p1,p2,p3,average I\n")

    random_grid = np.random.randint(
        3, size=(len(p_space), len(p_space), grid_size, grid_size)
    )

    for p1, super_grid in tqdm(zip(p_space, random_grid), total=len(p_space)):
        for p3, grid in zip(p_space, super_grid):

            # for each set of probabilities
            p_vals = [p1, p2, p3]
            inf_sites_list = np.zeros(1000)

            # run 1100 times to get measurements
            for k in range(1100):
                inf_sites = update_grid(None, None, grid, grid_size, p_vals)[1]
                if inf_sites == 0:
                    break

                # only take measurements past 100 sweeps
                if k >= 100:
                    inf_sites_list[k - 100] = inf_sites

            # save the evolution of that simulation
            np.savetxt(f"SIRS_data/infected.{p1}.{p2}.{p3}.dat", inf_sites_list[::10])

            # save the average infected for that simulation
            with open(f"SIRS_data/infected_plot.dat", "a") as f:
                f.write(f"{p1},{p2},{p3},{np.average(inf_sites_list)/grid_size**2}\n")


def plot_phase():
    data = np.genfromtxt(
        "SIRS_data/infected_plot.dat", delimiter=",", skip_header=1, dtype=float
    )
    d_l = int(len(data) ** (1 / 2))
    av_I = np.array(data[:, 3]).reshape((d_l, d_l))
    fig, ax = plt.subplots()
    ax.imshow(av_I, origin="lower", extent=(0, 1, 0, 1))
    ax.set_xlabel("P1")
    ax.set_ylabel("P3")
    plt.show()


def immunity(grid, grid_size, p_vals):

    # set the percentages of immune people to explore.
    percentage_range = np.arange(0.01, 1, 0.01)

    random_grids = np.random.randint(
        3, size=(len(percentage_range), grid_size, grid_size)
    )

    p1, p2, p3 = p_vals

    # how long should each simulation be?
    no_iterations = 1000

    with open(f"SIRS_data/immunity_plot.{p1}.{p2}.{p3}.dat", "w") as f:
        f.write("percentage, average I\n")

    for percentage, grid in tqdm(
        zip(percentage_range, random_grids), total=len(percentage_range)
    ):

        ### SINGLE SIMULATION ###

        # initialise immune cells.
        ijs = np.random.randint(
            0, grid_size, size=(int(grid_size**2 * percentage), 2)
        )
        for i, j in ijs:
            grid[i, j] = 3

        inf_sites_list = np.zeros(no_iterations)

        for k in range(no_iterations + 100):
            inf_sites = update_grid(None, None, grid, grid_size, p_vals)[1]
            if inf_sites == 0:
                break

            # only take measurements past 100 sweeps
            if k >= 100:
                inf_sites_list[k - 100] = inf_sites

        # save the evolution of this simulation
        np.savetxt(
            f"SIRS_data/immunity.{percentage}.{p1}.{p2}.{p3}.dat", inf_sites_list[::10]
        )

        # save only the average infected from this simulation
        with open(f"SIRS_data/immunity_plot.{p1}.{p2}.{p3}.dat", "a") as f:
            f.write(f"{percentage},{np.average(inf_sites_list)/grid_size**2}\n")

        ### END SINGLE SIMULATION ###


def get_bootstrap_error_immun(filename):

    inf_list = np.loadtxt(filename)

    k = 100
    n = int(0.8 * len(inf_list))

    av_infecteds = np.zeros(k)

    for i in range(k):
        infecteds = np.random.choice(inf_list, n)
        av_infecteds[i] = np.average(infecteds)

    error = np.std(av_infecteds) / 50**2
    return error


def plot_immunity():

    p1, p2, p3 = [0.5, 0.5, 0.5]
    # load in the variances.
    data = np.genfromtxt(
        f"SIRS_data/immunity_plot.{p1}.{p2}.{p3}.dat",
        delimiter=",",
        skip_header=1,
        dtype=float,
    )
    percentages = np.array(data[:, 0])
    inf_list = np.array(data[:, 1])

    errors = []

    # get the error for each percentage's average no. infected using the bootstrap method.
    for percentage in percentages:
        errors.append(
            get_bootstrap_error_immun(
                f"SIRS_data/immunity.{percentage}.{p1}.{p2}.{p3}.dat"
            )
        )

    # plot average infection as a function of percentage immune.
    fig, ax = plt.subplots()
    ax.errorbar(percentages, inf_list, yerr=errors, fmt="x")
    ax.set_xlabel("Proportion immune")
    ax.set_ylabel("Average proportion infected over 1000 iterations")
    plt.show()


def main():
    cmd_args = sys.argv
    if len(cmd_args) == 1:
        print("Usage SIRS.py <mode> <grid_size> <p1> <p2> <p3>")
        sys.exit()

    mode = cmd_args[1]

    if mode == "phase_plot":
        plot_phase()
        return
    elif mode == "var_plot":
        plot_variance()
        return
    elif mode == "immunity_plot":
        plot_immunity()
        return

    grid_size = int(cmd_args[2])

    grid = np.random.randint(3, size=(grid_size, grid_size))

    p_vals = [float(x) for x in cmd_args[3:]]

    if mode == "vis":
        visualisation(grid, grid_size, p_vals)
    elif mode == "phase":
        phase(grid_size)
    elif mode == "var":
        variance(grid_size)
    elif mode == "timer":
        t = timeit.Timer(
            lambda: update_grid(None, None, grid, grid_size, [0.5, 0.5, 0.5])
        )
        print(t.timeit(5))
    elif mode == "immunity":
        immunity(grid, grid_size, p_vals)


main()
