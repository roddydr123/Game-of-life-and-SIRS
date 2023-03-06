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

    # do a single sweep
    for step in range(grid_size**2):

        # fetch indices of point to update
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
    infected = np.sum(grid == 1)

    return im, infected


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


def colour_variance(grid_size):

    # set the range of probabilities to explore
    p_space = np.arange(0, 1, 0.05)
    p2 = 0.5

    # max. iterations per simulation
    simulation_length = 1000

    with open("SIRS_data/variance_line_plot.dat", "w") as f:
        f.write("p1, p2, p3, variance(I)/N\n")

    # make a set of starting grids for all sims
    random_grid = np.random.randint(
        3, size=(len(p_space), len(p_space), grid_size, grid_size)
    )

    # loop through all parameters to explore
    for p1, super_grid in tqdm(zip(p_space, random_grid), total=len(p_space)):
        for p3, grid in zip(p_space, super_grid):

            ### BEGIN SINGLE SIMULATION ###

            p_vals = [p1, p2, p3]
            infected_list = np.zeros(simulation_length)

            # evolve the simulation until it reaches zero or end
            for step in range(simulation_length + 100):
                infected = (
                    update_grid(None, None, grid, grid_size, p_vals)[1]
                )
                if infected == 0:
                    break

                # only take measurements past 100 sweeps
                if step >= 100:
                    infected_list[step - 100] = infected

            # only take every 10th measurement to avoid autocorrelation
            infected_list = infected_list[::10]

            # save total infected elements
            np.savetxt(
                f"SIRS_data/variance_line.{p1}.{p2}.{p3}.dat", infected_list
            )

            # save the variance infected for that simulation
            with open(f"SIRS_data/variance_line_plot.dat", "a") as f:
                f.write(f"{p1},{p2},{p3},{np.var(infected_list)}\n")

            ### SINGLE SIMULATION END ###


def line_variance(grid_size:int):
    """Same as colour_variance but takes a cut though at constant p2,p3

    Args:
        grid_size (int): N so that the grid is NxN.
    """

    p_space = np.arange(0.2, 0.5, 0.01)
    p2 = 0.5
    p3 = 0.5

    simulation_length = 10000

    with open("SIRS_data/variance_line_plot.dat", "w") as f:
        f.write("p1, p2, p3, variance(I)/N\n")

    random_grid = np.random.randint(
        3, size=(len(p_space), len(p_space), grid_size, grid_size)
    )

    for p1, grid in tqdm(zip(p_space, random_grid[0]), total=len(p_space)):

        ### BEGIN SINGLE SIMULATION ###

        p_vals = [p1, p2, p3]
        infected_list = np.zeros(simulation_length)

        for step in range(simulation_length + 100):
            infected = (
                update_grid(None, None, grid, grid_size, p_vals)[1]
            )
            if infected == 0:
                break

            # only take measurements past 100 sweeps
            if step >= 100:
                infected_list[step - 100] = infected

        infected_list = infected_list[::10]

        np.savetxt(
            f"SIRS_data/variance_line.{p1}.{p2}.{p3}.dat", infected_list
        )

        # save the variance infected for that simulation
        with open(f"SIRS_data/variance_line_plot.dat", "a") as f:
            f.write(f"{p1},{p2},{p3},{np.var(infected_list)}\n")

        ### SINGLE SIMULATION END ###


def plot_variance(colour:bool, grid_size:int):
    """Visualise the variance data, either in line plot or colour plot.

    Args:
        colour (bool): Colour plot (True) or line plot (False).
        grid_size (int): Grid dimension.
    """

    fig, ax = plt.subplots()

    # make a line plot
    if colour == "False":

        data = np.genfromtxt(
        "SIRS_data/variance_line_plot.dat", delimiter=",", skip_header=1, dtype=float
        )

        # make cut from data at p1,p2,p3 = p1,0.5,0.5
        selected_points = data[data[:, 2] == 0.5]

        errors = []

        # get the error for each variance using the bootstrap method.
        for p1, p2, p3 in zip(
            selected_points[:, 0], selected_points[:, 1], selected_points[:, 2]
        ):
            errors.append(get_jacknife_error(f"SIRS_data/variance_line.{p1}.{p2}.{p3}.dat"))

        # normalise by no. elements in grid
        errors = np.array(errors) / grid_size**2
        norm_variance = selected_points[:, 3] / grid_size**2
        # plot variance as a function of p1.
        p1s = np.array(selected_points[:, 0])
        ax.errorbar(p1s, norm_variance, yerr=errors, fmt="x")
        ax.set_xlabel("P1")
        ax.set_ylabel("normalised variance of infected")

    else:
        data = np.genfromtxt(
        "SIRS_data/variance_plot.dat", delimiter=",", skip_header=1, dtype=float
        )
        d_l = int(len(data) ** (1 / 2))
        # fetch var(I) and normalise
        av_I = np.array(data[:, 3]).reshape((d_l, d_l)).T / grid_size**2
        im = ax.imshow(av_I, origin="lower", extent=(0, 1, 0, 1))
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("Var(I) / N")
        ax.set_xlabel("P1")
        ax.set_ylabel("P3")

    plt.show()


def phase(grid_size:int):
    """Generate data for a phase plot of average infected vs p1 and p3.

    Args:
        grid_size (int): Grid dimension.
    """

    p_space = np.arange(0, 1, 0.01)
    p2 = 0.5

    simulation_length = 1000

    with open("SIRS_data/infected_plot.dat", "w") as f:
        f.write("p1,p2,p3,average I\n")

    random_grid = np.random.randint(
        3, size=(len(p_space), len(p_space), grid_size, grid_size)
    )

    for p1, super_grid in tqdm(zip(p_space, random_grid), total=len(p_space)):
        for p3, grid in zip(p_space, super_grid):

            # for each set of probabilities
            p_vals = [p1, p2, p3]
            proportion_infected_list = np.zeros(simulation_length)

            # run 1100 times to get measurements
            for step in range(simulation_length + 100):
                proportion_infected = (
                    update_grid(None, None, grid, grid_size, p_vals)[1] / grid_size**2
                )
                if proportion_infected == 0:
                    break

                # only take measurements past 100 sweeps
                if step >= 100:
                    proportion_infected_list[step - 100] = proportion_infected

            # save the evolution of that simulation
            np.savetxt(
                f"SIRS_data/infected.{p1}.{p2}.{p3}.dat", proportion_infected_list
            )

            # save the average infected for that simulation
            with open(f"SIRS_data/infected_plot.dat", "a") as f:
                f.write(f"{p1},{p2},{p3},{np.average(proportion_infected_list)}\n")


def plot_phase():
    """Make a colour plot of average infected vs p1 and p3.
    """
    data = np.genfromtxt(
        "SIRS_data/infected_plot.dat", delimiter=",", skip_header=1, dtype=float
    )
    d_l = int(len(data) ** (1 / 2))
    av_I = np.array(data[:, 3]).reshape((d_l, d_l)).T
    fig, ax = plt.subplots()
    im = ax.imshow(av_I, origin="lower", extent=(0, 1, 0, 1))
    ax.set_xlabel("P1")
    ax.set_ylabel("P3")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Average infected per grid element")
    plt.show()


def immunity(grid:np.ndarray, grid_size:int, p_vals:list):
    """Find average infected for different percentages of the population being immune.

    Args:
        grid (np.ndarray): An initial grid with no immune elements.
        grid_size (int): Grid dimensions.
        p_vals (list): Values of p1,p2,p3.
    """

    # set the percentages of immune people to explore.
    proportion_immune_range = np.arange(0.01, 1, 0.01)

    random_grids = np.random.randint(
        3, size=(len(proportion_immune_range), grid_size, grid_size)
    )

    p1, p2, p3 = p_vals

    # how long should each simulation be?
    simulation_length = 1000

    with open(f"SIRS_data/immunity_plot.{p1}.{p2}.{p3}.dat", "w") as f:
        f.write("immune %, <infected> %\n")

    for proportion_immune, grid in tqdm(
        zip(proportion_immune_range, random_grids), total=len(proportion_immune_range)
    ):

        ### BEGIN SINGLE SIMULATION ###

        # initialise immune cells.
        ijs = np.random.randint(
            0, grid_size, size=(int(grid_size**2 * proportion_immune), 2)
        )
        for i, j in ijs:
            grid[i, j] = 3

        proportion_infected_list = np.zeros(simulation_length)

        for step in range(simulation_length + 100):
            proportion_infected = (
                update_grid(None, None, grid, grid_size, p_vals)[1] / grid_size**2
            )
            if proportion_infected <= 1e-5:
                break

            # only take measurements past 100 sweeps
            if step >= 100:
                proportion_infected_list[step - 100] = proportion_infected

        proportion_infected_list = proportion_infected_list[::10]

        # save the evolution of this simulation, only every 10th measurement to prevent autocorrelation.
        np.savetxt(
            f"SIRS_data/immunity.{proportion_immune}.{p1}.{p2}.{p3}.dat",
            proportion_infected_list,
        )

        # save only the average infected from this simulation
        with open(f"SIRS_data/immunity_plot.{p1}.{p2}.{p3}.dat", "a") as f:
            f.write(f"{proportion_immune},{np.average(proportion_infected_list)}\n")

        ### END SINGLE SIMULATION ###


def get_jacknife_error(filename:str) -> float:
    """Use jacknife method to find errors on variance data points.

    Args:
        filename (str): File to load raw data to be resampled.

    Returns:
        float: the error on a single point.
    """

    inf_list = np.loadtxt(filename)

    c = np.var(inf_list)

    av_vars = np.zeros(len(inf_list))

    for i in range(len(inf_list)):
        other_var = np.var(inf_list[np.arange(len(inf_list)) != i])
        av_vars[i] = other_var

    error = np.sqrt(np.sum((av_vars - c) ** 2))
    return error


def plot_immunity():
    """Make a line plot of average infected vs immunity in the population.
    """

    p1, p2, p3 = [0.5, 0.5, 0.5]

    data = np.genfromtxt(
        f"SIRS_data/immunity_plot.{p1}.{p2}.{p3}.dat",
        delimiter=",",
        skip_header=1,
        dtype=float,
    )
    percentages = np.array(data[:, 0])
    inf_list = np.array(data[:, 1])

    # plot average infection as a function of proportion_immune immune.
    fig, ax = plt.subplots()
    ax.errorbar(percentages, inf_list, fmt="x")
    ax.set_xlabel("Proportion immune")
    ax.set_ylabel("Average proportion infected over 1000 iterations")
    plt.show()


def main():
    """Evaluate command line args to choose a function.
    """

    cmd_args = sys.argv
    if len(cmd_args) == 1:
        print("Usage SIRS.py <mode> <grid_size> <p1> <p2> <p3>")
        sys.exit()

    mode = cmd_args[1]

    if mode == "phase_plot":
        plot_phase()
        return
    elif mode == "var_plot":
        plot_variance(cmd_args[2], int(cmd_args[3]))
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
    elif mode == "line_var":
        line_variance(grid_size)
    elif mode == "colour_var":
        colour_variance(grid_size)
    elif mode == "timer":
        t = timeit.Timer(
            lambda: update_grid(None, None, grid, grid_size, [0.5, 0.5, 0.5])
        )
        print(t.timeit(5))
    elif mode == "immunity":
        immunity(grid, grid_size, p_vals)


main()
