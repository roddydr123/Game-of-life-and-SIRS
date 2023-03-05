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
    proportion_infected = np.sum(grid == 1)

    return im, proportion_infected


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


def get_bootstrap_error(filename):

    inf_list = np.loadtxt(filename)

    step = 100
    n = int(0.8 * len(inf_list))

    variances = np.zeros(step)

    for i in range(step):
        other_vars = np.random.choice(inf_list, n)
        variances[i] = np.var(other_vars)

    error = np.std(variances / 50**2)
    return error


def variance(grid_size):

    p_space = np.arange(0, 1, 0.05)
    p2 = 0.5

    simulation_length = 10000

    with open("SIRS_data/variance_plot.dat", "w") as f:
        f.write("p1, p2, p3, variance(I)/N\n")

    random_grid = np.random.randint(
        3, size=(len(p_space), len(p_space), grid_size, grid_size)
    )

    for p1, super_grid in tqdm(zip(p_space, random_grid), total=len(p_space)):
        for p3, grid in zip(p_space, super_grid):

            ### BEGIN SINGLE SIMULATION ###

            p_vals = [p1, p2, p3]
            proportion_infected_list = np.zeros(simulation_length)

            for step in range(simulation_length + 100):
                proportion_infected = update_grid(None, None, grid, grid_size, p_vals)[1] / grid_size**2
                if proportion_infected == 0:
                    break

                # only take measurements past 100 sweeps
                if step >= 100:
                    proportion_infected_list[step - 100] = proportion_infected

            proportion_infected_list = proportion_infected_list[::10]

            np.savetxt(f"SIRS_data/variance.{p1}.{p2}.{p3}.dat", proportion_infected_list)

            # save the variance infected for that simulation
            with open(f"SIRS_data/variance_plot.dat", "a") as f:
                f.write(f"{p1},{p2},{p3},{np.var(proportion_infected_list)}\n")
            

            ### SINGLE SIMULATION END ###


def plot_variance(colour):

    data = np.genfromtxt(
            "SIRS_data/variance_plot.dat", delimiter=",", skip_header=1, dtype=float
        )
    
    fig, ax = plt.subplots()

    if colour == 'False':

        # make cut from data at p1,p2,p3 = p1,0.5,0.5
        selected_points = data[data[:,2] == 0.5]

        errors = []

        # get the error for each variance using the bootstrap method.
        for p1, p2, p3 in zip(selected_points[:, 0], selected_points[:, 1], selected_points[:, 2]):
            errors.append(get_jacknife_error(f"SIRS_data/variance.{p1}.{p2}.{p3}.dat"))

        # plot variance as a function of p1.
        p1s = np.array(selected_points[:, 0])
        ax.errorbar(p1s, selected_points[:,3], yerr=errors)
        ax.set_xlabel("P1")
        ax.set_ylabel("normalised variance of infected")
        

    else:
        d_l = int(len(data) ** (1 / 2))
        av_I = np.array(data[:, 3]).reshape((d_l, d_l))
        ax.imshow(av_I, origin="lower", extent=(0, 1, 0, 1))
        ax.set_xlabel("P1")
        ax.set_ylabel("P3")
    
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
            proportion_infected_list = np.zeros(1000)

            # run 1100 times to get measurements
            for step in range(1100):
                proportion_infected = update_grid(None, None, grid, grid_size, p_vals)[1]
                if proportion_infected == 0:
                    break

                # only take measurements past 100 sweeps
                if step >= 100:
                    proportion_infected_list[step - 100] = proportion_infected

            # every 10th measurement to prevent autocorrelation
            proportion_infected_list = proportion_infected_list[::10]

            # save the evolution of that simulation
            np.savetxt(f"SIRS_data/infected.{p1}.{p2}.{p3}.dat", proportion_infected_list)

            # save the average infected for that simulation
            with open(f"SIRS_data/infected_plot.dat", "a") as f:
                f.write(f"{p1},{p2},{p3},{np.average(proportion_infected_list)}\n")


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

        ### SINGLE SIMULATION ###

        # initialise immune cells.
        ijs = np.random.randint(
            0, grid_size, size=(int(grid_size**2 * proportion_immune), 2)
        )
        for i, j in ijs:
            grid[i, j] = 3

        proportion_infected_list = np.zeros(simulation_length)

        for step in range(simulation_length + 100):
            proportion_infected = update_grid(None, None, grid, grid_size, p_vals)[1] / grid_size**2
            if proportion_infected <= 1E-5:
                break

            # only take measurements past 100 sweeps
            if step >= 100:
                proportion_infected_list[step - 100] = proportion_infected

        proportion_infected_list = proportion_infected_list[::10]

        # save the evolution of this simulation, only every 10th measurement to prevent autocorrelation.
        np.savetxt(
            f"SIRS_data/immunity.{proportion_immune}.{p1}.{p2}.{p3}.dat", proportion_infected_list
        )

        # save only the average infected from this simulation
        with open(f"SIRS_data/immunity_plot.{p1}.{p2}.{p3}.dat", "a") as f:
            f.write(f"{proportion_immune},{np.average(proportion_infected_list)}\n")

        ### END SINGLE SIMULATION ###


def get_jacknife_error(filename):

    inf_list = np.loadtxt(filename)

    c = np.var(inf_list)

    av_vars = np.zeros(len(inf_list))

    for i in range(len(inf_list)):
        other_var = np.var(inf_list[np.arange(len(inf_list))!=i])
        # other_vars = np.concatenate((inf_list[:i], inf_list[i:]))
        av_vars[i] = other_var

    error = np.sqrt(np.sum((av_vars - c)**2))
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

    # plot average infection as a function of proportion_immune immune.
    fig, ax = plt.subplots()
    ax.errorbar(percentages, inf_list, fmt="x")
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
        plot_variance(colour=cmd_args[2])
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
