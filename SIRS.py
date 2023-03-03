import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys



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

    im.set_data(grid)

    return im


def visualisation(grid, grid_size, p_vals):

    fig, ax = plt.subplots()

    # Create an image object to display the grid
    im = ax.imshow(grid)

    # Create the animation object
    ani = animation.FuncAnimation(fig, update_grid, frames=1000, fargs=(im, grid, grid_size, p_vals), interval=10)

    # Display the animation
    plt.show()


def main():
    cmd_args = sys.argv
    if len(cmd_args) != 5:
        print("Usage SIRS.py <grid_size> <p1> <p2> <p3>")
        sys.exit()

    grid_size = int(cmd_args[1])

    grid = np.random.randint(3, size=(grid_size, grid_size))

    p_vals = [float(x) for x in cmd_args[2:]]

    visualisation(grid, grid_size, p_vals)

main()
