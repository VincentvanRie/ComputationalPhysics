import numpy as np
import matplotlib.pyplot as plt

"""
In this part the data is analysed 
"""

##TODO:
# plot bins only till L/2


## Specify for which rho/T the analysis is done
T = 1
rho = 0.8

positions = np.load(f"data/positions{T}_{rho}_long.npy", allow_pickle=True)
velocities = np.load(f"data/velocities{T}_{rho}_long.npy", allow_pickle=True)

iterations = len(positions) - 1
num_pressures = 20  # number of times the pressure is measured

N = len(positions[0]["x"])

L = (N / rho) ** (1 / (3))
print("Box size is:", L)
h = 0.01  # 10**-15

bin_size = 0.01

pressure_times = np.linspace(0, iterations, num_pressures).astype(int)
print(pressure_times)


def two_part_distance(p1, p2, positions):
    """
    Take the position dictionary and 2 particles
    as input and return the distance between p1 and p1

    """

    delta_x = (positions["x"][p2] - positions["x"][p1] + L / 2) % L - L / 2
    delta_y = (positions["y"][p2] - positions["y"][p1] + L / 2) % L - L / 2
    delta_z = (positions["z"][p2] - positions["z"][p1] + L / 2) % L - L / 2
    r = ((delta_x) ** 2 + (delta_y) ** 2 + (delta_z) ** 2) ** 0.5

    return r


def pair_correlation_plot(positions):
    """
    input: position dictionary
    output: plot of the histogram
    """

    particle_distances = np.array([])

    for i in range(N):
        for j in range(i + 1, N):
            r = two_part_distance(i, j, positions)
            particle_distances = np.append(particle_distances, r)

    y = len(particle_distances)
    # print(particle_distances)

    # plt.scatter(np.arange(y),particle_distances)

    bins = L / bin_size

    plt.hist(particle_distances, np.arange(0, L, bin_size))
    plt.show()


def average_hist(all_positions, analyse_times, bin_size):
    """
    input: all positions and the times to be analysed
    output: pair correlation histogram of the average hist,bins
    """
    plot = True  # plots the avg_hist if true
    num_bins = int(L / bin_size)
    total_hist = np.zeros(num_bins)

    for t in analyse_times:
        positions = all_positions[t]
        particle_distances = np.array([])
        for i in range(N):
            for j in range(i + 1, N):
                r = two_part_distance(i, j, positions)

                particle_distances = np.append(particle_distances, r)

        hist, bins = np.histogram(particle_distances, bins=num_bins)
        total_hist += hist

    avg_hist = total_hist / len(analyse_times)

    # print(bins[:-1])
    pair_correlation_g = (
        (2 * L**3) / (N * (N - 1)) * avg_hist / (4 * np.pi * bins[:-1] ** 2 * bin_size)
    )

    if plot:
        plt.bar(bins[:-1], pair_correlation_g, width=np.diff(bins), edgecolor="black")
        plt.xlabel("Distance")
        plt.ylabel("Frequency")
        plt.title(f"Average Pair Correlation Histogram with rho = {rho},T = {T} ")
        plt.grid(True)
        plt.show()

    return pair_correlation_g, bins


def pressure(positions, rho):
    interaction_term = np.array([])

    for i in range(N):
        for j in range(i + 1, N):
            r = two_part_distance(i, j, positions)

            potential_der = 2 * (
                -24 * (1 / r) ** 13 + 12 * (1 / r) ** 7
            )  # TODO plus or min?
            # print(potential_der)

            interaction_term = np.append(interaction_term, r * potential_der)
            # print('r is ',r)
            # print('potential der:', potential_der)
    expec_value = np.sum(interaction_term)
    print(expec_value)
    # print(expec_value)
    pressure = (
        T * rho - rho / (6 * N) * expec_value
    )  # TODO something wrong with dimensionless units? expec value too low

    return pressure


def plot_fcc_lattice(positions):
    # Plot the FCC lattice
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    colors = np.random.rand(10, 3)

    """ keep this part chatgpt
    for i in random_indices:
    
        # Scatter plot for FCC lattice
        ax.scatter(positions['x'][i], positions['y'][i], positions['z'][i], c='b', marker='o', s=5, alpha=0.7) #, label='FCC lattice'
    
        plt.quiver(system.positions["x"][i], system.positions["y"][i],system.positions["z"][i],
                        system.velocities["x"][i], system.velocities["y"][i],system.velocities["z"][i], length = 0.1)
        """
    ax.scatter(
        positions["x"],
        positions["y"],
        positions["z"],
        c="b",
        marker="o",
        s=5,
        alpha=0.7,
    )  # , label='FCC lattice'

    # Set labels and title for the plot
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("FCC Lattice")

    """
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_zlim(0, L)
    """

    # Add legend and grid to the plot
    ax.legend()
    ax.grid(True, linestyle="dotted", linewidth=0.5, alpha=0.5)

    # Display the plot
    plt.show()


def main():
    # pair-correlation
    """
    plot_fcc_lattice(positions[1])

    hist_plots = np.arange(len(positions))


    times = np.arange(len(positions))
    times = np.arange(20,40)

    average_hist(positions,times,bin_size)
    """
    pressures = np.array([])
    for i in pressure_times:
        # pair_correlation_plot(positions[i])

        # print('pressure is:',i,pressure(positions[i], 0.8))
        pressures = np.append(pressures, pressure(positions[i], 0.8))

    plt.plot(pressure_times, pressures)
    # plt.ylim(0,np.max(pressures)+1)
    plt.title(f"rho is {rho} and T is {T}")
    plt.show()

    pressure_error = np.std(pressures)
    print("The error in the pressure is given by:", pressure_error)

    # pressure
    #  for i in hist_plots():
    #    print(pressure(positions[i], 0.8))


if __name__ == "__main__":
    main()
