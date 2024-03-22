import math
import random
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.animation import FuncAnimation
from functools import partial
import numpy as np
import copy
import sys

# class Box:
#     def __init__(self, L, N):
#         self.L = L
#         self.N = N

"""
Natan verandert: N vervangen door len(positions["x"]) waar handig 
                 Functie descripties erbij gezet
                 def main() opgeschoond, lattice initialised met velocity in lattice functie
                 aparte functie inside the class gemaakt om te plotten
                 2D weggehaald
                 
"""


"""
TO DO: 
    wat zijn de juiste waardes voor rescalig lambda?
    pressure and pair correlation in juiste waardes
    
"""

# TODO energy in terms of epsilon:


# %matplotlib inline


z_dimension = True
boltzmann = 1.38064852 * 10**-23
m = 39.948 * 1.66053906660 * 10**-27
epsilon = 119.8 * boltzmann
sigma = 3.405 * 10**-10

N = 108  # number of particles

T = 0.5
rho = 1.2

L = (N / rho) ** (1 / 3)
# L = 3
h = 0.001  # 10**-15

timesteps = 200  # number of iterations
t_end = timesteps * h

save_measurement = (
    True  # if this is true write the data to an external file for further analysis
)

# set the values for lambda renormalization
frequency_rescale = 10
number_of_rescales = 5

renorm_times = np.arange(
    0, t_end, frequency_rescale * h
)  # at these times lambda renormalization is performed
renorm_times = renorm_times[:number_of_rescales]
print(renorm_times)
show_system = False  # plot the system at each iteration

# num_


# print(np.arange(0, t_end, h))


class System_of_particles:
    def __init__(self, positions, velocities):
        self.positions = positions  # [[],[]],[],[]] # {"x": [2, 4, 13, 314 3, 234 2], "y": [[],[],[]]}
        self.velocities = velocities
        self.forces = {"x": [], "y": [], "z": []}

    def Change_positions(self):
        """
        Function that updates the system to the next timestep.
        """
        previous_positions = {"x": np.zeros(N), "y": np.zeros(N), "z": np.zeros(N)}
        dimensions = ["x", "y", "z"]

        for key in dimensions:
            previous_positions[key] = copy.deepcopy(self.positions[key])
            # self.previous_position = [self.particles[i]["position"] for i in range(len(self.particles))]

        current_force = Calculate_force(previous_positions)
        self.forces = current_force

        # current_force = [self.Calculate_force(self.particles[i]["position"]) for i in range(len(self.particles))]

        for key in dimensions:
            self.positions[key] = (
                self.positions[key]
                + self.velocities[key] * h
                + current_force[key] * h**2 / (2)
            )
            self.positions[key] = self.positions[key] % L

        # self.positions = [self.particles[i]["position"] + self.particles[i]["velocity"] * h + current_force[i] * h**2 / (2*m) for i in range(len(self.particles))]

        # x(t) = self.previous_position
        force = Calculate_force(self.positions)

        for key in dimensions:
            self.velocities[key] = self.velocities[key] + (
                force[key] + current_force[key]
            ) * h / (2)

        # print(time)

        print(time)
        if any(
            np.isclose(time, renorm_time) for renorm_time in renorm_times
        ):  # time == h % 10
            print("yeah i rescaled at time", time)
            velocities = np.sqrt(
                self.velocities["x"] ** 2
                + self.velocities["y"] ** 2
                + self.velocities["z"] ** 2
            )

            lambda_ = Lambda(velocities)
            for key in dimensions:
                self.velocities[key] = lambda_ * self.velocities[key]
            # print(self.velocities["y"])

        # return self

    def system_plotter(self, title):
        """
        Function to plot the current system. Input is the title for the plot
        and the output is a plotted figure.
        """
        positions = self.positions
        velocities = self.velocities
        # forces = self.forces

        if plt.fignum_exists(1):
            fig = plt.figure(1)
            fig.clf()
        else:
            fig = plt.figure(1)
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(
            positions["x"],
            positions["y"],
            positions["z"],
            c="b",
            marker="o",
            s=5,
            alpha=0.7,
        )  # , label='FCC lattice'
        plt.quiver(
            positions["x"],
            positions["y"],
            positions["z"],
            velocities["x"],
            velocities["y"],
            velocities["z"],
            length=0.5,
            normalize=True,
            color="red",
        )
        # Set labels and title for the plot
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title("FCC Lattice")

        ax.set_xlim(0, L)
        ax.set_ylim(0, L)
        ax.set_zlim(0, L)

        plt.title(title)

        # Add legend and grid to the plot
        # ax.legend()
        ax.grid(True, linestyle="dotted", linewidth=0.5, alpha=0.5)

        # Display the plot
        plt.pause(0.0001)


def Calculate_force(positions):
    """
    This function calculates the force that each particle experiences

    input: dictionary with x,y,z positions {x: ... y: ... z:...}
    output: dictionary with the forces in the x,y,z direction of
            each particle {x: ... y: ... z:...}

    """

    force_dict = {"x": np.zeros(N), "y": np.zeros(N), "z": np.zeros(N)}

    for i in range(N):  # i = p1, j = p2
        # cell_particle_i
        for j in range(N):
            # cell_particle_j
            if i != j:  # and cell_particle_i = cell_particle_j:
                delta_x = (positions["x"][j] - positions["x"][i] + L / 2) % L - L / 2
                delta_y = (positions["y"][j] - positions["y"][i] + L / 2) % L - L / 2
                delta_z = (positions["z"][j] - positions["z"][i] + L / 2) % L - L / 2

                r = ((delta_x) ** 2 + (delta_y) ** 2 + (delta_z) ** 2) ** 0.5

                if r == 0:
                    print("Zero")

                force = -24 * (2 * (1 / r) ** 13 - (1 / r) ** 7)

                force_dict["x"][i] += force * delta_x / r
                force_dict["y"][i] += force * delta_y / r
                force_dict["z"][i] += force * delta_z / r

    return force_dict


def Calculate_potential(positions):
    """
    input: Dictionary with x,y,z positions {x: ... y: ... z:...}
    This function calcuates the potential between all the particles
    outptu: Array with the potential energy of every particle
    """

    potential_list = []

    N = len(positions["x"])

    for i in range(N):
        # cell_particle_i
        potential = 0

        for j in range(i + 1, N):
            # cell_particle_j
            r = two_part_distance(i, j, positions)
            if np.isclose(r, 0):
                print("Zero")

            potential += 4 * ((1 / r) ** 12 - (1 / r) ** 6)

        potential_list.append(potential)

    return potential_list


def Calculate_kinetic(velocities):
    """
    input: Dictionary with x,y,z velocities {x: ... y: ... z:...} of every particle
    This function calcuates the kinetic energy
    outptu: Total kinetic energy of the sytem
    """

    velocities = np.sqrt(
        velocities["x"] ** 2 + velocities["y"] ** 2 + velocities["z"] ** 2
    )

    return 0.5 * np.sum(velocities**2)


def Lambda(velocities):
    """
    Rescalement factor to get the correct equilibrium
    """
    # print("target velocity: " + str(np.sqrt((3 * boltzmann * T) / (epsilon))))

    return np.sqrt(
        3 * (N - 1) * T / np.sum(velocities**2)
    )  # np.sqrt((3 * (N - 1) * boltzmann * T) / (np.sum(velocities**2) * m))


def PrepareLattice():
    """
    Function to set the lattice positions and velocities of the particles
    output: dictionary with particle positions {"x": ..., "y": ..., "z": ...},
            dictionary with particle velocities {"x": ..., "y": ..., "z": ...}
    """

    # Define lattice positions in 3D
    x1 = np.array([0, 1 / 3 * L, 2 / 3 * L] * 9)
    y1 = np.array(([0] * 3 + [1 / 3 * L] * 3 + [2 / 3 * L] * 3) * 3)
    z1 = np.repeat([0, 1 / 3 * L, 2 / 3 * L], 9)

    # Compute other lattice positions
    x2 = x1 + np.repeat(1 / 6 * L, 27)
    y2 = y1 + np.repeat(1 / 6 * L, 27)
    z2 = z1

    x3 = x1 + np.repeat(1 / 6 * L, 27)
    y3 = y1
    z3 = z1 + np.repeat(1 / 6 * L, 27)

    x4 = x1
    y4 = y1 + np.repeat(1 / 6 * L, 27)
    z4 = z1 + np.repeat(1 / 6 * L, 27)

    # Concatenate positions
    x = np.concatenate((x1, x2, x3, x4))
    y = np.concatenate((y1, y2, y3, y4))
    z = np.concatenate((z1, z2, z3, z4))

    positions = {"x": x, "y": y, "z": z}

    # Initialize velocities

    velocities = np.array(
        [random.normalvariate(0, 2 * T) for _ in range(N)]
    )  # TODO is dit goed??
    velocities = np.random.normal(0, T, N)
    velocities = Lambda(velocities) * velocities

    velocity_XYdirections = np.array([random.uniform(-np.pi, np.pi) for _ in range(N)])
    velocity_Zdirections = np.array([random.uniform(0, np.pi) for _ in range(N)])

    # Compute 3D velocities
    velocities = {
        "x": np.array(
            [
                np.cos(velocity_XYdirections[i])
                * np.cos(velocity_Zdirections[i])
                * velocities[i]
                for i in range(N)
            ]
        ),
        "y": np.array(
            [
                np.sin(velocity_XYdirections[i])
                * np.cos(velocity_Zdirections[i])
                * velocities[i]
                for i in range(N)
            ]
        ),
        "z": np.array(
            [np.sin(velocity_Zdirections[i]) * velocities[i] for i in range(N)]
        ),
    }

    # testing:
    # positions = {"y": np.array([2,2,2]),"x": np.array([1,1,1]), "z": np.array([1,3,2])}
    # velocities = {"y": np.array([2,2,2]),"x": np.array([0,0,0]), "z": np.array([1,-1,0])}
    return positions, velocities


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


def main():
    # prepare the initial state of the lattice
    lattice_positions = PrepareLattice()[0]
    lattice_velocities = PrepareLattice()[1]

    Argon_system = System_of_particles(
        positions=lattice_positions, velocities=lattice_velocities
    )

    positions_over_time = [Argon_system.positions]
    velocities_over_time = [Argon_system.velocities]

    potential_over_time = []
    kinetic_over_time = []
    """
    fig = plt.figure(figsize=(9, 9))
    if not z_dimension:
        ax = fig.add_subplot(projection="3d")
    else:
        ax = fig.add_subplot()
    # Set limits
    """

    # start the cycle of simulations
    global time
    for time in np.arange(0, t_end, h):
        Argon_system.Change_positions()

        if show_system:
            Argon_system.system_plotter(f"Time: {time}")

        # calculate the potential and kinetic energy at each time step

        # print(Lambda(velocities))

        """
        # quite programm if explotion occured 
        if np.sum(Calculate_potential(Argon_system.positions) + 
                  Calculate_kinetic(Argon_system.velocities)) > 10**10:
            print("Explotion")
            sys.exit()
        """

        positions = copy.deepcopy(Argon_system.positions)
        velocities = copy.deepcopy(Argon_system.velocities)
        positions_over_time.append(positions)
        velocities_over_time.append(velocities)

        potential_over_time.append(np.sum(Calculate_potential(Argon_system.positions)))
        kinetic_over_time.append(Calculate_kinetic(Argon_system.velocities))

        velocities = np.sqrt(
            Argon_system.velocities["x"] ** 2
            + Argon_system.velocities["y"] ** 2
            + Argon_system.velocities["z"] ** 2
        )

    total_energy = np.array(potential_over_time) + np.array(kinetic_over_time)

    plt.figure()
    plt.plot(
        np.arange(0, t_end, h),
        potential_over_time,
        label="Potential energy",
    )
    plt.plot(
        np.arange(0, t_end, h),
        kinetic_over_time,
        label="Kinetic energy",
    )
    # plot the equilibrium energy:
    E_target = (N - 1) * 3 / 2 * T
    plt.hlines(E_target, 0, t_end, label="equilibrium energy", color="purple")

    plt.plot(
        np.arange(0, t_end, h),
        total_energy,
        label="Total energy",
        linestyle="-",
    )
    plt.legend()
    plt.xlabel("time (t tilde)")
    plt.ylabel("energy (\u03B5)")
    plt.title(f"Energies during the simulation with rho = {rho} and T = {T}")
    plt.show()

    # save the data for further analysis
    if save_measurement:
        np.save(f"positions{T}_{rho}_long", np.array(positions_over_time))
        np.save(f"velocities{T}_{rho}_long", np.array(velocities_over_time))


if __name__ == "__main__":
    main()


# %%


# %% test

n = np.array([1, 5, -7, 2])
r = np.array([1, 2, 3, 4])
g = n / (r**2 * 0.1)
print(n % 3)