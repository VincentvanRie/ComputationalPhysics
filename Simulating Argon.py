"""
Molecular Dynamics Simulation of Argon

This script simulates the molecular dynamics of Argon atoms in a 3D box.
The atoms are placed in a face-centered cubic (FCC) lattice and are given random velocities.
The system is evolved in time using the Verlet algorithm. The system is periodic in all directions.
The Lennard-Jones potential is used to calculate the forces between the particles.

The script calculates the potential, kinetic, and total energy of the system at each time step.


Authors:
- Natan van Steenis
- Vincent van Rie

Date: 22/03/2024

"""

import random
import matplotlib.pyplot as plt
import numpy as np
import copy

# %matplotlib inline
"""
Set parameters for the simulation
"""

"""Main parameters for the simulation"""
T = 1.0
rho = 0.8

"""Additional parameters for the simulation"""
N = 1372  # number of particles for big system
N = 108  # number of particles for small system

# parameter to determine how many radial cells to check for particles.
# Where if equal to sim_size, all particles are checked
perimeter_parameter = 2

# number of iterations
timesteps = 100

"""Boolean parameters for the simulation"""
show_system = True  # plot the system at each iteration
save_measurement = (
    False  # if this is true write the data to an external file for further analysis
)

"""
End set parameters of the simulation
"""
# timestep
h = 0.01

t_end = timesteps * h

sim_size = round((N / 4) ** (1 / 3))
L = (N / rho) ** (1 / 3)


# set the values for lambda renormalization
frequency_rescale = 20
number_of_rescales = 5

renorm_times = np.arange(
    0, t_end, frequency_rescale * h
)  # at these times lambda renormalization is performed
renorm_times = renorm_times[:number_of_rescales]


class System_of_particles:
    def __init__(self, positions, velocities):
        self.positions = positions
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

        current_force = Calculate_force(previous_positions)
        self.forces = current_force

        for key in dimensions:
            self.positions[key] = (
                self.positions[key]
                + self.velocities[key] * h
                + current_force[key] * h**2 / (2)
            )
            self.positions[key] = self.positions[key] % L

        force = Calculate_force(self.positions)

        for key in dimensions:
            self.velocities[key] = self.velocities[key] + (
                force[key] + current_force[key]
            ) * h / (2)

        print(time)
        if any(np.isclose(time, renorm_time) for renorm_time in renorm_times):
            print("yeah i rescaled at time", time)
            velocities = np.sqrt(
                self.velocities["x"] ** 2
                + self.velocities["y"] ** 2
                + self.velocities["z"] ** 2
            )

            lambda_ = Lambda(velocities)
            for key in dimensions:
                self.velocities[key] = lambda_ * self.velocities[key]

    def system_plotter(self, title):
        """
        Function to plot the current system. Input is the title for the plot
        and the output is a plotted figure.
        """
        positions = self.positions
        velocities = self.velocities

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

    First the function creates a grid of cells and assigns each particle to a cell.
    and assigns each cell a list of particles that are in the cell.

    input: dictionary with x,y,z positions {x: ... y: ... z:...}
    output: dictionary with the forces in the x,y,z direction of
            each particle {x: ... y: ... z:...}

    """

    # Initialize the size of cells to be about 1 sigma
    number_of_cells = round(L + 0.5)
    cells_with_particles = np.full(
        (number_of_cells, number_of_cells, number_of_cells), False, dtype=object
    )

    particles_with_cell_index = [[] for _ in range(N)]

    for i in range(N):
        x = int(positions["x"][i] // (L / number_of_cells))
        y = int(positions["y"][i] // (L / number_of_cells))
        z = int(positions["z"][i] // (L / number_of_cells))

        # Save the cell coordinates of the particle
        particles_with_cell_index[i] = [x, y, z]

        # And save the particle index in the cell
        if not cells_with_particles[
            x, y, z
        ]:  # if the cell is empty add the particle index to the cell
            cells_with_particles[x, y, z] = [i]
        else:  # if the cell is not empty append the particle index to the list
            cells_with_particles[x, y, z].append(i)

    force_dict = {"x": np.zeros(N), "y": np.zeros(N), "z": np.zeros(N)}

    for i in range(N):  # i = particle 1, j = particle 2
        x_cell_i, y_cell_i, z_cell_i = particles_with_cell_index[i]

        perimeter_parameter = 1

        # Retrieve all particles in the perimeter of the cell of particle i
        nearby_particles = []
        for n in range(-perimeter_parameter, perimeter_parameter + 1):
            for j in range(-perimeter_parameter, perimeter_parameter + 1):
                for k in range(-perimeter_parameter, perimeter_parameter + 1):
                    if cells_with_particles[
                        (x_cell_i + n) % number_of_cells,
                        (y_cell_i + j) % number_of_cells,
                        (z_cell_i + k) % number_of_cells,
                    ]:
                        for particle in cells_with_particles[
                            (x_cell_i + n) % number_of_cells,
                            (y_cell_i + j) % number_of_cells,
                            (z_cell_i + k) % number_of_cells,
                        ]:
                            nearby_particles.append(particle)

        for j in nearby_particles:
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
    This function calculates the potential energy of given positions

    First the function creates a grid of cells and assigns each particle to a cell.
    and assigns each cell a list of particles that are in the cell.

    input: Dictionary with x,y,z positions {x: ... y: ... z:...}
    This function calcuates the potential between all the particles
    output: Array with the potential energy of every particle
    """

    potential_list = []

    N = len(positions["x"])

    # Initialize the size of cells to be about 1 sigma
    number_of_cells = round(L + 0.5)
    cells_with_particles = np.full(
        (number_of_cells, number_of_cells, number_of_cells), False, dtype=object
    )

    particles_with_cell_index = [[] for _ in range(N)]

    for i in range(N):
        x = int(positions["x"][i] // (L / number_of_cells))
        y = int(positions["y"][i] // (L / number_of_cells))
        z = int(positions["z"][i] // (L / number_of_cells))

        # Save the cell coordinates of the particle
        particles_with_cell_index[i] = [x, y, z]

        # And save the particle index in the cell
        if not cells_with_particles[
            x, y, z
        ]:  # if the cell is empty add the particle index to the cell
            cells_with_particles[x, y, z] = [i]
        else:  # if the cell is not empty append the particle index to the list
            cells_with_particles[x, y, z].append(i)

    for i in range(N):  # i = particle 1, j = particle 2
        x_cell_i, y_cell_i, z_cell_i = particles_with_cell_index[i]

        # Retrieve all particles in the perimeter of the cell of particle i
        nearby_particles = []
        for n in range(-perimeter_parameter, perimeter_parameter + 1):
            for j in range(-perimeter_parameter, perimeter_parameter + 1):
                for k in range(-perimeter_parameter, perimeter_parameter + 1):
                    if cells_with_particles[
                        (x_cell_i + n) % number_of_cells,
                        (y_cell_i + j) % number_of_cells,
                        (z_cell_i + k) % number_of_cells,
                    ]:
                        for particle in cells_with_particles[
                            (x_cell_i + n) % number_of_cells,
                            (y_cell_i + j) % number_of_cells,
                            (z_cell_i + k) % number_of_cells,
                        ]:
                            if particle > i:
                                nearby_particles.append(particle)

        potential = 0

        for j in nearby_particles:
            r = two_part_distance(i, j, positions)

            if np.isclose(r, 0):
                print("Zero")

            potential += 4 * ((1 / r) ** 12 - (1 / r) ** 6)

        potential_list.append(potential)

    return potential_list


def Calculate_kinetic(velocities):
    """
    This function calculates the kinetic energy given the velocities.

    input: Dictionary with x,y,z velocities {x: ... y: ... z:...} of every particle
    This function calcuates the kinetic energy
    output: Total kinetic energy of the sytem
    """

    velocities = np.sqrt(
        velocities["x"] ** 2 + velocities["y"] ** 2 + velocities["z"] ** 2
    )

    return 0.5 * np.sum(velocities**2)


def Lambda(velocities):
    """
    Rescalement factor to get the correct equilibrium
    """

    return np.sqrt(3 * (N - 1) * T / np.sum(velocities**2))


def PrepareLatticeBlock():
    """
    Function to set the lattice positions of the particles in a unit block

    output: dictionary with particle positions {"x": ..., "y": ..., "z": ...}
    """
    x = np.array([0, 0, 1 / 2, 1 / 2])
    y = np.array([0, 1 / 2, 0, 1 / 2])
    z = np.array([0, 1 / 2, 1 / 2, 0])

    positions = {"x": x, "y": y, "z": z}

    return positions


def PrepareLattice():
    """
    Function to set the lattice positions and velocities of the particles

    output: dictionary with particle positions {"x": ..., "y": ..., "z": ...},
            dictionary with particle velocities {"x": ..., "y": ..., "z": ...}
    """

    positions = {"x": [], "y": [], "z": []}

    for k in range(sim_size):
        for j in range(sim_size):
            for i in range(sim_size):
                positions["x"] = np.concatenate(
                    (
                        positions["x"],
                        PrepareLatticeBlock()["x"] + i,
                    )
                )
                positions["y"] = np.concatenate(
                    (
                        positions["y"],
                        PrepareLatticeBlock()["y"] + j,
                    )
                )
                positions["z"] = np.concatenate(
                    (
                        positions["z"],
                        PrepareLatticeBlock()["z"] + k,
                    )
                )

    for key in positions:
        positions[key] = positions[key] * L / sim_size

    # Initialize velocities
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

    return positions, velocities


def two_part_distance(p1, p2, positions):
    """
    Take the position dictionary and 2 particles
    as input and return the distance between p1 and p1

    input: position dictionary and 2 particles
    output: distance between p1 and p2
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

    # start the cycle of simulations
    # Time is global so that the lambda renormalization can be done at the correct times.
    # We could give it as a parameter to the function, but this is easier.
    global time
    for time in np.arange(0, t_end, h):
        Argon_system.Change_positions()

        if show_system:
            Argon_system.system_plotter(f"Time: {time}")

        # calculate the potential and kinetic energy at each time step
        positions = copy.deepcopy(Argon_system.positions)
        velocities = copy.deepcopy(Argon_system.velocities)
        positions_over_time.append(positions)
        velocities_over_time.append(velocities)

        potential_over_time.append(np.sum(Calculate_potential(Argon_system.positions)))
        kinetic_over_time.append(Calculate_kinetic(Argon_system.velocities))

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
    plt.ylabel("Energy (\u03B5)")
    plt.title(f"Energies during the simulation with rho = {rho} and T = {T}")
    plt.xlabel(r"Time [$\left(m \sigma^2/\epsilon\right)^{-\frac{1}{2}}]$")
    plt.legend()
    plt.savefig("solid_energies.pdf", dpi=1200)
    plt.show()

    # save the data for further analysis
    if save_measurement:
        np.save(f"positions{T}_{rho}_long", np.array(positions_over_time))
        np.save(f"velocities{T}_{rho}_long", np.array(velocities_over_time))


if __name__ == "__main__":
    main()
