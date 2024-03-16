import math
import random
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.animation import FuncAnimation
from functools import partial
import numpy as np
import copy

# class Box:
#     def __init__(self, L, N):
#         self.L = L
#         self.N = N

'''
Natan verandert: N vervangen door len(positions["x"]) waar handig 
                 Functie descripties erbij gezet
                 def main() opgeschoond, lattice initialised met velocity in lattice functie
                 aparte functie inside the class gemaakt om te plotten
'''

%matplotlib inline


z_dimension = True
boltzmann = 1.38064852 * 10**-23
m = 39.948 * 1.66053906660 * 10**-27
epsilon = 119.8 * boltzmann
sigma = 3.405 * 10**-10

z_dimension = True
N = 108 if z_dimension else 18

T = 100
rho = 0.8

L = (N / rho) ** (1 / (3 if z_dimension else 2))
h = 0.001  # 10**-15

timesteps = 50
t_end = timesteps * h

save_measurement = True #if this is true write the data to an external file for further analysis 



class System_of_particles:
    def __init__(self, positions, velocities):
        self.positions = positions  # [[],[]],[],[]] # {"x": [2, 4, 13, 314 3, 234 2], "y": [[],[],[]]}
        self.velocities = velocities
        self.forces = {"x": [], "y": [], "z":[]}

    def Change_positions(self):
        previous_positions = {"x": np.zeros(N),"y": np.zeros(N),"z": np.zeros(N)}

        previous_positions["x"] = copy.deepcopy(self.positions["x"])
        previous_positions["y"] = copy.deepcopy(self.positions["y"])
        if z_dimension:
            previous_positions["z"] = copy.deepcopy(self.positions["z"])
        # self.previous_position = [self.particles[i]["position"] for i in range(len(self.particles))]

        current_force = Calculate_force(previous_positions)
        self.forces = current_force

        # current_force = [self.Calculate_force(self.particles[i]["position"]) for i in range(len(self.particles))]

        self.positions["x"] = (self.positions["x"] + self.velocities["x"] * h
            + current_force["x"] * h**2 / (2))
        self.positions["y"] = (self.positions["y"]+ self.velocities["y"] * h
            + current_force["y"] * h**2 / (2))
        if z_dimension:
            self.positions["z"] = (self.positions["z"] + self.velocities["z"] * h
                + current_force["z"] * h**2 / (2))
        # self.positions = [self.particles[i]["position"] + self.particles[i]["velocity"] * h + current_force[i] * h**2 / (2*m) for i in range(len(self.particles))]

        for i in range(N):
            self.positions["x"][i] = self.positions["x"][i] % L
            self.positions["y"][i] = self.positions["y"][i] % L
            if z_dimension:
                self.positions["z"][i] = self.positions["z"][i] % L

        # x(t) = self.previous_position
        force = Calculate_force(self.positions)

        self.velocities["x"] = self.velocities["x"] + (
            force["x"] + current_force["x"]
        ) * h / (2)
        self.velocities["y"] = self.velocities["y"] + (
            force["y"] + current_force["y"]
        ) * h / (2)
        if z_dimension:
            self.velocities["z"] = self.velocities["z"] + (
                force["z"] + current_force["z"]
            ) * h / (2)

        if time == h % 10 and time < (50 * h) and False:
            velocities = np.sqrt(self.velocities["x"] ** 2 + self.velocities["y"] ** 2)
            if z_dimension:
                velocities = np.sqrt(
                    self.velocities["x"] ** 2
                    + self.velocities["y"] ** 2
                    + self.velocities["z"] ** 2
                )
            lambda_ = Lambda(velocities)
            self.velocities["x"] = lambda_ * self.velocities["x"]
            self.velocities["y"] = lambda_ * self.velocities["y"]
            if z_dimension:
                self.velocities["z"] = lambda_ * self.velocities["z"]
        #print(self.velocities["y"])

        # return self
        
    def system_plotter(self, title):
        positions = self.positions
        velocities = self.velocities
        #forces = self.forces
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(positions['x'], positions['y'], positions['z'], c='b', marker='o', s=5, alpha=0.7) #, label='FCC lattice'
        plt.quiver(positions["x"], positions["y"],positions["z"],
                   velocities["x"], velocities["y"],velocities["z"], length = 0.1) 
        # Set labels and title for the plot
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('FCC Lattice')
        
        ax.set_xlim(0, L)
        ax.set_ylim(0, L)
        ax.set_zlim(0, L)
        
        plt.title(title)
    
        
        # Add legend and grid to the plot
        # ax.legend()
        ax.grid(True, linestyle='dotted', linewidth=0.5, alpha=0.5)
    
        # Display the plot
        plt.show()            


def Calculate_force(positions):
    
    '''
    This function calculates the force that each particle experiences
    
    input: dictionary with x,y,z positions {x: ... y: ... z:...}
    output: dictionary with the forces in the x,y,z direction of 
            each particle {x: ... y: ... z:...}
    
    '''
    
    force_dict = {"x": np.zeros(N),"y": np.zeros(N),"z": np.zeros(N)}

    for i in range(N):
        # cell_particle_i
        for j in range(N):
            # cell_particle_j
            if i != j:  # and cell_particle_i = cell_particle_j:
                delta_x = (positions["x"][j] - positions["x"][i] + L / 2) % L - L / 2
                delta_y = (positions["y"][j] - positions["y"][i] + L / 2) % L - L / 2
                r = ((delta_x) ** 2 + (delta_y) ** 2) ** 0.5
                if z_dimension:
                    delta_z = (positions["z"][j] - positions["z"][i] + L / 2) % L - L / 2
                    r = ((r) ** 2 + (delta_z) ** 2) ** 0.5
                    Zangle = math.atan2(delta_z, r)
                if r == 0:
                    print("Zero")
                XYangle = math.atan2(delta_y, delta_x)

                force = (-24 * (2 * (1 / r) ** 13 - (1 / r) ** 7))

                force_dict["x"][i] += force * math.cos(XYangle)
                force_dict["y"][i] += force * math.sin(XYangle)
                if z_dimension:
                    force_dict["z"][i] += force * math.sin(Zangle)

    return force_dict


def Calculate_potential(positions):
    '''
    input: Dictionary with x,y,z positions {x: ... y: ... z:...}
    This function calcuates the potential between all the particles 
    outptu: Array with the potential energy of every particle  
    '''

    potential_list = []
    
    N = len(positions["x"]) 

    for i in range(N):
        # cell_particle_i
        potential = 0

        for j in range(i + 1, N):
            # cell_particle_j

            delta_x = (positions["x"][j] - positions["x"][i] + L / 2) % L - L / 2
            delta_y = (positions["y"][j] - positions["y"][i] + L / 2) % L - L / 2
            r = ((delta_x) ** 2 + (delta_y) ** 2) ** 0.5
            if z_dimension:
                delta_z = (positions["z"][j] - positions["z"][i] + L / 2) % L - L / 2
                r = ((r) ** 2 + (delta_z) ** 2) ** 0.5
            if np.isclose(r, 0):
                print("Zero")

            potential += 4 * ((1 / r) ** 12 - (1 / r) ** 6)

        potential_list.append(potential)

    return potential_list


def Calculate_kinetic(velocities):
    '''
    input: Dictionary with x,y,z velocities {x: ... y: ... z:...} of every particle
    This function calcuates the kinetic energy
    outptu: Total kinetic energy of the sytem
    '''
    
    if z_dimension:
        velocities = np.sqrt(
            velocities["x"] ** 2 + velocities["y"] ** 2 + velocities["z"] ** 2
        )
    else:
        velocities = np.sqrt(velocities["x"] ** 2 + velocities["y"] ** 2)

    return 0.5 * np.sum(velocities**2)


def Lambda(velocities):
    '''
    Rescalement factor to get the correct equilibrium
    '''
    # print("target velocity: " + str(np.sqrt((3 * boltzmann * T) / (epsilon))))

    return np.sqrt(
        3 * (N - 1) * T / np.sum(velocities**2)
    )  # np.sqrt((3 * (N - 1) * boltzmann * T) / (np.sum(velocities**2) * m))


def PrepareLattice():
    '''
    Function to set the lattice positions and velocities of the particles
    output: dictionary with particle positions {"x": ..., "y": ..., "z": ...}, 
            dictionary with particle velocities {"x": ..., "y": ..., "z": ...}
    '''

    z_dimension = True  # Corrected variable name

    if z_dimension:
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
    else:
        # Define lattice positions in 2D
        x1 = np.array([0, 1 / 3 * L, 2 / 3 * L] * 3)
        y1 = np.array([0] * 3 + [1 / 3 * L] * 3 + [2 / 3 * L] * 3)

        # Compute other lattice positions
        x2 = x1 + np.repeat(1 / 6 * L, 9)
        y2 = y1 + np.repeat(1 / 6 * L, 9)

        # Concatenate positions
        x = np.concatenate((x1, x2))
        y = np.concatenate((y1, y2))
        z = None

    # Initialize velocities

    velocities = np.array([random.normalvariate(0, 2 * T) for _ in range(N)])
    velocities = Lambda(velocities) * velocities

    velocity_XYdirections = np.array([random.uniform(-np.pi, np.pi) for _ in range(N)])
    velocity_Zdirections = np.array([random.uniform(0, np.pi) for _ in range(N)])
    

    if z_dimension:
        # Compute 3D velocities
        velocities = {
            "x": np.array([np.cos(velocity_XYdirections[i]) * np.cos(velocity_Zdirections[i]) * velocities[i] for i in range(N)]),
            "y": np.array([np.sin(velocity_XYdirections[i]) * np.cos(velocity_Zdirections[i]) * velocities[i] for i in range(N)]),
            "z": np.array([np.sin(velocity_Zdirections[i]) * velocities[i] for i in range(N)])
        }
    else:
        # Compute 2D velocities
        velocities = {
            "x": np.array([np.cos(velocity_XYdirections[i]) * velocities[i] for i in range(N)]),
            "y": np.array([np.sin(velocity_XYdirections[i]) * velocities[i] for i in range(N)])
        }

    return {"x": x, "y": y, "z": z}, velocities

def pair_correlation(system):
    '''
    input: position dictionary
    '''
    positions = system.positions
    
    particle_distances = np.array([])
    
    for i in range(N):
        for j in range(i + 1, N):
            delta_x = (positions["x"][j] - positions["x"][i] + L / 2) % L - L / 2
            delta_y = (positions["y"][j] - positions["y"][i] + L / 2) % L - L / 2
            r = ((delta_x) ** 2 + (delta_y) ** 2) ** 0.5
            if z_dimension:
                delta_z = (positions["z"][j] - positions["z"][i] + L / 2) % L - L / 2
                r = ((r) ** 2 + (delta_z) ** 2) ** 0.5
            if np.isclose(r, 0):
                print("Zero")
                
            particle_distances = np.append(particle_distances,r)
            
            
            
    y = len(particle_distances)
    #print(particle_distances)
    
    #plt.scatter(np.arange(y),particle_distances)
    
    bin_size = 0.1
    plt.hist(particle_distances,np.arange(0,L,bin_size))
     
    
    plt.show()
    
def pressure(system,rho):
    positions  = system.positions
    
    potentials = np.array([])
    
    for i in range(N):
        for j in range(i + 1, N):
            delta_x = (positions["x"][j] - positions["x"][i] + L / 2) % L - L / 2
            delta_y = (positions["y"][j] - positions["y"][i] + L / 2) % L - L / 2
            r = ((delta_x) ** 2 + (delta_y) ** 2) ** 0.5
            if z_dimension:
                delta_z = (positions["z"][j] - positions["z"][i] + L / 2) % L - L / 2
                r = ((r) ** 2 + (delta_z) ** 2) ** 0.5
            if np.isclose(r, 0):
                print("Zero")
            
            potential_der = 2 * (-24*(1 / r) ** 13 + 12* (1 / r) ** 7)
            #print(potential_der)
            
            potentials = np.append(potentials,r*potential_der)
    expec_value = np.mean(potentials)
    pressure = epsilon*rho*(1 - 1/(3*len(positions["x"])*epsilon)*expec_value) #TODO something wrong with dimensionless units?
    
    return pressure



def main():

    #prepare the initial state of the lattice
    lattice_positions = PrepareLattice()[0]
    lattice_velocities = PrepareLattice()[1]

    Argon_system = System_of_particles(
        positions=lattice_positions,velocities = lattice_velocities)

    positions_over_time = [Argon_system.positions]
    velocities_over_time = [Argon_system.velocities]

    potential_over_time = []
    kinetic_over_time = []

    fig = plt.figure(figsize=(9, 9))
    if not z_dimension:
        ax = fig.add_subplot(projection="3d")
    else:
        ax = fig.add_subplot()
    # Set limits

    #start the cycle of simulations
    global time
    for time in np.arange(0, t_end, h):
        Argon_system.Change_positions()
        
        Argon_system.system_plotter(f"Time: {time}")

        #calculate the potential and kinetic energy at each time step
        potential_over_time.append(np.sum(Calculate_potential(Argon_system.positions)))
        kinetic_over_time.append(Calculate_kinetic(Argon_system.velocities))

        positions = copy.deepcopy(Argon_system.positions)
        velocities = copy.deepcopy(Argon_system.velocities)
        positions_over_time.append(positions)
        velocities_over_time.append(velocities)
        
        
        
        
        '''
        ax.cla()
        ax.set_title(f"Time: {time}")
        ax.set_xlim(0, L)
        ax.set_ylim(0, L)
        if not z_dimension:
            plt.title('why?')
            ax.set_zlim(0, L)
            ax.quiver(
                Argon_system.positions["x"],
                Argon_system.positions["y"],
                Argon_system.positions["z"],
                Argon_system.velocities["x"],
                Argon_system.velocities["y"],
                Argon_system.velocities["z"],
                length=0.1,
            )
        else:
            plt.scatter(
                Argon_system.positions["x"],
                Argon_system.positions["y"], s=5
                # color=["red", "blue"],
            )
            ax.quiver(
                Argon_system.positions["x"],
                Argon_system.positions["y"],
                Argon_system.forces["x"],
                Argon_system.forces["y"],
            )
        plt.pause(0.000005)
    
    plt.show()
    
    plt.figure()
    plt.xlim(0, L)
    plt.ylim(0, L)
    for i in range(N):
        plt.plot(
            [positions_over_time[j]["x"][i] for j in range(len(positions_over_time))],
            [positions_over_time[j]["y"][i] for j in range(len(positions_over_time))],
            color="blue",
            alpha=0.5,
        )
        plt.title("what is this?")
    plt.show()
    '''

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
    plt.plot(
        np.arange(0, t_end, h),
        np.array(potential_over_time) + np.array(kinetic_over_time),
        label="Total energy",
        linestyle="--",
    )
    plt.legend()
    plt.show()
    
    # save the data for further analysis
    if save_measurement:
        np.save('positions',np.array(positions_over_time))
        np.save('velocities',np.array(velocities_over_time))
    

    # for time in np.arange(h, t_end, h):
    #     print(time)
    #     particles = copy.deepcopy(
    #         positions_over_time[
    #             round(time - h, 10)  # round to avoid floating point errors
    #         ]
    #     )  # Copy the argons from the previous time step.

    #     for argon in particles:
    #         argon.F = {"x": 0, "y": 0}

    #         for argon2 in particles:
    #             if argon != argon2:
    #                 argon.F["x"] += argon.Potential(argon2)["x"]
    #                 argon.F["y"] += argon.Potential(argon2)["y"]

    #     for argon in particles:
    #         argon.Evolve()

    #     positions_over_time[
    #         round(time, 10)
    #     ] = particles  # round to avoid floating point errors

    # for time in np.arange(dt, t_end, dt):
    #     time = round(time, 10)
    #     plt.clf()
    #     plt.xlim(0, L)
    #     plt.ylim(0, L)
    #     plt.quiver(
    #         [[positions_over_time[time][i_argon].position["x"] for i_argon in range(N)],
    #         [positions_over_time[time][i_argon].position["y"] for i_argon in range(N)]],
    #         [positions_over_time[time][i_argon].velocity["x"] for i_argon in range(N)],
    #         [positions_over_time[time][i_argon].velocity["y"] for i_argon in range(N)],
    #     )
    #     plt.scatter(
    #         [
    #             positions_over_time[round(time - dt, 10)][i_argon].position["x"]
    #             for i_argon in range(N)
    #         ],
    #         [
    #             positions_over_time[round(time - dt, 10)][i_argon].position["y"]
    #             for i_argon in range(N)
    #         ],
    #         alpha=0.5,
    #     )
    #     plt.pause(0.1)


if __name__ == "__main__":
    main()
    
    
#%%

import numpy as np
'''
In this part the data is analysed 
'''


positions = np.load('positions.npy', allow_pickle=True)
velocities = np.load('velocities.npy', allow_pickle=True)


def pair_correlation(positions):
    '''
    input: position dictionary
    '''
    
    particle_distances = np.array([])
    
    for i in range(N):
        for j in range(i + 1, N):
            delta_x = (positions["x"][j] - positions["x"][i] + L / 2) % L - L / 2
            delta_y = (positions["y"][j] - positions["y"][i] + L / 2) % L - L / 2
            delta_z = (positions["z"][j] - positions["z"][i] + L / 2) % L - L / 2
            r = ((delta_x) ** 2 + (delta_y)**2 + (delta_z) ** 2) ** 0.5

                
            particle_distances = np.append(particle_distances,r)
            
            
            
    y = len(particle_distances)
    #print(particle_distances)
    
    #plt.scatter(np.arange(y),particle_distances)
    
    bin_size = 0.01
    plt.hist(particle_distances,np.arange(0,L,bin_size))
     
    
    plt.show()
    
def pressure(positions,rho):
    
    interaction_term = np.array([])
    
    for i in range(N):
        for j in range(i + 1, N):
            delta_x = (positions["x"][j] - positions["x"][i] + L / 2) % L - L / 2
            delta_y = (positions["y"][j] - positions["y"][i] + L / 2) % L - L / 2
            delta_z = (positions["z"][j] - positions["z"][i] + L / 2) % L - L / 2
            
            r = ((delta_x) ** 2 + (delta_y) ** 2 + (delta_z) ** 2) ** 0.5
        
            potential_der = 2 * (-24*(1 / r) ** 13 + 12* (1 / r) ** 7)
            #print(potential_der)
            
            interaction_term = np.append(interaction_term,r*potential_der)
    expec_value = np.mean(interaction_term)
    pressure = epsilon*rho*(1 - 1/(3*len(positions["x"])*epsilon)*expec_value) #TODO something wrong with dimensionless units?
    
    return pressure
    


def main():
    #pair-correlation
    hist_plots = np.arange(len(positions))
    pressures = np.array([])
    for i in hist_plots:
        pair_correlation(positions[i])
        print('pressure is:',i,pressure(positions[i], 0.8))
        pressures = np.append(pressures,pressure(positions[i],0.8))
        
    plt.plot(pressures)
    
    # pressure
    #  for i in hist_plots():
    #    print(pressure(positions[i], 0.8))
    
    
    

if __name__ == "__main__":
    main()

