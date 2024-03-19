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

'''
Natan verandert: N vervangen door len(positions["x"]) waar handig 
                 Functie descripties erbij gezet
                 def main() opgeschoond, lattice initialised met velocity in lattice functie
                 aparte functie inside the class gemaakt om te plotten
                 2D weggehaald
                 
'''

%matplotlib inline


z_dimension = True
boltzmann = 1.38064852 * 10**-23
m = 39.948 * 1.66053906660 * 10**-27
epsilon = 119.8 * boltzmann
sigma = 3.405 * 10**-10

N = 108 #number of particles

T = 0.5
rho = 1.2

L = (N / rho) ** (1 / 3)
h = 0.01  # 10**-15

timesteps = 300 #nujmber of iterations
t_end = timesteps * h

save_measurement = True #if this is true write the data to an external file for further analysis 



class System_of_particles:
    def __init__(self, positions, velocities):
        self.positions = positions  # [[],[]],[],[]] # {"x": [2, 4, 13, 314 3, 234 2], "y": [[],[],[]]}
        self.velocities = velocities
        self.forces = {"x": [], "y": [], "z":[]}

    def Change_positions(self):
        '''
        Function that updates the system to the next timestep. 
        '''
        previous_positions = {"x": np.zeros(N),"y": np.zeros(N),"z": np.zeros(N)}
        dimensions = ["x", "y", "z"]
        
        
        for key in dimensions:
            previous_positions[key] = copy.deepcopy(self.positions[key])
            # self.previous_position = [self.particles[i]["position"] for i in range(len(self.particles))]

        current_force = Calculate_force(previous_positions)
        self.forces = current_force

        # current_force = [self.Calculate_force(self.particles[i]["position"]) for i in range(len(self.particles))]

        for key in dimensions: 
            self.positions[key] = (self.positions[key] +    self.velocities[key] * h
                + current_force[key] * h**2 / (2)) 
            self.positions[key] = self.positions[key] % L
        
        # self.positions = [self.particles[i]["position"] + self.particles[i]["velocity"] * h + current_force[i] * h**2 / (2*m) for i in range(len(self.particles))]


        # x(t) = self.previous_position
        force = Calculate_force(self.positions)
        
        for key in dimensions:
            self.velocities[key] = self.velocities[key] + (force[key] + current_force[key]) * h / (2)


        if time == h % 10 and time < (50 * h):
            velocities = np.sqrt(self.velocities["x"] ** 2 + self.velocities["y"] ** 2
                                 + self.velocities["z"]**2)

            lambda_ = Lambda(velocities)
            for key in dimensions:
                self.velocities[key] = lambda_ * self.velocities[key]
            #print(self.velocities["y"])

        # return self
        
    def system_plotter(self, title):
        '''
        Function to plot the current system. Input is the title for the plot
        and the output is a plotted figure. 
        '''
        positions = self.positions
        velocities = self.velocities
        #forces = self.forces
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(positions['x'], positions['y'], positions['z'], c='b', marker='o', s=5, alpha=0.7) #, label='FCC lattice'
        plt.quiver(positions["x"], positions["y"],positions["z"],
                   velocities["x"], velocities["y"],velocities["z"], length = 0.5, normalize = True,color = 'red') 
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
        plt.pause(0.05)
        plt.close()

def Calculate_force(positions):
    
    '''
    This function calculates the force that each particle experiences
    
    input: dictionary with x,y,z positions {x: ... y: ... z:...}
    output: dictionary with the forces in the x,y,z direction of 
            each particle {x: ... y: ... z:...}
    
    '''
    
    force_dict = {"x": np.zeros(N),"y": np.zeros(N),"z": np.zeros(N)}

    for i in range(N): #i = p1, j = p2
        # cell_particle_i
        for j in range(N):
            # cell_particle_j
            if i != j:  # and cell_particle_i = cell_particle_j:
                delta_x = (positions["x"][j] - positions["x"][i] + L / 2) % L - L / 2
                delta_y = (positions["y"][j] - positions["y"][i] + L / 2) % L - L / 2
                delta_z = (positions["z"][j] - positions["z"][i] + L / 2) % L - L / 2
                
                r = ((delta_x) ** 2 + (delta_y) ** 2 + (delta_z) ** 2) ** 0.5
                
                Zangle = math.atan2(delta_z, r)
                if r == 0:
                    print("Zero")
                
                XYangle = math.atan2(delta_y, delta_x)

                force = (-24 * (2 * (1 / r) ** 13 - (1 / r) ** 7))

                force_dict["x"][i] += force * math.cos(XYangle)
                force_dict["y"][i] += force * math.sin(XYangle)
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
            r = two_part_distance(i, j, positions)
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
    
    positions = {"x": x,"y": y,"z": z}

    # Initialize velocities

    velocities = np.array([random.normalvariate(0, 2 * T) for _ in range(N)])
    velocities = Lambda(velocities) * velocities

    velocity_XYdirections = np.array([random.uniform(-np.pi, np.pi) for _ in range(N)])
    velocity_Zdirections = np.array([random.uniform(0, np.pi) for _ in range(N)])
    
    # Compute 3D velocities
    velocities = {
        "x": np.array([np.cos(velocity_XYdirections[i]) * np.cos(velocity_Zdirections[i]) * velocities[i] for i in range(N)]),
        "y": np.array([np.sin(velocity_XYdirections[i]) * np.cos(velocity_Zdirections[i]) * velocities[i] for i in range(N)]),
        "z": np.array([np.sin(velocity_Zdirections[i]) * velocities[i] for i in range(N)])}
    return positions, velocities

def pair_correlation(system):
    '''
    input: position dictionary
    '''
    positions = system.positions
    
    particle_distances = np.array([])
    
    for i in range(N):
        for j in range(i + 1, N):
            r = two_part_distance(i, j, positions)
                
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
            r = two_part_distance(i, j, positions)
            
            potential_der = 2 * (-24*(1 / r) ** 13 + 12* (1 / r) ** 7)
            #print(potential_der)
            
            potentials = np.append(potentials,r*potential_der)
    expec_value = np.mean(potentials)
    pressure = epsilon*rho*(1 - 1/(3*len(positions["x"])*epsilon)*expec_value) #TODO something wrong with dimensionless units?
    
    return pressure


def two_part_distance(p1,p2,positions):
    '''
    Take the position dictionary and 2 particles 
    as input and return the distance between p1 and p1
    
    '''

    delta_x = (positions["x"][p2] - positions["x"][p1] + L / 2) % L - L / 2
    delta_y = (positions["y"][p2] - positions["y"][p1] + L / 2) % L - L / 2
    delta_z = (positions["z"][p2] - positions["z"][p1] + L / 2) % L - L / 2
    r = ((delta_x) ** 2 + (delta_y) ** 2 + (delta_z) ** 2) ** 0.5
    
    return r


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
    '''
    fig = plt.figure(figsize=(9, 9))
    if not z_dimension:
        ax = fig.add_subplot(projection="3d")
    else:
        ax = fig.add_subplot()
    # Set limits
    '''

    #start the cycle of simulations
    global time
    for time in np.arange(0, t_end, h):
        Argon_system.Change_positions()
        
        Argon_system.system_plotter(f"Time: {time}")

        #calculate the potential and kinetic energy at each time step
        potential_over_time.append(np.sum(Calculate_potential(Argon_system.positions)))
        kinetic_over_time.append(Calculate_kinetic(Argon_system.velocities))

        
        '''
        # quite programm if explotion occured 
        if np.sum(Calculate_potential(Argon_system.positions) + 
                  Calculate_kinetic(Argon_system.velocities)) > 10**10:
            print("Explotion")
            sys.exit()
        '''    

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
        np.save(f'positions{T}_{rho}',np.array(positions_over_time))
        np.save(f'velocities{T}_{rho}',np.array(velocities_over_time))
    

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
import matplotlib.pyplot as plt
'''
In this part the data is analysed 
'''
## Specify for which rho/T the analysis is done
T = 0.5
rho = 1.2


positions = np.load(f'positions{T}_{rho}.npy', allow_pickle=True)
velocities = np.load(f'velocities{T}_{rho}.npy', allow_pickle=True)



N = len(positions[0]["x"])
L = (N / rho) ** (1 / (3))
print(L)
h = 0.001  # 10**-15
boltzmann = 1.38064852 * 10**-23
epsilon = 119.8 * boltzmann
bin_size = 0.1






def two_part_distance(p1,p2,positions):
    '''
    Take the position dictionary and 2 particles 
    as input and return the distance between p1 and p1
    
    '''

    delta_x = (positions["x"][p2] - positions["x"][p1] + L / 2) % L - L / 2
    delta_y = (positions["y"][p2] - positions["y"][p1] + L / 2) % L - L / 2
    delta_z = (positions["z"][p2] - positions["z"][p1] + L / 2) % L - L / 2
    r = ((delta_x) ** 2 + (delta_y) ** 2 + (delta_z) ** 2) ** 0.5
    
    return r


def pair_correlation_plot(positions):
    '''
    input: position dictionary
    output: plot if the histogram
    '''
    
    particle_distances = np.array([])
    
    for i in range(N):
        for j in range(i + 1, N):
            r = two_part_distance(i, j, positions)
            particle_distances = np.append(particle_distances,r)
            
            
            
    y = len(particle_distances)
    #print(particle_distances)
    
    #plt.scatter(np.arange(y),particle_distances)
    
    
    bins = L/bin_size
    
    plt.hist(particle_distances,np.arange(0,L,bin_size))
    plt.show()
    

    
def average_hist(all_positions,analyse_times,bin_size):
    '''
    input: all positions and the times to be analysed
    output: pair correlation histogram of the average hist,bins
    '''
    plot = True  #plots the avg_hist if true
    num_bins = int(L/bin_size)
    total_hist = np.zeros(num_bins)
    
    
    
    for t in analyse_times:
        positions = all_positions[t]
        particle_distances = np.array([])
        for i in range(N):
            for j in range(i + 1, N):
                r = two_part_distance(i, j, positions)
    
                    
                particle_distances = np.append(particle_distances,r)
        
        hist,bins = np.histogram(particle_distances,bins = num_bins)
        total_hist += hist
        
    
    avg_hist = total_hist/len(analyse_times)
    
    #print(bins[:-1])
    pair_correlation_g = (2*L**3)/(N*(N-1))*avg_hist/(4*np.pi*bins[:-1]**2*bin_size)
    
    if plot:
        plt.bar(bins[:-1], pair_correlation_g, width=np.diff(bins), edgecolor='black')
        plt.xlabel('Distance')
        plt.ylabel('Frequency')
        plt.title(f'Average Pair Correlation Histogram with rho = {rho},T = {T} ')
        plt.grid(True)
        plt.show()
    
    return pair_correlation_g, bins
    
    
def pressure(positions,rho):
    
    interaction_term = np.array([])
    
    for i in range(N):
        for j in range(i + 1, N):
            r = two_part_distance(i, j, positions)
        
            potential_der = 2 * (-24*(1 / r) ** 13 + 12* (1 / r) ** 7)
            #print(potential_der)
            
            interaction_term = np.append(interaction_term,r*potential_der)
    expec_value = np.mean(interaction_term)
    pressure = epsilon*rho*(1 - 1/(3*len(positions["x"])*epsilon)*expec_value) #TODO something wrong with dimensionless units?
    
    return pressure

def plot_fcc_lattice(positions):

    # Plot the FCC lattice
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    
    
    colors = np.random.rand(10, 3)
    
    ''' keep this part chatgpt
    for i in random_indices:
    
        # Scatter plot for FCC lattice
        ax.scatter(positions['x'][i], positions['y'][i], positions['z'][i], c='b', marker='o', s=5, alpha=0.7) #, label='FCC lattice'
    
        plt.quiver(system.positions["x"][i], system.positions["y"][i],system.positions["z"][i],
                        system.velocities["x"][i], system.velocities["y"][i],system.velocities["z"][i], length = 0.1)
        '''
    ax.scatter(positions['x'], positions['y'], positions['z'], c='b', marker='o', s=5, alpha=0.7) #, label='FCC lattice'

    # Set labels and title for the plot
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('FCC Lattice')
    
    '''
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_zlim(0, L)
    '''
    
    # Add legend and grid to the plot
    ax.legend()
    ax.grid(True, linestyle='dotted', linewidth=0.5, alpha=0.5)

    # Display the plot
    plt.show()   


def main():
    #pair-correlation
    plot_fcc_lattice(positions[1])
    
    hist_plots = np.arange(len(positions))
    pressures = np.array([])
    
    times = np.arange(len(positions))
    times = np.arange(20,40)
    
    average_hist(positions,times,bin_size)
    
    '''
    for i in hist_plots:
        pair_correlation_plot(positions[i])
        
        print('pressure is:',i,pressure(positions[i], 0.8))
        pressures = np.append(pressures,pressure(positions[i],0.8))
        
    plt.plot(pressures)
    plt.show()
    '''
    

                


    # pressure
    #  for i in hist_plots():
    #    print(pressure(positions[i], 0.8))

    

    
    

if __name__ == "__main__":
    main()
    
#%% test

n = np.array([1,5,-7,2])
r = np.array([1,2,3,4])
g = n/(r**2*0.1)
print(n % 3)


