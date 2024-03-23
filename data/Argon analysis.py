"""
Molecular Dynamics Simulation Analysis

This script analyzes the data obtained from a Molecular Dynamics (MD) 
simulation of a system of particles. The script should be run in the same 
folder as the data files.

Various analyses are performed on the data, including pair correlation function 
calculation and pressure measurement.

Authors:
- Natan van Steenis
- Vincent van Rie

Date: 22/03/2024

"""

import numpy as np
import matplotlib.pyplot as plt

# Specify which Temperature density configuration will be analysed 
T = 1
rho = 0.8


# Import the data 
positions = np.load(f'positions{T}_{rho}_long.npy', allow_pickle=True)
velocities = np.load(f'velocities{T}_{rho}_long.npy', allow_pickle=True)



# Define the system properties 
N = len(positions[0]["x"]) #number of particles
L = (N / rho) ** (1 / (3))
iterations = len(positions)  # Total number of time steps in the data



# Set the parameters for the pressure measurement
analyse_pressure = True  
num_pressures = 50 # number of times the pressure is measured
pressure_times = np.linspace(0,iterations-1,num_pressures).astype(int) #at this timesteps the pressure is measured

# Set the parameters for the hist measurement
analyse_pair_correlation = False
bin_size = 0.01
num_hist_measurements = 10
hist_times = np.linspace(500,len(positions)-1,num_hist_measurements).astype(int)






def two_part_distance(p1,p2,positions):
    
    '''
    Calculate the distance between two particles.

    Parameters:
        p1 (int): Index of first particle.
        p2 (int): Index of second particle.
        positions (dict): Dictionary containing particle positions.

    Returns:
        r (float): Distance between the particles.
    '''
    delta_x = (positions["x"][p2] - positions["x"][p1] + L / 2) % L - L / 2
    delta_y = (positions["y"][p2] - positions["y"][p1] + L / 2) % L - L / 2
    delta_z = (positions["z"][p2] - positions["z"][p1] + L / 2) % L - L / 2
    r = ((delta_x) ** 2 + (delta_y) ** 2 + (delta_z) ** 2) ** 0.5
    
    return r


def pair_correlation_plot(positions,title):
    
    '''
    Plot pair correlation histogram.

    Parameters:
        positions (dict): Dictionary containing particle positions.
        title (str): Title of the plot.
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
    plt.title(f'histogram at timestap {title}')
    plt.show()
    

    
def average_hist(all_positions,analyse_times,bin_size):
    '''
    Calculate average pair correlation histogram.

    Parameters:
        all_positions (list): List of position dictionaries at different timesteps.
        analyse_times (array-like): Timesteps to be analyzed.
        bin_size (float): Size of histogram bins.

    Returns:
        pair_correlation_g (array): Pair correlation function values.
        bins (array): Bin values.
    '''
    
    plot = True #plots the avg_hist if true
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
    
    pair_correlation_g = (2*L**3)/(N*(N-1))*avg_hist/(4*np.pi*bins[:-1]**2*bin_size)
    
    if plot:
        # Plot only the first half of the bins
        plt.bar(bins[:-1][:len(bins)//2], pair_correlation_g[:len(bins)//2],width=np.diff(bins[:len(bins)//2+1]), edgecolor='black')
        plt.xlabel('Distance $(\sigma)$')
        plt.ylabel('g(r)')
        plt.title(f'Average Pair Correlation Histogram with rho = {rho}, T = {T}')
 
        plt.savefig("pair_correlation_solid.pdf", dpi=1200)
        plt.show()
    
    return pair_correlation_g, bins
    
    
def pressure(positions,rho):
    '''
    Calculate pressure of the system.
 
    Parameters:
     positions (dict): Dictionary containing particle positions.
     rho (float): Density of the system.

     Returns:
     pressure (float): Pressure of the system.
 '''
    interaction_term = np.array([])
    
    for i in range(N):
        for j in range(i + 1, N):
            r = two_part_distance(i, j, positions)
        
            potential_der = 2 * (-24*(1 / r) ** 13 + 12* (1 / r) ** 7) 
            
            interaction_term = np.append(interaction_term,r*potential_der)
    expec_value = np.sum(interaction_term)

    pressure = T*rho - rho/(6*N) * expec_value      

    return pressure



def main():
    #pair-correlation analysis:
    if analyse_pair_correlation:
        hist_plots = np.arange(len(positions))
        average_hist(positions,hist_times,bin_size)
                
    # pressure analysis:
    if analyse_pressure:   
        pressures = np.array([])
        for i in pressure_times:
            #pair_correlation_plot(positions[i])
            pressures = np.append(pressures,pressure(positions[i],0.8))
            
        plt.plot(pressure_times,pressures)
        plt.title(f"rho is {rho} and T is {T}")
        np.save('pressures_liquid',pressures)
        plt.show()
        
        pressure_error = np.std(pressures)
        print('The error in the pressure is given by:',pressure_error)

    
    

if __name__ == "__main__":
    main()
    
#%% analyse the pressures
pres_gas = np.load('pressures_gas.npy')[500:]
pres_liquid = np.load('pressures_liquid.npy')
pres_solid = np.load('pressures_solid.npy')[500:]

averages = []



for i in range(0, len(pres_gas), 100):
    # Take the average of every 5 values
    avg = np.mean(pres_gas[i:i+5])
    # Append the average to the list
    averages.append(avg)

print('mean', np.mean(averages))
print('std', np.std(averages))


plt.plot(np.arange(500,iterations-1)*0.01,pres_gas, label = 'gas')
plt.plot(np.arange(500,iterations-1)*0.01,pres_solid, label = 'solid')
#plt.plot(np.arange(500,iterations-1)*0.01,pres_liquid, label = 'liquid')
plt.xlabel(r"Time [$\left(m \sigma^2/\epsilon\right)^{-\frac{1}{2}}]$")
plt.ylabel(r"Pressure $[\epsilon/\sigma^3]$")
plt.title('Pressures for the different argon phases')
plt.legend(loc='upper right', bbox_to_anchor=(1, 0.8))
plt.savefig("pressures.pdf", dpi=1200)
plt.show()



#%% plot the kinetic/potential/total energies
# Only works for the solid phase
import numpy as np
import matplotlib.pyplot as plt

T = 0.5
rho = 1.2


#import the data
all_positions = np.load(f'positions{T}_{rho}_long.npy', allow_pickle=True)
all_velocities = np.load(f'velocities{T}_{rho}_long.npy', allow_pickle=True)

#print(len(positions))

#print(positions[1])
'''
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(positions[1]["x"], positions[1]["y"], positions[1]["z"])
'''

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


def Calculate_potential(positions):
    """
    First the function creates a grid of cells and assigns each particle to a cell.
    and assigns each cell a list of particles that are in the cell.

    input: Dictionary with x,y,z positions {x: ... y: ... z:...}
    This function calcuates the potential between all the particles
    output: Array with the potential energy of every particle
    """
    potential_list = []
    
    N = len(positions["x"]) 

    for i in range(N):
        # cell_particle_i
        potential = 0

        for j in range(i + 1, N):
            # cell_particle_j
            r = two_part_distance(i, j, positions)
            #print(r)
            if np.isclose(r, 0):
                print("Zero")

            potential += 4  * ((1 / r) ** 12 - (1 / r) ** 6)
            #print(potential)

        potential_list.append(potential)

    return potential_list


def Calculate_kinetic(velocities):
    """
    input: Dictionary with x,y,z velocities {x: ... y: ... z:...} of every particle
    This function calcuates the kinetic energy
    outptu: Total kinetic energy of the sytem
    """
    
    
    velocities = np.sqrt(velocities["x"] ** 2 + 
                         velocities["y"] ** 2 + velocities["z"] ** 2)
   
    
    return  0.5 * np.sum(velocities**2)





#these arrays will contain the potential and kinetic energy at the different timestamps
potential_over_time = []
kinetic_over_time = []
total_over_time = []


timestamps = np.linspace(0,len(positions)-1,2000).astype(int)

for time in timestamps:
    print(time)
    #print(time)
    
    #kinetic energy
    E_kin = Calculate_kinetic(all_velocities[time])
    kinetic_over_time.append(E_kin)
    #print(E_kin)
    
    # potential energy
    E_pot = np.sum(Calculate_potential(all_positions[time]))

    #print(E_pot)
    potential_over_time.append(E_pot)
    
    total_over_time.append(E_kin+E_pot)
    
print(potential_over_time)

E_target = (N-1)*3/2*T
plt.hlines(E_target,0,len(positions)*0.01,label='equilibrium energy',color = 'purple')
plt.plot(timestamps*0.01,kinetic_over_time,label = 'E_kin')
plt.plot(timestamps*0.01,potential_over_time, label = 'E_pot')
plt.plot(timestamps*0.01,total_over_time, label = 'E_total')
#plt.ylim(0,1500)
plt.ylabel('Energy (\u03B5)')
plt.title(f'Energies during the simulation with rho = {rho} and T = {T}')
plt.xlabel(r"Time [$\left(m \sigma^2/\epsilon\right)^{-\frac{1}{2}}]$")
plt.legend()
plt.savefig("solid_energies.pdf", dpi=1200)
plt.show()