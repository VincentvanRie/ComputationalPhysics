import math
import random
import matplotlib.pyplot as plt
import numpy as np


# class Box:
#     def __init__(self, L, N):
#         self.L = L
#         self.N = N

#define the consants
''' 
boltzmann = 1.38064852 * 10**-23
T = 20
m = 39.948 * 1.66053906660 * 10**-27
epsilon = 119.8 * boltzmann
sigma = 3.405 * 10**-10
'''
  
N = 18
L = 2 # in units of sigma

h = 0.0001
reps = 500
m=1


class System_of_particles:
    def __init__(self, positions, velocities):
        self.positions = positions  # positions = {x: [,,,], y:[,,,]}
        self.velocities = velocities # velocities = {x: [,,,], y:[,,,]}
        #self.previous_positions = positions
        
        
    def Ek(self): #take the system class as input and calculate the total kinetic energy
        E_kinetic = np.array([])
        #print(self.velocities['x'][0])
        for i in np.arange(N): #calculate the kinetic energy of each particle
            
            vel = np.sqrt(self.velocities['x'][i]**2 + self.velocities['y'][i]**2)
            E_kin = 1/2 * m * vel**2
            #print(E_kin)
            
            E_kinetic = np.append(E_kinetic,E_kin)
            
        total_Ek = np.sum(E_kinetic) 
        total_Ek
        return total_Ek
    
    def Ep(self): #calculate the potential energy
        E_potential = np.array([])
        
        for i in np.arange(N): #calculate potential energy of each particle     #TODO potential energy j = i+1
            potential_i = np.array([]) # calculate the potential between particle i and all other particles
            for j in range(N):
                if i != j:
                    delta_x = (self.positions["x"][i] - self.positions["x"][j] + L/2) % L - L/2
                    delta_y = (self.positions["y"][i] - self.positions["y"][j] + L/2) % L - L/2
                    
                    r = ((delta_x) ** 2 + (delta_y) ** 2) ** 0.5
                
                    E_pot = ljpot_til(r)
                
                
                    E_potential = np.append(E_potential,E_pot)
                    #print(E_pot)
            total_i= np.sum(potential_i)
            E_potential = np.append(E_potential,total_i)
            
        total_potential = np.sum(E_potential)
        return total_potential/2
    

def initial_setup(): #function to define the initial setup returns the position and velocity dictionary
    '''
    starting_position = {
        "x": np.linspace(0, L, N, endpoint=False),
        "y": np.linspace(0, L, N, endpoint=False)}
    '''

    x1 = np.array([0, 1 / 3 * L, 2 / 3 * L] * 3)
    print(x1)
    y1 = np.array([0] * 3 + [1 / 3 * L] * 3 + [2 / 3 * L] * 3)

    x2 = x1 + np.repeat(1 / 6 * L, 9)
    y2 = y1 + np.repeat(1 / 6 * L, 9)

    x = np.concatenate((x1, x2))
    y = np.concatenate((y1, y2))
    
    starting_position = {"x": x, "y":y}
    
    angles_in_vel = np.array([random.uniform(-np.pi, np.pi) for _ in range(N)])
    
    
    starting_velocity = {"x": 10*np.array([np.cos(i) for i in angles_in_vel]),
                         "y": 10*np.array([np.sin(i) for i in angles_in_vel])} # nog aanpassen vincent code #TODO
    return starting_position, starting_velocity

def ljpot_til(r): #Lennard jones potential for the dimensionless units
    potential = 4*(r**(-12)-r**(-6))    
    return potential

def force_til(r): #calculate the force for the new variables, input is difference between two particles
    force = m * 24*(2*r**(-13)-r**(-7)) 
    return force


def Calculate_force(positions):
    z_dimension = False
    if "z" in positions:
        z_dimension = True

    force_dict = {
        "x": np.zeros(N),
        "y": np.zeros(N),
        "z": np.zeros(N),
    }

    for i in range(N):
        for j in range(N):
            if i != j:
                delta_x = (positions["x"][i] - positions["x"][j] + L/2) % L - L/2
                delta_y = (positions["y"][i] - positions["y"][j] + L/2) % L - L/2
                
                
                r = np.sqrt((delta_x) ** 2 + (delta_y) ** 2)
                if z_dimension:
                    delta_z = positions["z"][i] - positions["z"][j]
                    r = ((r) ** 2 + (delta_z) ** 2) ** 0.5
                angle = math.atan2(delta_y, delta_x)

                force = force_til(r)

                force_dict["x"][i] += force * math.cos(angle)
                force_dict["y"][i] += force * math.sin(angle)
                if z_dimension:
                    force_dict["z"][i] += force * math.sin(angle)  # TODO

    return force_dict
'''
Natans work using numpy vectors
def Calculate_force(positions): #input should be system.positions which is a dict 
    force_dict = {'x': np.zeros(N), 'y': np.zeros(N)}
    
    for i in range(N):
        x_pos = positions['x'][i]
        y_pos = positions['y'][i]
        #calculate the distance to the each closest particle in the min. image conv.
        
        delta_x = (x_pos - positions['x'] + L/2) % L - L/2
        delta_y = (y_pos - positions['y'] + L/2) % L - L/2
        #print(closest_x)
        #print(np.shape(closest_y))
        
        r = ((delta_x) ** 2 + (delta_y) ** 2) ** 0.5
        
        
        #print(f'distance between particle {i}: {r}')
        #r is now an array with the distance from particle i to each particle 
        

        angle = math.atan2(delta_y, delta_x)

        potential = 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)
        force = (100 * 24 * epsilon * (2 * (sigma / r) ** 13 - (sigma / r) ** 7) / sigma )

        force_dict["x"][i] += force * math.cos(angle)
        force_dict["y"][i] += force * math.sin(angle)


    return force_dict
 '''   
#This function takes the system class, updates it and then returns the updated system
def updated_system(system):
    x_pos = system.positions['x']
    y_pos = system.positions['y']
    x_vel = system.velocities['x']
    y_vel = system.velocities['y']
    
    #calculate the force on all the particles:
    current_force = Calculate_force(system.positions)
    #print(current_force['x'])
    
    #calculate the new positions:
    new_x_pos = x_pos + h*x_vel + h**2/(2*m)*current_force['x']
    new_y_pos = y_pos + h*y_vel + h**2/(2*m)*current_force['y']
       
    # enfore boundary conditions: 
    new_x_pos = new_x_pos % L
    new_y_pos = new_y_pos % L
    
    new_positions = {'x': new_x_pos, 'y': new_y_pos}
    
    #calculate the updated velocity(v(t+h)): 
    #first we need the force of the updated particles (t+h)
    updated_force = Calculate_force(new_positions)
    new_x_vel = x_vel + h/(2*m)*(updated_force['x'] + current_force['x'])
    new_y_vel = y_vel + h/(2*m)*(updated_force['y'] + current_force['y'])
    
    # make a new dict for the updated velocities
    new_velocities = {'x': new_x_vel, 'y': new_y_vel}
    
    updated_system = System_of_particles(new_positions, new_velocities)
    
    return updated_system 

def system_plotter(system): #takes the system as input and plots it
        plt.scatter(system.positions['x'],system.positions['y'])
        
        plt.xlim(0, L)
        plt.ylim(0, L)
        
        plt.quiver(system.positions["x"], system.positions["y"],
                system.velocities["x"], system.velocities["y"],)
        plt.xlabel('x position in (\u03C3)')
        plt.ylabel('y position (\u03C3)')
        plt.show()
        
    


    
def main():
    system = System_of_particles(initial_setup()[0], initial_setup()[1])      #create the intitial system
    #print(system.velocities['x'])
    x_positions = np.array([system.positions['x']])
    y_positions = np.array([system.positions['y']])
    x_vel = np.array([system.velocities['x']])
    y_vel = np.array([system.velocities['y']])
    
    E_kin = np.array([]) #create array of kinetic energies
    E_pot = np.array([]) #create array of potential energies

    
    for i in np.arange(reps):
        if i == 0:
            system_plotter(system)

        else:
            system = updated_system(system)
            system_plotter(system) 
            #print(system.Ek())
            E_kin = np.append(E_kin, system.Ek())
            E_pot = np.append(E_pot, system.Ep())
            #print(system.Ek())
    #print(E_pot)
    
    E_tot = E_kin + E_pot
    
    plt.scatter(np.arange(reps-1),E_kin,label = 'kinetic',s=1)
    plt.scatter(np.arange(reps-1),E_pot, label = 'potential',s=1)
    plt.scatter(np.arange(reps-1),E_tot, label = 'total',s=1)
    plt.title('Energies of the system')
    plt.legend()
    plt.xlabel('x position in (\u03C3)')
    plt.ylabel('y position in (\u03C3)')
    plt.show()
    
if __name__ == "__main__":
    main()

#%% make fcc lattice 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def fcc_lattice(a):
    # Define the positions of particles in a single FCC unit cell
    ground_layer = np.array([[0, 0, 0], [a, 0, 0], [0, a, 0], [a, a, 0], [a/2, a/2, 0]])
    middle_layer = np.array([[a/2, 0, a/2], [0, a/2, a/2], [a, a/2, a/2], [a/2, a, a/2]])
    upper_layer = np.array([[0, 0, a], [a, 0, a], [0, a, a], [a, a, a], [a/2, a/2, a]])

    # Combine ground and upper layers to get the complete FCC lattice
    positions = {
        'x': np.concatenate((ground_layer[:, 0], middle_layer[:,0], upper_layer[:, 0])),
        'y': np.concatenate((ground_layer[:, 1],middle_layer[:,1], upper_layer[:, 1])),
        'z': np.concatenate((ground_layer[:, 2],middle_layer[:,2], upper_layer[:, 2]))
    }

    return positions

def create_multi_box_lattice(a, repetitions):
    # Generate positions for a multi-box FCC lattice by repeating the single FCC unit cell
    lattice = fcc_lattice(a)
    multi_box_positions = {'x': [],'y': [],'z': []
    }

    for i in range(repetitions[0]):
        for j in range(repetitions[1]):
            for k in range(repetitions[2]):
                # Calculate offsets for each repetition
                offset_x = i * a
                offset_y = j * a
                offset_z = k * a

                # Add the positions of the current FCC unit cell with offset to the multi-box lattice
                multi_box_positions['x'] = np.concatenate((multi_box_positions['x'], lattice['x'] + offset_x))
                multi_box_positions['y'] = np.concatenate((multi_box_positions['y'], lattice['y'] + offset_y))
                multi_box_positions['z'] = np.concatenate((multi_box_positions['z'], lattice['z'] + offset_z))

    return multi_box_positions

def plot_fcc_lattice(positions):
    # Plot the FCC lattice
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Scatter plot for FCC lattice
    ax.scatter(positions['x'], positions['y'], positions['z'], c='b', marker='o', s=2, alpha=0.7, label='FCC lattice')

    # Set labels and title for the plot
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('FCC Lattice')

    # Add legend and grid to the plot
    ax.legend()
    ax.grid(True, linestyle='dotted', linewidth=0.5, alpha=0.5)

    # Display the plot
    plt.show()

# Example usage:
reps = 3
repetitions = (reps,reps ,reps )
multi_box_positions = create_multi_box_lattice(L, repetitions)

# Plot the multi-box FCC lattice
plot_fcc_lattice(multi_box_positions)


