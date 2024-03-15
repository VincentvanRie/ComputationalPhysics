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
  
N = 108
L = 2.2 # in units of sigma the size of the system

v_start = 10000
h = 0.000001
reps = 50

#time in units of \sqrt{\frac{m\sigma^2}{\epsilon}}

'''
units: 
    x in [\sigma]
    t in [(m\sigma^2)^1/2]
    

'''



class System_of_particles:
    def __init__(self, positions, velocities):
        self.positions = positions  # positions = {x: [,,,], y:[,,,], z:[,,,]}
        self.velocities = velocities # velocities = {x: [,,,], y:[,,,], z:[,,,]}
        #self.previous_positions = positions
        
        
    def Ek(self): #take the system class as input and calculate the total kinetic energy
        E_kinetic = np.array([])
        #print(self.velocities['x'][0])
        for i in np.arange(N): #calculate the kinetic energy of each particle
            
            vel = np.sqrt(self.velocities['x'][i]**2 + self.velocities['y'][i]**2 + self.velocities['z'][i]**2)
            E_kin = 1/2 * vel**2
            #print(E_kin)
            
            E_kinetic = np.append(E_kinetic,E_kin)
            
        total_Ek = np.sum(E_kinetic) 
        total_Ek
        return total_Ek
    
    def Ep(self): #calculate the potential energy
        E_potential = np.array([])
        
        
        for i in np.arange(N): #calculate potential energy of each particle     #TODO potential energy j = i+1
            potential_i = np.array([]) # calculate the potential between particle i and all other particles
            for j in range(i,N):
                if i != j:
                    delta_x = (self.positions["x"][i] - self.positions["x"][j] + L/2) % L - L/2
                    delta_y = (self.positions["y"][i] - self.positions["y"][j] + L/2) % L - L/2
                    delta_z = (self.positions["z"][i] - self.positions["z"][j] + L/2) % L - L/2
                    
                    r = ((delta_x) ** 2 + (delta_y) ** 2 + (delta_z) ** 2) ** 0.5
                    #print('delta x =', delta_x, 'delta y =', delta_y, 'delta z =', delta_z, )
                    #print('distance =', r)
                
                    E_pot = ljpot_til(r)
                
                
                    E_potential = np.append(E_potential,E_pot)
                    #print(E_pot)
            total_i= np.sum(potential_i)
            E_potential = np.append(E_potential,total_i)
        
        print(E_potential)
            
        total_potential = np.sum(E_potential)
        return total_potential/2


def PrepateLattice():
    x1 = np.array([0, 1 / 3 * L, 2 / 3 * L] * 9)
    y1 = np.array(([0] * 3 + [1 / 3 * L] * 3 + [2 / 3 * L] * 3) * 3)
    z1 = np.repeat([0, 1 / 3 * L, 2 / 3 * L], 9)

    x2 = x1 + np.repeat(1 / 6 * L, 27)
    y2 = y1 + np.repeat(1 / 6 * L, 27)
    z2 = z1

    x3 = x1 + np.repeat(1 / 6 * L, 27)
    y3 = y1
    z3 = z1 + np.repeat(1 / 6 * L, 27)

    x4 = x1
    y4 = y1 + np.repeat(1 / 6 * L, 27)
    z4 = z1 + np.repeat(1 / 6 * L, 27)

    x = np.concatenate((x1, x2, x3, x4))
    y = np.concatenate((y1, y2, y3, y4))
    z = np.concatenate((z1, z2, z3, z4))
    
    positions = {"x": x, "y": y, "z": z}
    
    
    
    vel_XY = np.array([random.uniform(-np.pi, np.pi) for _ in range(N)])
    vel_Z = np.array([random.uniform(-np.pi, np.pi) for _ in range(N)])

    starting_velocity = {"x": v_start * np.array([np.cos(i) for i in vel_XY]),
                         "y": v_start * np.array([np.sin(i) for i in vel_XY]),
                         "z": v_start * np.array([np.sin(i) for i in vel_Z])
                         }
    
    #initialize some test positions
    #x_test = np.array[0,0.3*L,0.6*L,0.8*L]
    #y_test = np.array[0,0.4*L,0.7*L,0.8*L]
    #z_test = np.array[0,0.3*L,0.6*L,0.9*L]
    
    #test_positions = {"x": }
    
    


    return positions, starting_velocity


def ljpot_til(r): #Lennard jones potential for the dimensionless unit r in unit epsilon
    potential = 4*(r**(-12)-r**(-6))    
    return potential

def force_til(r): #calculate the force for the new variables, input is difference between two particles
    force = 24*(2*r**(-13)-r**(-7)) 
    return force


def Calculate_force(positions):
    
    part_pos = np.array([(x, y, z) for x, y, z in zip(positions["x"], positions["y"], positions["z"])])
    #print(part_pos)
    
    force_dict = {"x": np.zeros(N),"y": np.zeros(N),"z": np.zeros(N)}
    
    for i in range(len(positions["x"])):
        position = (positions["x"][i], positions["y"][i],positions["z"][i])
        #print('particle position = ', position )


    #print("positions = ", positions)
    
    for i in range(N):
        #print('force_dict = ',i, force_dict)
        for j in range(i,N):
            if i != j:
                delta_x = (positions["x"][i] - positions["x"][j] + L/2) % L - L/2
                #print("delta_x is ",delta_x)
                delta_y = (positions["y"][i] - positions["y"][j] + L/2) % L - L/2
                delta_z = (positions["z"][i] - positions["z"][j] + L/2) % L - L/2
                
                #print(delta_x,delta_y,delta_z)

                
                r = np.sqrt((delta_x) ** 2 + (delta_y) ** 2+ (delta_z) ** 2)
                
                #print("r = ", "i,j = ", i,j, "r=", r)
                
                XY_angle = math.atan2(delta_y, delta_x)
                Z_angle = math.atan2(delta_z, r)   # TODO shouldn't this be different?
                
                force = force_til(r)

                force_dict["x"][i] += force * math.cos(XY_angle)
                force_dict["y"][i] += force * math.sin(XY_angle)
                force_dict["z"][i] += force * math.sin(Z_angle)
              
    #print('force_dict = ', force_dict)

    return force_dict
 
#This function takes the system class, updates it and then returns the updated system
def updated_system(system):
    x_pos = system.positions['x']
    y_pos = system.positions['y']
    z_pos = system.positions['z']
    
    #print("pos = ", x_pos)
    
    
    x_vel = system.velocities['x']
    y_vel = system.velocities['y']
    z_vel = system.velocities['z']
    
    #calculate the force on all the particles:
    current_force = Calculate_force(system.positions)
    #print(current_force)
    #print(current_force)
    
    #calculate the new positions:
    new_x_pos = x_pos + h*x_vel + h**2/(2)*current_force['x']
    new_y_pos = y_pos + h*y_vel + h**2/(2)*current_force['y']
    new_z_pos = z_pos + h*z_vel + h**2/(2)*current_force['z']
       
    # enforce boundary conditions: 
    new_x_pos = new_x_pos % L
    new_y_pos = new_y_pos % L
    new_z_pos = new_z_pos % L
    
    new_positions = {'x': new_x_pos, 'y': new_y_pos, 'z': new_z_pos}
    
    #calculate the updated velocity(v(t+h)): 
    #first we need the force of the updated particles (t+h)
    #print('now updated force:')
    updated_force = Calculate_force(new_positions)
    new_x_vel = x_vel + h/(2)*(updated_force['x'] + current_force['x'])
    new_y_vel = y_vel + h/(2)*(updated_force['y'] + current_force['y'])
    new_z_vel = z_vel + h/(2)*(updated_force['z'] + current_force['z'])
    
    # make a new dict for the updated velocities
    new_velocities = {'x': new_x_vel, 'y': new_y_vel, 'z': new_z_vel}
    
    updated_system = System_of_particles(new_positions, new_velocities)
    
    return updated_system 

    
random_indices = [np.random.randint(0, 101, size=10)]  
#random_indices = np.arange(10)  
    
def plot_fcc_lattice(system):
    
    positions = system.positions
    
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
    plt.quiver(system.positions["x"], system.positions["y"],system.positions["z"],
               system.velocities["x"], system.velocities["y"],system.velocities["z"], length = 0.00001)    

    # Set labels and title for the plot
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('FCC Lattice')
    
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_zlim(0, L)

    
    # Add legend and grid to the plot
    ax.legend()
    ax.grid(True, linestyle='dotted', linewidth=0.5, alpha=0.5)

    # Display the plot
    plt.show()


    
def main():
    
    system = System_of_particles(PrepateLattice()[0], PrepateLattice()[1])      #create the intitial system
    #print(system.velocities['x'])

    
    E_kin = np.array([]) #create array of kinetic energies
    E_pot = np.array([]) #create array of potential energies
    
    

    
    for i in np.arange(reps):
        if i == 0:
            #system_plotter(system)
            plot_fcc_lattice(system)
        else:
            system = updated_system(system)
            plot_fcc_lattice(system) 
            #print(system.Ek())
            E_kin = np.append(E_kin, system.Ek())
            E_pot = np.append(E_pot, system.Ep())
            #print(system.Ek())
    #print(E_pot)
    
    E_tot = E_kin + E_pot
    
    
    b = 5 # spot size
    
    plt.scatter(np.arange(reps-1),E_kin,label = 'kinetic',s=b)
    plt.scatter(np.arange(reps-1),E_pot, label = 'potential',s=b)
    plt.scatter(np.arange(reps-1),E_tot, label = 'total',s=b)
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

reps = 3
a = 10
L = a 

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
    #for i in np.arange(13):
        #print(positions["x"][i], positions["y"][i], positions["z"][i])
    #print(len(positions["x"]))

    return positions

def create_multi_box_lattice(a, repetitions):
    # Generate positions for a multi-box FCC lattice by repeating the single FCC unit cell
    lattice = fcc_lattice(a)
    multi_box_positions = {'x': [],'y': [],'z': []}

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
    
    
    part_pos = np.array([(x, y, z) for x, y, z in zip(multi_box_positions["x"], multi_box_positions["y"], multi_box_positions["z"])])
    
    # Remove duplicates from part_pos
    unique_positions = np.unique(part_pos, axis=0)

    # Reconstruct starting_positions using unique positions
    multi_box_positions = {
        'x': unique_positions[:, 0],
        'y': unique_positions[:, 1],
        'z': unique_positions[:, 2]
    }
    
    #print(unique_positions)
    
    #print(len(multi_box_positions["z"]))
    
    
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
    #ax.legend()
    ax.grid(True, linestyle='dotted', linewidth=0.5, alpha=0.5)

    # Display the plot
    plt.show()

# Example usage:

repetitions = (reps,reps ,reps )
multi_box_positions = create_multi_box_lattice(L, repetitions)

# Plot the multi-box FCC lattice
#plot_fcc_lattice(multi_box_positions)


def PrepateLattice():
    x1 = np.array([0, 1 / 3 * L, 2 / 3 * L] * 9)
    y1 = np.array(([0] * 3 + [1 / 3 * L] * 3 + [2 / 3 * L] * 3) * 3)
    z1 = np.repeat([0, 1 / 3 * L, 2 / 3 * L], 9)

    x2 = x1 + np.repeat(1 / 6 * L, 27)
    y2 = y1 + np.repeat(1 / 6 * L, 27)
    z2 = z1

    x3 = x1 + np.repeat(1 / 6 * L, 27)
    y3 = y1
    z3 = z1 + np.repeat(1 / 6 * L, 27)

    x4 = x1
    y4 = y1 + np.repeat(1 / 6 * L, 27)
    z4 = z1 + np.repeat(1 / 6 * L, 27)

    x = np.concatenate((x1, x2, x3, x4))
    y = np.concatenate((y1, y2, y3, y4))
    z = np.concatenate((z1, z2, z3, z4))
    
    part_pos = np.array([(x, y, z) for x, y, z in zip(x, y, z)])
    
    unique_positions = np.unique(part_pos, axis=0)

    # Reconstruct starting_positions using unique positions
    positions = {
        'x': unique_positions[:, 0],
        'y': unique_positions[:, 1],
        'z': unique_positions[:, 2]
    }
    
    #create initial velocities: 
    
    
    
    return positions

plot_fcc_lattice(PrepateLattice())


#%% testing
for i in np.arange(10):
    for j in np.arange(i,10):
        print(i,j)



#%% old functions 

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
 
 '''
 replaced by PrepateLattice 
 def initial_setup():
     number_of_cells = 3
     starting_positions = create_multi_box_lattice(3, (number_of_cells, number_of_cells, number_of_cells))

     vel_XY = np.array([random.uniform(-np.pi, np.pi) for _ in range(N)])
     vel_Z = np.array([random.uniform(-np.pi, np.pi) for _ in range(N)])

     starting_velocity = {"x": 10 * np.array([np.cos(i) for i in vel_XY]),
                          "y": 10 * np.array([np.sin(i) for i in vel_XY]),
                          "z": 10 * np.array([np.sin(i) for i in vel_Z])
                          }

     part_pos = np.array([(x, y, z) for x, y, z in zip(starting_positions["x"], starting_positions["y"], starting_positions["z"])])
     
     # Remove duplicates from part_pos
     unique_positions = np.unique(part_pos, axis=0)

     # Reconstruct starting_positions using unique positions
     starting_positions = PrepateLattice()


     return starting_positions, starting_velocity
 '''
 
 '''    
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
     
     #for i in np.arange(13):
       #  print(positions["x"][i], positions["y"][i], positions["z"][i])

     return positions

 def create_multi_box_lattice(a, repetitions):
     # Generate positions for a multi-box FCC lattice by repeating the single FCC unit cell
     lattice = fcc_lattice(a)
     multi_box_positions = {'x': [],'y': [],'z': []}

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

 '''
