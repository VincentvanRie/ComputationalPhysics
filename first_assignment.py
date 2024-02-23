import numpy as np
import matplotlib.pyplot as plt

#set the parameters
L = 100        #box size
h = 1           #time stamp
N = 4   #number of particles
m = 1           #particle mass 
epsilon = 1         # value for epsilon
sigma = 1       # value for sigma


def random_pos():
    return np.random.uniform(0,L,1)
def random_v():
    return np.random.uniform(-10,10,1)

#create the argon class 
class argon:
    def __init__(self, x,y,vx,vy):
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy

#Make the function the create N particles
def create_particles():
    particles = []
    for _ in range(N):
        x = random_pos()
        y = random_pos()
        vx = random_v()
        vy = random_v()

        particle = argon(x, y, vx, vy)
        particles.append(particle)

    return np.array(particles)

particles = create_particles()

#calculate the force between particles
def Force(r):
    F = -48*epsilon*sigma**12*r**(-13)+24*epsilon*sigma**6*r**(-7)
    return F

# Calculate the distance between two particles with particles as input 
def distance(p1,p2):
    r = np.sqrt((p1.x-p2.x)**2+(p1.y-p2.y)**2)
    #print('distance =',r)
    return r

# find the nearest particle 
def find_nearest_particle(particle, particles):
    positions = np.array([[p.x, p.y] for p in particles])
    #print(positions)
    particle_position = np.array([particle.x, particle.y])

    distances = np.linalg.norm(positions - particle_position, axis=1)
    min_distance_index = np.argmin(distances[distances > 0])
    
    
    return particles[min_distance_index]



# This function calculates the position and velocity for the next time stamp
def evolve(particle_array,h):
    particles = particle_array
    for particle in particles:
        # new positions
        xt = particle.x + particle.vx*h
        yt = particle.y + particle.vy*h
        
        #set boundary conditions:
        xt = particle.x % L
        yt = particle.y % L
        
        # new velocities
        nearest = find_nearest_particle(particle, particles)
        vx2 = particle.vx - (1/m)*Force(distance(particle,nearest))
        vy2 = particle.vy - (1/m)*Force(distance(particle,nearest))

        particle.x = xt
        particle.y = yt
        particle.vx = vx2
        particle.vy = vy2
        
    return particles  
        

        
        
# Function to plot the particles
def plot_particles(particles):
    plt.figure(figsize=(8, 8))
    plt.scatter([particle.x for particle in particles], [particle.y for particle in particles], marker='o')
    plt.title('Particle Positions')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.xlim(0, L)
    plt.ylim(0, L)
    plt.grid(True)
    plt.show()

def plotter(particle_array,iterations):
    for i in np.arange(iterations):
        updated_particles = evolve(particle_array,h) # update the particles
        for particle in updated_particles: #plot the particles
                plt.scatter(particle.x,particle.y)
        plt.title(i)
        plt.show()
        


def main(): 
    
    particles = create_particles()
    
    for i in np.arange(10):
        for particle in particles: #plot the particles
            plt.scatter(particle.x,particle.y)
        plt.title(i)
        plt.show()
        
        #update the particles list
        particles = evolve(particles, h)
    
    #for particle in particles:
        #print(f"Particle at ({particle.x}, {particle.y}) with velocity ({particle.vx}, {particle.vy})")
    
    # Plot particle positions
    #plot_particles(particles)
    
    #plotter(particles,10)
    
    return 0
    
main()


'''    
    
#%%
def random_pos():
    return np.random.uniform(0,L,1)
def random_v():
    return np.random.uniform(-np.pi,np.pi,1)


class argon:
    def __init__(self, x,y,vx,xy):
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
        
for i in np.arange(3):
    i = argon(random_pos(),random_pos(),)

    
def main():
    
'''

