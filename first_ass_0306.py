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

#define the consants
boltzmann = 1.38064852 * 10**-23
T = 20
m = 39.948 * 1.66053906660 * 10**-27
epsilon = 119.8 * boltzmann
sigma = 3.405 * 10**-10

class System_of_particles:
    def __init__(self, positions, velocities):
        self.positions = positions  # positions = {x: [,,,], y:[,,,]}
        self.velocities = velocities # velocities = {x: [,,,], y:[,,,]}
        #self.previous_positions = positions
        

def initial_setup(): #function to define the initial setup returns the position and velocity dictionary
    starting_position = {
        "x": np.array([random.uniform(0, L) for _ in range(N)]),
        "y": np.array([random.uniform(0, L) for _ in range(N)])}
    starting_velocity = {"x": np.array([np.cos(random.uniform(-np.pi, np.pi)) for i in range(N)]),
                         "y": np.array([np.sin(random.uniform(-np.pi, np.pi)) for i in range(N)])}
    return starting_position, starting_velocity


        

def Calculate_force(positions): #input should be system.positions which is a dict 
    force_dict = {'x': np.zeros(N), 'y': np.zeros(N)}
    
    for i in range(N):
        for j in range(N):
            if i != j:
                delta_x = positions["x"][i] - positions["x"][j]
                delta_y = positions["y"][i] - positions["y"][j]
                r = ((delta_x) ** 2 + (delta_y) ** 2) ** 0.5

                if r == 0:
                    print("Zero")

                angle = math.atan2(delta_y, delta_x)

                potential = 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)
                force = (100 * 24 * epsilon * (2 * (sigma / r) ** 13 - (sigma / r) ** 7) / sigma )

                force_dict["x"][i] += force * math.cos(angle)
                force_dict["y"][i] += force * math.sin(angle)


    return force_dict
    
#This function takes the system class, updates it and then returns the updated system
def updated_system(system):
    x_pos = system.positions['x']
    y_pos = system.positions['y']
    x_vel = system.velocities['x']
    y_vel = system.velocities['y']
    
    #calculate the force on all the particles:
    current_force = Calculate_force(system.positions)
    print(current_force['x'])
    
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
        plt.show()

    

    


    
def main():
    N = 5 #number of particles
    reps = 100 #number of iterations of the simulation
    h = 1
    system = System_of_particles(initial_setup()[0], initial_setup()[1])            #create the intitial system
    

    for i in np.arange(reps):
        if i == 0:
            system_plotter(system)
        else:
            system = updated_system(system)
            system_plotter(system)
    
    
if __name__ == "__main__":
    main()

#%% 
    def Change_positions(self):
        z_dimension = False

        self.previous_positions["x"] = self.positions["x"]
        self.previous_positions["y"] = self.positions["y"]
        if "z" in self.positions:
            z_dimension = True
            self.previous_positions["z"] = self.positions["z"]
        # self.previous_position = [self.particles[i]["position"] for i in range(len(self.particles))]

        current_force = Calculate_force(self.previous_positions)

        # current_force = [self.Calculate_force(self.particles[i]["position"]) for i in range(len(self.particles))]

        self.positions["x"] = (self.positions["x"] + self.velocities["x"] * h
            + current_force["x"] * h**2 / (2 * m))
        self.positions["y"] = (self.positions["y"] + self.velocities["y"] * h
            + current_force["y"] * h**2 / (2 * m))
        if z_dimension:
            self.positions["z"] = (
                self.positions["z"]
                + self.velocities["z"] * h
                + current_force["z"] * h**2 / (2 * m)
            )
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
        ) * h / (2 * m)
        self.velocities["y"] = self.velocities["y"] + (
            force["y"] + current_force["y"]
        ) * h / (2 * m)
        if z_dimension:
            self.velocities["z"] = self.velocities["z"] + (
                force["z"] + current_force["z"]
            ) * h / (2 * m)

        if time < 500 * h:
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

        return self.positions, self.velocities


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
                delta_x = positions["x"][i] - positions["x"][j]
                delta_y = positions["y"][i] - positions["y"][j]
                r = ((delta_x) ** 2 + (delta_y) ** 2) ** 0.5
                if z_dimension:
                    delta_z = positions["z"][i] - positions["z"][j]
                    r = ((r) ** 2 + (delta_z) ** 2) ** 0.5
                if r == 0:
                    print("Zero")
                angle = math.atan2(delta_y, delta_x)

                potential = 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)
                force = (100 * 24 * epsilon * (2 * (sigma / r) ** 13 - (sigma / r) ** 7) / sigma )

                force_dict["x"][i] += force * math.cos(angle)
                force_dict["y"][i] += force * math.sin(angle)
                if z_dimension:
                    force_dict["z"][i] += force * math.sin(angle)  # TODO

    return force_dict


def Lambda(velocities):
    return np.sqrt((3 * (N - 1) * boltzmann * T) / (np.sum(velocities**2) * m))


def main():
    global m, epsilon, sigma, L, h, boltzmann, T, N
    boltzmann = 1.38064852 * 10**-23
    T = 20
    m = 39.948 * 1.66053906660 * 10**-27
    epsilon = 119.8 * boltzmann
    sigma = 3.405 * 10**-10
    L = 100 * sigma
    h = 10**-14

    N = 5
    v = 401.0484
    timesteps = 5000
    t_end = timesteps * h

    velocities = np.array([random.uniform(0, v) for _ in range(N)])

    velocities = Lambda(velocities) * velocities

    Argon_system = System_of_particles(
        positions={
            "x": np.array([random.uniform(0, L) for _ in range(N)]),
            "y": np.array([random.uniform(0, L) for _ in range(N)]),
        },
        velocities={
            "x": np.array(
                [
                    np.cos(random.uniform(-np.pi, np.pi)) * velocities[i]
                    for i in range(N)
                ]
            ),
            "y": np.array(
                [
                    np.sin(random.uniform(-np.pi, np.pi)) * velocities[i]
                    for i in range(N)
                ]
            ),
        },
    )

    positions_over_time = [Argon_system.positions]
    velocities_over_time = [Argon_system.velocities]

    global time
    for time in np.arange(0, t_end, h):
        Argon_system.Change_positions()
        positions = Argon_system.positions
        velocities = Argon_system.velocities
        positions_over_time.append(positions)
        velocities_over_time.append(velocities)
        plt.clf()
        plt.xlim(0, L)
        plt.ylim(0, L)
        plt.quiver(
            Argon_system.positions["x"],
            Argon_system.positions["y"],
            Argon_system.velocities["x"],
            Argon_system.velocities["y"],
        )
        plt.scatter(
            Argon_system.previous_positions["x"],
            Argon_system.previous_positions["y"],
            alpha=0.5,
        )
        plt.pause(0.05)

    plt.show()
    plt.figure()
    plt.xlim(0, L)
    plt.ylim(0, L)
    for i in range(N):
        plt.scatter(
            [positions_over_time[j]["x"][i] for j in range(len(positions_over_time))],
            [positions_over_time[j]["y"][i] for j in range(len(positions_over_time))],
            color="blue",
            alpha=0.5,
        )
    plt.show()

    plt.figure()
    plt.xlim(0, L)
    plt.ylim(0, L)
    plt.quiver(
        Argon_system.positions["x"],
        Argon_system.positions["y"],
        Argon_system.velocities["x"],
        Argon_system.velocities["y"],
    )
    plt.scatter(
        Argon_system.previous_positions["x"],
        Argon_system.previous_positions["y"],
        alpha=0.5,
    )
    plt.show()

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