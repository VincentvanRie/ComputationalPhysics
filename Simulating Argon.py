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


class Argon:
    def __init__(self, position, velocity, F, time=[0]):
        self.position = position
        self.velocity = velocity
        self.F = F
        self.t = time

    def Evolve(self):
        # self.position["x"] += self.velocity["x"] * dt
        # self.position["y"] += self.velocity["y"] * dt

        # self.velocity["x"] += self.F["x"] * dt / m
        # self.velocity["y"] += self.F["y"] * dt / m

        self.position["x"] += (
            2 * self.velocity["x"] * dt + 0.5 * self.F["x"] * dt**2 / m
        )

        if self.position["x"] > L or self.position["x"] < 0:
            # if self.position["x"] > 2 * L or self.position["x"] < -L:
            #     print(
            #         f"Particle has left the adjacent boxes {round(self.position["x"]/L,2)}. Please reduce the time interval or increase the box size."
            #     )
            self.position["x"] = (self.position["x"] % L + L) % L

        if self.position["y"] > L or self.position["y"] < 0:
            # if self.position["y"] > 2 * L or self.position["y"] < -L:
            #     print(
            #         f"Particle has left the adjacent boxes {round(self.position["y"]/L)}. Please reduce the time interval or increase the box size."
            #     )
            self.position["y"] = (self.position["y"] % L + L) % L

    def Potential(self, argon_companion):
        delta_x = argon_companion.position["x"] - self.position["x"]
        delta_y = argon_companion.position["y"] - self.position["y"]
        r = ((delta_x) ** 2 + (delta_y) ** 2) ** 0.5
        if r == 0:
            return {"x": 0, "y": 0}
        angle = math.atan2(delta_y, delta_x)

        potential = 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)
        force = 100 * 24 * epsilon * (2 * (sigma / r) ** 13 - (sigma / r) ** 7) / sigma

        return {"x": force * math.cos(angle), "y": force * math.sin(angle)}


def main():
    global m, epsilon, sigma, L, dt
    boltzmann = 1.38064852 * 10**-23
    m = 1  # 39.948 * 1.66053906660 * 10**-27
    epsilon = 1  # 119.8 * boltzmann
    sigma = 1  # 3.405 * 10**-10
    L = 100  # 100 * sigma
    dt = 0.01

    N = 100
    v = 1  # 10 * 10**-50
    t_end = 100 * dt

    argons = [
        Argon(
            position={"x": random.random() * L, "y": random.random() * L},
            velocity={
                "x": v * (random.random() * 2 - 1),
                "y": v * (random.random() * 2 - 1),
            },
            F={"x": 0, "y": 0},
        )
        for _ in range(N)
    ]

    positions_over_time = {
        0: argons,
    }

    for time in np.arange(dt, t_end, dt):
        print(time)
        particles = copy.deepcopy(
            positions_over_time[
                round(time - dt, 10)  # round to avoid floating point errors
            ]
        )  # Copy the argons from the previous time step.

        for argon in particles:
            argon.F = {"x": 0, "y": 0}

            for argon2 in particles:
                if argon != argon2:
                    argon.F["x"] += argon.Potential(argon2)["x"]
                    argon.F["y"] += argon.Potential(argon2)["y"]

        for argon in particles:
            argon.Evolve()

        positions_over_time[
            round(time, 10)
        ] = particles  # round to avoid floating point errors

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

    # Make an animation of the positions of the particles over time
    fig, ax = plt.subplots()
    line1 = ax.quiver([], [], [], [], [], clim=[-np.pi, np.pi])

    def init():
        ax.set_xlim(0, L)
        ax.set_ylim(0, L)
        return (line1,)

    def update(frame, ln, x, y):
        x.append(frame)
        y.append(frame)
        ln.set_data(x, y)
        return (ln,)

    def animate(i):
        time = round(i * dt, 10)
        x = [positions_over_time[time][i_argon].position["x"] for i_argon in range(N)]
        y = [positions_over_time[time][i_argon].position["y"] for i_argon in range(N)]

        cos = [positions_over_time[time][i_argon].velocity["x"] for i_argon in range(N)]
        sin = [positions_over_time[time][i_argon].velocity["y"] for i_argon in range(N)]
        orient = [0 for _ in range(N)]

        line1.set_offsets(np.c_[x, y])
        line1.set_UVC(np.c_[cos], np.c_[sin], np.c_[orient])

        return (line1,)

    ani = FuncAnimation(
        fig,
        animate,
        frames=np.arange(0, int(t_end / dt)),
        init_func=init,
        blit=True,
    )
    plt.show()

    exit()
    # Make an animation of the positions of the particles over time
    fig, ax = plt.subplots()
    line1 = ax.quiver([], [], [], [], [], clim=[-np.pi, np.pi])

    def init():
        ax.set_xlim(0, L)
        ax.set_ylim(0, L)
        return (line1,)

    def update(frame, ln, x, y):
        x.append(frame)
        y.append(frame)
        ln.set_data(x, y)
        return (ln,)

    def animate(i):
        time = round(i * dt, 10)
        x = [positions_over_time[time][i_argon].position["x"] for i_argon in range(N)]
        y = [positions_over_time[time][i_argon].position["y"] for i_argon in range(N)]

        cos = [positions_over_time[time][i_argon].velocity["x"] for i_argon in range(N)]
        sin = [positions_over_time[time][i_argon].velocity["y"] for i_argon in range(N)]
        orient = [0 for _ in range(N)]

        line1.set_offsets(np.c_[x, y])
        line1.set_UVC(np.c_[cos], np.c_[sin], np.c_[orient])

        return (line1,)

    ani = FuncAnimation(
        fig,
        animate,
        frames=np.arange(0, int(t_end / dt)),
        init_func=init,
        blit=True,
    )

    plt.show()


if __name__ == "__main__":
    main()
