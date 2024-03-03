import numpy as np

init_dict = np.array([0, 1, 2, 3])

squared = init_dict**2

if any(squared == 0):
    print("Zero")
print(squared)
print(squared == 0)

squared[squared == 0] = 1

print(squared)


positions = {"x": np.ones(3), "y": [], "z": []}
force = {"x": [], "y": [], "z": []}

for key in positions:
    print(positions[key])


# If z exists in the dictionary, print "Z exists"
if "z" in positions:
    print("Z exists")

positions["x"][1] = 0
print(positions["x"])

if any(positions["x"] == 0):
    print("Zero")


boltzmann = 1.38064852 * 10**-23
T = 260
m = 39.948 * 1.66053906660 * 10**-27
N = 108


v = np.sqrt(3 * boltzmann * T * (N - 1) / (N * m))

print(v)
