# Simulating Argon by Natan van Steenis and Vincent van Rie
This Python script simulates the molecular dynamics of Argon atoms in a 3D box using the Verlet algorithm. The atoms are placed in a face-centered cubic (FCC) lattice and are given random velocities. The system is periodic in all directions, and the Lennard-Jones potential is used to calculate the forces between the particles.

## Start-up

 - Clone this repo in your desired directory
 - Go to "Simulating Argon.py"
 - Set the desired parameters at the top of the file
   - Set the density and temperature
   - Set the number of particles (only works for N = 4 * n^3; n an integer)
   - Set the perimeter_parameter to set how many neighbouring cells the force is calculated
   - Set the number of iterations
   - Set the bool to display the system at each timestep
   - Set the bool to save the data
 - Optionally one can use %matplotlib inline/qt when running as a module


## Authors
- Natan van Steenis
- Vincent van Rie
## Date
22/03/2024

## Dependencies
- Python 3.x
- Matplotlib
- NumPy
