import numpy as np
import matplotlib.pyplot as plt
print("USE SI UNITS ONLY FOR ALL INPUTS")

#initial conditions and parameters
v0 = float(input("Enter initial velocity of the body: "))
if v0 >= 299792458 or v0 <= -299792458:
    print("Velocity cannot exceed the speed of light. Please enter a valid velocity.")
    exit()
    
x0 = float(input("Enter initial position of the body: "))
m = float(input("Enter mass of the body: "))
q = float(input("Enter charge of the body: "))
E = float(input("Enter magnitude of the electric field: "))
steps = int(input("Enter the number of time steps for the simulation: "))*100
dt = 0.01

#Arrays
t = np.linspace(0, (steps-1) * dt, steps)
x = np.zeros(len(t))
v = np.zeros(len(t))
gamma = np.zeros(len(t))
a = np.zeros(len(t))

#Initial Values
x[0] = x0
v[0] = v0

#RK4 Integration Method
for i in range(1, len(t)):
    gamma[i-1] = 1 / np.sqrt(1 - (v[i-1] / 299792458) ** 2)
    a[i-1] = (q * E) / (m * gamma[i-1] ** 3)
    
    k1_v = a[i-1] * dt
    k1_x = v[i-1] * dt
    
    k2_v = ((q * E) / (m * (gamma[i-1] + 0.5 * k1_v) ** 3)) * dt
    k2_x = (v[i-1] + 0.5 * k1_v) * dt
    
    k3_v = ((q * E) / (m * (gamma[i-1] + 0.5 * k2_v) ** 3)) * dt
    k3_x = (v[i-1] + 0.5 * k2_v) * dt
    
    k4_v = ((q * E) / (m * (gamma[i-1] + k3_v) ** 3)) * dt
    k4_x = (v[i-1] + k3_v) * dt
    
    v[i] = v[i-1] + (k1_v + 2*k2_v + 2*k3_v + k4_v) / 6
    x[i] = x[i-1] + (k1_x + 2*k2_x + 2*k3_x + k4_x) / 6
    
#Plotting of x-t graph of particle motion
plt.plot(t, x)
plt.title("Relativistic Charged Particle in Electric Field Constrained to Move in One Dimension")
plt.xlabel("Time (s)")
plt.ylabel("Position (m)")
plt.grid()
plt.show()


#Plotting of v-t graph of particle motion
plt.plot(t, v)
plt.title("Relativistic Charged Particle in Electric Field Constrained to Move in One Dimension")
plt.xlabel("Time (s)")
plt.ylabel("Velocity (m/s)")
plt.grid()
plt.show()