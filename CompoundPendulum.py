import numpy as np
import matplotlib.pyplot as plt

print("This program simulates the motion of a pendulum.")
print("=================================================================")
print("Assumptions: ")
print("1. The pendulum is a rigid rod constrained to raotate about a pivot with a spherical mass at the end.")
print("2. The mass of the rod and the bob is distributed uniformly.")
print("3. The motion is planar.")
print("4. The system is subject to a uniform gravitational field.")
print("5. The system is isolated and there is no friction or air resistance.")
print("=================================================================")

print("USE SI UNITS ONLY FOR ALL INPUTS")
# Parameters
length = float(input("Enter the length of pendulum's rod: "))
R = float(input("Enter the radius of the pendulum's rod: "))
m = float(input("Enter the mass of the pendulum's bob: "))
M = float(input("Enter the mass of the pendulum's rod: "))
r = float(input("Enter the radius of the pendulum's bob: "))
I = (1/3)*M*length**2 + m*(length+r)**2
g = 9.81
dt = 0.01
steps = int(input("Enter the total time for the simulation: "))*100

# Arrays
x = np.zeros(steps)
v = np.zeros(steps)
p = np.zeros(steps)
t = np.linspace(0, dt * (steps - 1), steps)
assert len(t) == len(x)

# Initial conditions
x[0] = float(input("Enter the initial angular position: "))
v[0] = float(input("Enter the initial angular velocity: "))
p[0] = I*v[0]

#Verlet Integration
for i in range(steps - 1):
    a = - (m*(length+r) + M*length/2)/I *g* np.sin(x[i])
    x[i + 1] = x[i] + v[i] * dt + 0.5 * a * dt**2
    a_next = - (m*(length+r) + M*length/2)/I * g * np.sin(x[i + 1])
    v[i + 1] = v[i] + 0.5 * (a + a_next) * dt
    p[i + 1] = I*v[i + 1]


#Plot of phase space
plt.figure()
plt.plot(x, p)
plt.xlabel("Angular Position (radians)")
plt.ylabel("Angular Momentum (kg m^2/s)")
plt.title("Phase Space Plot for Pendulum")
plt.show()

#Animation of the pendulum's motion using vpython
from vpython import *   
scene = canvas(title="Compound Pendulum Simulation", width=1000, height=1000) 
rod = cylinder(pos=vector(0, 0, 0), axis=vector(length, 0, 0), radius=R, color=color.blue)
bob = sphere(pos=vector(length + r, 0, 0), radius=r, color=color.green)
for i in range(steps):
    rate(100)
    rod.axis = vector(length * np.sin(x[i]), -length * np.cos(x[i]), 0)
    bob.pos = vector((length + r) * np.sin(x[i]), -(length + r) * np.cos(x[i]), 0)


