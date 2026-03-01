import numpy as np
import matplotlib.pyplot as plt

print("This program plots the phase space of a damped harmonic oscillator using the Runge-Kutta 4th order method.")
print ("=================================================================================")
print("Assumptions: ")
print("1. The system is a simple harmonic oscillator subjected to a damping force proportional to velocity.")
print("2. The natural angular frequency and the mass is constant and does not change with time.")
print("3. The damping coefficient is constant and does not change with time.")
print("4. The oscillator is isolated and there are no external forces acting on it except the damping force.")
print("==================================================================================")
print("USE SI UNITS ONLY FOR ALL INPUTS")
# Parameters
omega = float(input("Enter the natural angular frequency: "))
b = float(input("Enter the damping coefficient: "))
m = float(input("Enter the mass: "))
gamma = b / m
dt = 0.01
steps = int(input("Enter the number of time steps: "))*100

# Arrays
x = np.zeros(steps)
v = np.zeros(steps)
p = np.zeros(steps)  
t = np.linspace(0, dt * (steps - 1), steps)
assert len(t) == len(x)


# Initial conditions
x[0] = float(input("Enter the initial position: "))
v[0] = float(input("Enter the initial velocity: "))
p[0] = m * v[0]

# RK4 Integration Method
for i in range(steps - 1):
    k1v = - (omega**2 * x[i] + gamma * v[i]) * dt
    k1x = v[i] * dt

    k2v = - (omega**2 * (x[i] + k1x/2) + gamma * (v[i] + k1v/2)) * dt
    k2x = (v[i] + k1v/2) * dt

    k3v = - (omega**2 * (x[i] + k2x/2) + gamma * (v[i] + k2v/2)) * dt
    k3x = (v[i] + k2v/2) * dt

    k4v = - (omega**2 * (x[i] + k3x) + gamma * (v[i] + k3v)) * dt
    k4x = (v[i] + k3v) * dt

    v[i+1] = v[i] + (k1v + 2*k2v + 2*k3v + k4v)/6
    x[i+1] = x[i] + (k1x + 2*k2x + 2*k3x + k4x)/6
    p[i+1] = m * v[i+1]  


#Plot of position vs time graph
import matplotlib.pyplot as pt
pt.figure()
pt.plot(t, x)
pt.xlabel("Time (s)")
pt.ylabel("Position (m)")
pt.title("Position vs Time Graph for Damped Harmonic Oscillator")
pt.show()

#Plot of velocity vs time graph
pt.figure()
pt.plot(t, v)
pt.xlabel("Time (s)")
pt.ylabel("Velocity (m/s)")
pt.title("Velocity vs Time Graph for Damped Harmonic Oscillator")
pt.show()

#Plot of phase space graph
pt.figure() 
pt.plot(x, p)
pt.xlabel("Position (m)")
pt.ylabel("Momentum (Ns)")
pt.title("Phase Space Graph for Damped Harmonic Oscillator")
pt.show()   
