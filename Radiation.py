import numpy as np
import matplotlib.pyplot as plt
print("This program plots the temperature change of a substance over time due to heat transfer by radiation. It is based on the Stefan-Boltzmann Law.")
print("=================================================================")
print("Assumptions: ")
print("1. The heat transfer is only due to radiation.")
print("2. The ambient temperature is constant.")
print("3. The heat capacity is constant throughout the process.")
print("==================================================================")
print("USE SI UNITS ONLY FOR ALL INPUTS")
# Parameters and Initial Conditions
E = float(input("Enter the emissivity of the substance: "))
if E < 0 or E > 1:
    print("Emissivity must be between 0 and 1. Please enter a valid value.")
    exit()
T_s = float(input("Enter the constant ambient temperature: "))
A = float(input("Enter the effective surface area of the substance: "))
T_0 = float(input("Enter the initial temperature of the substance: "))
m = float(input("Enter the mass of the substance: "))
c = float(input("Enter the specific heat capacity of the substance: "))
steps = int(input("Enter the time limit: "))*100
dt = 0.01
sigma = 5.67e-8 
p = A*sigma*E

#Arrays
t = np.linspace(0, dt*(steps-1), steps)
T = np.zeros(steps)
T[0] = T_0

#RK4 Method
for i in range(1, steps):
    dT_dt =  p*(T_s**4-T[i-1]**4) / (m*c)
    k1 = dt*dT_dt
    dT_dt =  p*(T_s**4-(T[i-1]+k1/2)**4) / (m*c)
    k2 = dt*dT_dt
    dT_dt =  p*(T_s**4-(T[i-1]+k2/2)**4) / (m*c)
    k3 = dt*dT_dt
    dT_dt =  p*(T_s**4-(T[i-1]+k3)**4) / (m*c)
    k4 = dt*dT_dt
    T[i] = T[i-1] + (k1 + 2*k2 + 2*k3 + k4)/6

# Plot
plt.plot(t, T, label="Temperature vs Time", color='green')
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.title("Heat Transfer by Radiation")
plt.show()
