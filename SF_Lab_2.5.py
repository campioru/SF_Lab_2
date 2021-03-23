# Exercise 5

# This script plots the phase portraits of the damped, driven pendulum by
# calculating theta and omega for various values of the driving coefficient
# using the Runge-Kutta method and plotting them against each other


import numpy as np
import matplotlib.pyplot as plt

k = 0.5
phi = 0.66667
A = [0.90, 1.07, 1.35, 1.47, 1.5]
a = [0.0, 0.1, 0.0, 0.0, 0.0]
colours = [(62.0/255.0, 118.0/255.0, 236.0/255.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 206.0/255.0, 1.0/255.0), (23.0/255.0, 154.0/255.0, 19.0/255.0)]

dt = 0.01
timesteps = 25000
transient = 10000

def f_nonlinear(theta_, omega_, t_, A_):
    return - np.sin(theta_) - k * omega_ + A_ * np.cos(phi * t_)

def theta_omega_rungekutta(theta_, omega_, t_, A_, a_):
    t_array = [t_]
    theta_array = [theta_]
    omega_array = [omega_]
    
    for i in range(timesteps):
        k1a = omega_ * dt
        k1b = f_nonlinear(theta_, omega_, t_, A_) * dt
        k2a = (omega_ + k1b / 2) * dt
        k2b = f_nonlinear(theta_ + k1a / 2, omega_ + k1b / 2, t_ + dt / 2, A_) * dt
        k3a = (omega_ + k2b / 2) * dt
        k3b = f_nonlinear(theta_ + k2a / 2, omega_ + k2b / 2, t_ + dt / 2, A_) * dt
        k4a = (omega_ + k3b) * dt
        k4b = f_nonlinear(theta_ + k3a, omega_ + k3b, t_ + dt, A_) * dt
        
        theta_ += (k1a + 2 * k2a + 2 * k3a + k4a) / 6
        omega_ += (k1b + 2 * k2b + 2 * k3b + k4b) / 6
        t_ += dt
        
        if theta_ < -np.pi + a_ or theta_ > np.pi + a_:
            theta_ -= 2 * np.pi * np.abs(theta_) / theta_
    
        t_array.append(t_)
        theta_array.append(theta_)
        omega_array.append(omega_)
    
    return [theta_array, omega_array]

for i in range(5):
    theta = 1.0
    omega = 0.0
    t = 0.0
    [theta_array, omega_array] = theta_omega_rungekutta(theta, omega, t, A[i], a[i])
    plt.scatter(theta_array[transient:], omega_array[transient:], s = 0.2, color = colours[i])
    plt.xlabel(r"$\theta$, in rad")
    plt.ylabel(r"$\omega$, in rad s$^{-1}$")
    plt.title(r"Phase portrait for $A = %s$" % A[i])
    plt.xlim(-np.pi + a[i], np.pi + a[i])
    plt.savefig("Phase portrait for A = %s.pdf" % A[i])
    plt.show()