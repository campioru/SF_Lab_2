# Exercise 4


# This script plots theta and omega for the damped, nondriven pendulum using
# the Runge-Kutta method


import numpy as np
import matplotlib.pyplot as plt

k = 0.5
phi = 0.66667
A = 0.0

dt = 0.01
timesteps = 3000

def f_nonlinear(theta_, omega_, t_):
    return - np.sin(theta_) - k * omega_ + A * np.cos(phi * t_)

def theta_plot_rungekutta(theta_, omega_, t_):
    t_array = [t_]
    theta_array = [theta_]
    
    for i in range(timesteps):
        k1a = omega_ * dt
        k1b = f_nonlinear(theta_, omega_, t_) * dt
        k2a = (omega_ + k1b / 2) * dt
        k2b = f_nonlinear(theta_ + k1a / 2, omega_ + k1b / 2, t_ + dt / 2) * dt
        k3a = (omega_ + k2b / 2) * dt
        k3b = f_nonlinear(theta_ + k2a / 2, omega_ + k2b / 2, t_ + dt / 2) * dt
        k4a = (omega_ + k3b) * dt
        k4b = f_nonlinear(theta_ + k3a, omega_ + k3b, t_ + dt) * dt
        
        theta_ += (k1a + 2 * k2a + 2 * k3a + k4a) / 6
        omega_ += (k1b + 2 * k2b + 2 * k3b + k4b) / 6
        t_ += dt
    
        t_array.append(t_)
        theta_array.append(theta_)

    plt.plot(t_array, theta_array, color = "blue")

def omega_plot_rungekutta(theta_, omega_, t_):
    t_array = [t_]
    omega_array = [omega_]
    
    for i in range(timesteps):
        k1a = omega_ * dt
        k1b = f_nonlinear(theta_, omega_, t_) * dt
        k2a = (omega_ + k1b / 2) * dt
        k2b = f_nonlinear(theta_ + k1a / 2, omega_ + k1b / 2, t_ + dt / 2) * dt
        k3a = (omega_ + k2b / 2) * dt
        k3b = f_nonlinear(theta_ + k2a / 2, omega_ + k2b / 2, t_ + dt / 2) * dt
        k4a = (omega_ + k3b) * dt
        k4b = f_nonlinear(theta_ + k3a, omega_ + k3b, t_ + dt) * dt
        
        theta_ += (k1a + 2 * k2a + 2 * k3a + k4a) / 6
        omega_ += (k1b + 2 * k2b + 2 * k3b + k4b) / 6
        t_ += dt
    
        t_array.append(t_)
        omega_array.append(omega_)

    plt.plot(t_array, omega_array, color = "red")

theta_plot_rungekutta(3.0, 0.0, 0.0)
plt.xlabel("$t$, in s")
plt.ylabel(r"$\theta$, in rad")
plt.title(r"Plot of $\theta(t)$ for $\theta(0) = 3.0$, $\omega(0) = 0.0$")
plt.xlim(0.0, timesteps * dt)
plt.ylim(-1.0, 3.1)
plt.savefig("Damped theta(t) for theta(0) = 3.0, omega(0) = 0.0.pdf")
plt.show()
omega_plot_rungekutta(3.0, 0.0, 0.0)
plt.xlabel("$t$, in s")
plt.ylabel(r"$\omega$, in rad s$^{-1}$")
plt.title(r"Plot of $\omega(t)$ for $\theta(0) = 3.0$, $\omega(0) = 0.0$")
plt.xlim(0.0, timesteps * dt)
plt.savefig("Damped omega(t) for theta(0) = 3.0, omega(0) = 0.0.pdf")
plt.show()