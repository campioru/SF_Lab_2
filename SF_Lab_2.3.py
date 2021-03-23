# Exercise 3


# This script compares the trapezoidal rule and the Runge-Kutta method by
# plotting theta for the nonlinear pendulum using each method


import numpy as np
import matplotlib.pyplot as plt

k = 0.0
phi = 0.66667
A = 0.0

dt = 0.01
timesteps = 10000

def f_nonlinear(theta_, omega_, t_):
    return - np.sin(theta_) - k * omega_ + A * np.cos(phi * t_)

def theta_plot_trapezoidal(theta_, omega_, t_):
    t_array = [t_]
    theta_array = [theta_]
    
    for i in range(timesteps):
        k1a = omega_ * dt
        k1b = f_nonlinear(theta_, omega_, t_) * dt
        k2a = (omega_ + k1b) * dt
        k2b = f_nonlinear(theta_ + k1a, omega_ + k1b, t_ + dt) * dt
        
        t_ += dt
        theta_ += (k1a + k2a) / 2
        omega_ += (k1b + k2b ) / 2
    
        t_array.append(t_)
        theta_array.append(theta_)

    plt.plot(t_array, theta_array, label = "Trapezoidal rule", color = "blue")

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

    plt.plot(t_array, theta_array, label = "Runge-Kutta algorithm", color = "orange")

theta_plot_trapezoidal(3.14, 0.0, 0.0)
theta_plot_rungekutta(3.14, 0.0, 0.0)
plt.xlabel("$t$, in s")
plt.ylabel(r"$\theta$, in rad")
plt.title(r"Plot of $\theta(t)$ for $\theta(0) = 3.14$, $\omega(0) = 0.0$")
plt.xlim(0.0, timesteps * dt)
plt.legend(loc = "best")
plt.savefig("Trapezoidal vs Runge-Kutta for theta(0) = 3.14, omega(0) = 0.0.pdf")
plt.show()