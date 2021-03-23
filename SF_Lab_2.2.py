# Exercise 2


# This script compares the linear and nonlinear pendulum by plotting theta
# and omega for each pendulum, using the trapezoidal rule


import numpy as np
import matplotlib.pyplot as plt

k = 0.0
phi = 0.66667
A = 0.0

dt = 0.01
timesteps = 3500

def f_linear(theta_, omega_, t_):
    return - theta_ - k * omega_ + A * np.cos(phi * t_)

def f_nonlinear(theta_, omega_, t_):
    return - np.sin(theta_) - k * omega_ + A * np.cos(phi * t_)

def theta_plot_linear(theta_, omega_, t_):
    t_array = [t_]
    theta_array = [theta_]
    
    for i in range(timesteps):
        k1a = omega_ * dt
        k1b = f_linear(theta_, omega_, t_) * dt
        k2a = (omega_ + k1b) * dt
        k2b = f_linear(theta_ + k1a, omega_ + k1b, t_ + dt) * dt
        
        t_ += dt
        theta_ += (k1a + k2a) / 2
        omega_ += (k1b + k2b ) / 2
    
        t_array.append(t_)
        theta_array.append(theta_)

    plt.plot(t_array, theta_array, label = r"Linear $\theta(t)$", color = "red")

def omega_plot_linear(theta_, omega_, t_):
    t_array = [t_]
    omega_array = [omega_]
    
    for i in range(timesteps):
        k1a = omega_ * dt
        k1b = f_linear(theta_, omega_, t_) * dt
        k2a = (omega_ + k1b) * dt
        k2b = f_linear(theta_ + k1a, omega_ + k1b, t_ + dt) * dt
        
        t_ += dt
        theta_ += (k1a + k2a) / 2
        omega_ += (k1b + k2b ) / 2
    
        t_array.append(t_)
        omega_array.append(omega_)

    plt.plot(t_array, omega_array, label = r"Linear $\omega(t)$", color = "black")

def theta_plot_nonlinear(theta_, omega_, t_):
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

    plt.plot(t_array, theta_array, label = r"Non-linear $\theta(t)$", color = "blue")

def omega_plot_nonlinear(theta_, omega_, t_):
    t_array = [t_]
    omega_array = [omega_]
    
    for i in range(timesteps):
        k1a = omega_ * dt
        k1b = f_nonlinear(theta_, omega_, t_) * dt
        k2a = (omega_ + k1b) * dt
        k2b = f_nonlinear(theta_ + k1a, omega_ + k1b, t_ + dt) * dt
        
        t_ += dt
        theta_ += (k1a + k2a) / 2
        omega_ += (k1b + k2b ) / 2
    
        t_array.append(t_)
        omega_array.append(omega_)

    plt.plot(t_array, omega_array, label = r"Non-linear $\omega(t)$", color = "purple")

thetas = [0.2, 1.0, 3.14, 0.0]
omegas = [0.0, 0.0, 0.0, 1.0]

for i in range(4):
    theta_plot_linear(thetas[i], omegas[i], 0.0)
    theta_plot_nonlinear(thetas[i], omegas[i], 0.0)
    plt.xlabel("$t$, in s")
    plt.ylabel(r"$\theta$, in rad")
    plt.title(r"Plot of $\theta(t)$ for $\theta(0) = %s$, $\omega(0) = %s$" % (thetas[i], omegas[i]))
    plt.xlim(0.0, timesteps * dt)
    plt.legend(loc = "upper right", fontsize = "small")
    plt.savefig("Linear and non-linear theta(t) for theta(0) = %s and omega(0) = %s.pdf" % (thetas[i], omegas[i]))
    plt.show()
    omega_plot_linear(thetas[i], omegas[i], 0.0)
    omega_plot_nonlinear(thetas[i], omegas[i], 0.0)
    plt.xlabel("$t$, in s")
    plt.ylabel(r"$\omega$, in rad s$^{-1}$")
    plt.title(r"Plot of $\omega(t)$ for $\theta(0) = %s$, $\omega(0) = %s$" % (thetas[i], omegas[i]))
    plt.xlim(0.0, timesteps * dt)
    plt.legend(loc = "upper right", fontsize = "small")
    plt.savefig("Linear and non-linear omega(t) for theta(0) = %s and omega(0) = %s.pdf" % (thetas[i], omegas[i]))
    plt.show()