# Exercise 1


# This script plots theta and omega of the linear pendulum as functions of
# time using the trapezoidal rule


import numpy as np
import matplotlib.pyplot as plt

k = 0.0
phi = 0.66667
A = 0.0

theta = 0.2
omega = 0.0
t = 0.0
dt = 0.01
timesteps = 1000

def f_linear(theta_, omega_, t_):
    return - theta_ - k * omega_ + A * np.cos(phi * t_)

t_array = [t]
theta_array = [theta]
omega_array = [omega]
for i in range(timesteps):
    k1a = omega * dt
    k1b = f_linear(theta, omega, t) * dt
    k2a = (omega + k1b) * dt
    k2b = f_linear(theta + k1a, omega + k1b, t + dt) * dt
    
    t += dt
    theta += (k1a + k2a) / 2
    omega += (k1b + k2b ) / 2
    
    t_array.append(t)
    theta_array.append(theta)
    omega_array.append(omega)

plt.plot(t_array, theta_array, color = "blue")
plt.xlabel("$t$, in s")
plt.ylabel(r"$\theta$, in rad")
plt.title(r"Plot of $\theta(t)$ for $\theta(0) = 0.2$ and $\omega(0) = 0.0$")
plt.xlim(0.0, timesteps * dt)
plt.savefig("Linear theta(t) for theta(0) = 0.2 and omega(0) = 0.0.pdf")
plt.show()
plt.plot(t_array, omega_array, color = "red")
plt.xlabel("$t$, in s")
plt.ylabel(r"$\omega$, in rad s$^{-1}$")
plt.title(r"Plot of $\omega(t)$ for $\theta(0) = 0.2$ and $\omega(0) = 0.0$")
plt.xlim(0.0, timesteps * dt)
plt.savefig("Linear omega(t) for theta(0) = 0.2 and omega(0) = 0.0.pdf")
plt.show()

def theta_plot(theta_, omega_, t_):
    theta0 = theta_
    omega0 = omega_
    
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

    plt.plot(t_array, theta_array, label = r"$\theta(0) = %s$, $\omega(0) = %s$" % (theta0, omega0))

def omega_plot(theta_, omega_, t_):
    theta0 = theta_
    omega0 = omega_
    
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

    plt.plot(t_array, omega_array, label = r"$\theta(0) = %s$, $\omega(0) = %s$" % (theta0, omega0))

thetas = [0.2, 1.0, 3.14, 0.0]
omegas = [0.0, 0.0, 0.0, 1.0]

for i in range(4):
    theta_plot(thetas[i], omegas[i], 0.0)
plt.xlabel("$t$, in s")
plt.ylabel(r"$\theta$, in rad")
plt.title(r"Plot of $\theta(t)$ for various $\theta(0)$ and $\omega(0)$")
plt.xlim(0.0, timesteps * dt)
plt.legend(loc = "upper left", fontsize = "small")
plt.savefig("Linear theta(t).pdf")
plt.show()
for i in range(4):
    omega_plot(thetas[i], omegas[i], 0.0)
plt.xlabel("$t$, in s")
plt.ylabel(r"$\omega$, in rad s$^{-1}$")
plt.title(r"Plot of $\omega(t)$ for various $\theta(0)$ and $\omega(0)$")
plt.xlim(0.0, timesteps * dt)
plt.legend(loc = "upper left", fontsize = "small")
plt.savefig("Linear omega(t).pdf")
plt.show()