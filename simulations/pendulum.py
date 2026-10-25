import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation

# 1. Physical parameters
g = 9.81     # Gravity (m/s²)
L = 1.0      # Length of the pendulum (m)

# 2. Initial conditions
theta0 = np.pi / 4  # Initial angle: 45°
omega0 = 0.0        # Initial angular velocity: 0

# 3. Equations of the pendulum
# theta: angle (rad)
# theta' = omega: angular velocity (rad/s)
# theta'' = omega' = -(g / L) * sin(theta): angular acceleration (rad/s²)
def pendulum(t, y):
    theta, omega = y

    dtheta = omega
    domega = -(g / L) * np.sin(theta)

    return [dtheta, domega]

# 4. Simulation time
# Time interval for the simulation (start, end, seconds)
t = np.linspace(0, 10, 500)

# 5. Résolution numérique
solution = solve_ivp(
    pendulum, # Physical equations of the pendulum
    [0, 10], # Simulation time interval (start, end, seconds)
    [theta0, omega0], # Initial conditions (angle, angular velocity)
    t_eval=t # Time points at which to store the computed solution
)

# 6. Get angles from numerical simulation
theta = solution.y[0]

# 7. Convert angles to Cartesian coordinates
x = L * np.sin(theta)
y = -L * np.cos(theta)

# 8. Create figure
fig, ax = plt.subplots()
ax.set_xlim(-1.5 * L, 1.5 * L)
ax.set_ylim(-1.5 * L, 0.5 * L)
ax.set_aspect("equal") # Equal aspect ratio to ensure the pendulum looks correct
ax.grid()

# 9. Create pendulum graphics
line, = ax.plot([], [], "k-", linewidth=2)
bob, = ax.plot([], [], "bo", markersize=15)

# 10. Update pendulum for each frame
def update(frame):
    line.set_data([0, x[frame]], [0, y[frame]])
    bob.set_data([x[frame]], [y[frame]])
    return line, bob

# 11. Create animation
fps = 30  # Frames per second

animation = FuncAnimation(
    fig,
    update,
    frames=len(theta),
    interval=1000/fps, # Time between frames in milliseconds
    blit=True # Optimize drawing by only updating changed parts of the plot
)

# 12. Save animation as GIF
animation.save("videos/pendulum.gif", writer="pillow", fps=fps)

plt.close(fig)