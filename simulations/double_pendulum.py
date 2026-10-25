import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation

# 1. Physical parameters
g = 9.81     # Gravity (m/s²)
L1 = 1.0      # Length of the pendulum (m)
L2 = 1.0      # Length of the second pendulum (m)
L = L1 + L2  # Total length of the double pendulum (m)
M1 = 1.0      # Mass of the first pendulum (kg)
M2 = 1.0      # Mass of the second pendulum (kg)

# 2. Initial conditions
theta1_0 = np.radians(120)  # Initial angle: 120°
theta2_0 = np.radians(-10)   # Initial angle: -10°
omega1_0 = 0.0        # Initial angular velocity: 0
omega2_0 = 0.0        # Initial angular velocity: 0

# 3. Equations of the pendulum
def double_pendulum(t, y):
    # State variables
    theta1, omega1, theta2, omega2 = y

    # Difference between the angles
    delta = theta1 - theta2

    # Coupled acceleration coefficients
    A = np.array([
        [(M1 + M2) * L1, M2 * L2 * np.cos(delta)],
        [L1 * np.cos(delta), L2]
    ])

    # Gravity and velocity-dependent terms
    B = np.array([
        -M2 * L2 * np.sin(delta) * omega2**2
        - (M1 + M2) * g * np.sin(theta1),

        L1 * np.sin(delta) * omega1**2
        - g * np.sin(theta2)
    ])

    # Solve for both angular accelerations
    alpha1, alpha2 = np.linalg.solve(A, B)

    # Return derivatives of the state
    return [omega1, alpha1, omega2, alpha2]

# 4. Simulation time
# Time interval for the simulation (start, end, seconds)
t = np.linspace(0, 10, 500)

# 5. Résolution numérique
solution = solve_ivp(
    double_pendulum,
    [0, 15],
    [theta1_0, omega1_0, theta2_0, omega2_0],
    t_eval=t,
    rtol=1e-9,
    atol=1e-9
)

# 6. Get angles from numerical simulation
theta1 = solution.y[0]
theta2 = solution.y[2]

# 7. Convert angles to Cartesian coordinates
x1 = L1 * np.sin(theta1)
y1 = -L1 * np.cos(theta1)
x2 = x1 + L2 * np.sin(theta2)
y2 = y1 - L2 * np.cos(theta2)

# 8. Create figure
fig, ax = plt.subplots()
ax.set_xlim(-1.2 * L, 1.2 * L)
ax.set_ylim(-1.2 * L, 1.2 * L)
ax.set_aspect("equal") # Equal aspect ratio to ensure the pendulum looks correct
ax.grid()

# 9. Create pendulum graphics
line1, = ax.plot([], [], "k-", linewidth=2)
line2, = ax.plot([], [], "k-", linewidth=2)
bob1, = ax.plot([], [], "bo", markersize=15)
bob2, = ax.plot([], [], "ro", markersize=15)

# 10. Update pendulum for each frame
def update(frame):
    line1.set_data([0, x1[frame]], [0, y1[frame]])
    line2.set_data([x1[frame], x2[frame]], [y1[frame], y2[frame]])
    bob1.set_data([x1[frame]], [y1[frame]])
    bob2.set_data([x2[frame]], [y2[frame]])
    return line1, line2, bob1, bob2

# 11. Create animation
fps = 30  # Frames per second

animation = FuncAnimation(
    fig,
    update,
    frames=len(theta1),
    interval=1000/fps, # Time between frames in milliseconds
    blit=True # Optimize drawing by only updating changed parts of the plot
)

# 12. Save animation as GIF
animation.save("videos/double_pendulum.gif", writer="pillow", fps=fps)

plt.close(fig)