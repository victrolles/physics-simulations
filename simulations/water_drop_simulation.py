"""2D falling water drop made of 1000 interacting particles.

Run: python water_drop_simulation.py
Controls: Space = pause/resume; R = reset.

Educational double-density-relaxation fluid solver inspired by Clavet et al. (2005).
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.spatial import cKDTree

# Screen/world parameters, measured in display units.
WIDTH, HEIGHT = 800, 600
FLOOR = 570
PARTICLE_SPACING = 6.5
INTERACTION_RADIUS = 18.0
GRAVITY = 0.60                       # pixels per simulation step squared
VISCOSITY = 0.018                    # velocity smoothing between neighbors
STIFFNESS = 0.22                     # pressure stiffness
NEAR_STIFFNESS = 0.32                # prevents clumping, allows splashes
SOLVER_ITERATIONS = 3
STEPS_PER_FRAME = 1


def initial_positions():
    rng = np.random.default_rng(7)
    radius = 117.0
    x = np.arange(-radius, radius + PARTICLE_SPACING, PARTICLE_SPACING)
    xx, yy = np.meshgrid(x, x)
    points = np.column_stack((xx.ravel(), yy.ravel()))
    points = points[np.sum(points**2, axis=1) < radius**2]
    points = points[rng.choice(len(points), size=1000, replace=False)]
    points += rng.uniform(-0.25, 0.25, points.shape)
    points += (WIDTH / 2, 165)
    return points

class Fluid:
    def __init__(self):
        self.pos = initial_positions()
        self.vel = np.zeros_like(self.pos)
        self.time = 0
        self.rest_density = self.find_rest_density()

    def neighbors(self):
        pairs = cKDTree(self.pos).query_pairs(INTERACTION_RADIUS, output_type="ndarray")
        if not len(pairs):
            return None
        i, j = pairs[:, 0], pairs[:, 1]
        delta = self.pos[j] - self.pos[i]
        dist = np.maximum(np.linalg.norm(delta, axis=1), 1e-6)
        q = 1.0 - dist / INTERACTION_RADIUS
        return i, j, delta, dist, q

    def find_rest_density(self):
        neighbors = self.neighbors()
        i, j, _, _, q = neighbors
        rho = np.zeros(len(self.pos))
        np.add.at(rho, i, q*q)
        np.add.at(rho, j, q*q)
        return float(np.percentile(rho, 60))

    def step(self):
        old_pos = self.pos.copy()

        # Smooth nearby particles' relative velocities (viscosity).
        neighbors = self.neighbors()
        if neighbors is not None:
            i, j, _, _, q = neighbors
            vdelta = (self.vel[j] - self.vel[i]) * (VISCOSITY * q)[:, None]
            vchange = np.zeros_like(self.vel)
            np.add.at(vchange, i, vdelta)
            np.add.at(vchange, j, -vdelta)
            self.vel += vchange

        # Gravity + predicted particle positions.
        self.vel[:, 1] += GRAVITY
        self.pos += self.vel

        # Pressure / near-pressure corrections make the points behave as fluid.
        for _ in range(SOLVER_ITERATIONS):
            neighbors = self.neighbors()
            if neighbors is None:
                break
            i, j, delta, dist, q = neighbors
            rho = np.zeros(len(self.pos))
            near_rho = np.zeros(len(self.pos))
            np.add.at(rho, i, q**2)
            np.add.at(rho, j, q**2)
            np.add.at(near_rho, i, q**3)
            np.add.at(near_rho, j, q**3)

            pressure = STIFFNESS * (rho - self.rest_density)
            near_pressure = NEAR_STIFFNESS * near_rho

            # Positive pressure pushes particles apart; negative pressure
            # provides weak cohesion and helps keep droplets together.
            amount = 0.5 * (
                (pressure[i] + pressure[j]) * q
                + (near_pressure[i] + near_pressure[j]) * q**2
            )
            amount = np.clip(amount, -2.0, 2.0)
            correction = (amount / dist)[:, None] * delta
            shift = np.zeros_like(self.pos)
            np.add.at(shift, i, -correction)
            np.add.at(shift, j, correction)
            self.pos += shift

            # Keep all particles inside the container.
            self.pos[:, 0] = np.clip(self.pos[:, 0], 8, WIDTH - 8)
            self.pos[:, 1] = np.clip(self.pos[:, 1], 8, FLOOR)

        # Reconstruct velocities after the pressure solver and collisions.
        self.vel = self.pos - old_pos
        grounded = self.pos[:, 1] >= FLOOR - 0.01
        self.vel[grounded, 0] *= 0.93
        self.vel[grounded, 1] = np.minimum(self.vel[grounded, 1], 0)
        self.time += 1

fluid = Fluid()

fig, ax = plt.subplots(figsize=(10, 7.5))
fig.patch.set_facecolor("#101827")
ax.set_facecolor("#17243b")
ax.set_xlim(0, WIDTH)
ax.set_ylim(HEIGHT, 0)  # Positive y points down, as in screen coordinates.
ax.set_aspect('equal')
ax.axis('off')
ax.axhline(FLOOR + 3, color="#d0d9e6", linewidth=3)
ax.vlines([8, WIDTH - 8], ymin=100, ymax=FLOOR + 3, color="#65758d", linewidth=1, alpha=0.35)
scat = ax.scatter(fluid.pos[:,0], fluid.pos[:,1], s=28, color="#39c8ff", edgecolor='none')
ax.text(18, 32, "2D water drop  |  %s particles" % len(fluid.pos), color="white", fontsize=14)
time_label = ax.text(18, 53, "", color="#b8c5d7", fontsize=10)
paused = False

def key_press(event):
    global fluid, paused
    if event.key == " ":
        paused = not paused
    elif event.key and event.key.lower() == "r":
        fluid = Fluid()
        paused = False
fig.canvas.mpl_connect("key_press_event", key_press)

def update(_frame):
    if not paused:
        for _ in range(STEPS_PER_FRAME):
            fluid.step()
    scat.set_offsets(fluid.pos)
    time_label.set_text(f"Simulation step: {fluid.time}   |   SPACE pause   R reset")
    return scat, time_label

fps = 30  # Frames per second

animation = FuncAnimation(
    fig,
    update,
    interval=30,
    blit=True,
    cache_frame_data=False
    )

animation.save("videos/water_drop_simulation.gif", writer="pillow", fps=fps)

plt.close(fig)