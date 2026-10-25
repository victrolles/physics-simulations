import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# 1. Parameters
N = 50                  # Grid size
alpha = 1.0             # Thermal diffusivity (arbitrary units)
dx = 1.0                # Spatial step
dt = 0.1                # Time step

# 2. Initial temperature
T = np.full((N, N), 20.0)
T[20:30, 20:30] = 100.0  # Hot center

# 3. Create figure
fig, ax = plt.subplots()

im = ax.imshow(T, cmap="hot", vmin=20, vmax=100)
plt.colorbar(im, ax=ax, label="Temperature (°C)")
ax.set_title("Heat diffusion")

# 4. Update temperature
def update(frame):
    global T

    T_new = T.copy()

    T_new[1:-1, 1:-1] = T[1:-1, 1:-1] + alpha * dt / dx**2 * (
        T[2:, 1:-1] +
        T[:-2, 1:-1] +
        T[1:-1, 2:] +
        T[1:-1, :-2] -
        4 * T[1:-1, 1:-1]
    )

    T = T_new
    im.set_data(T)

    return (im,)

# 5. Create animation
fps = 30  # Frames per second

animation = FuncAnimation(
    fig,
    update,
    frames=300,
    interval=1000/fps,
    blit=True
)

# 6. Save animation as GIF
animation.save("videos/heat_diffusion.gif", writer="pillow", fps=fps)

plt.close(fig)