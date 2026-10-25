# Water Drop: Simulating 1,000 Interacting Particles

[← Back to README](../README.md) · [View Python source](../simulations/water_drop_simulation.py)

The previous simulations track a few degrees of freedom (pendulum angles) or a regular grid (temperature). A falling water drop is much more complicated: it **deforms, flows, collides, and splashes**. Here we approximate that behavior with **1,000 interacting particles**, rather than solving the full Navier–Stokes equations.

This tutorial walks through the actual implementation in `simulations/water_drop_simulation.py`. The algorithm is a simplified **double-density-relaxation (DDR)** or position-based particle method inspired by Clavet, Beaudoin, and Poulin, *Particle-based Viscoelastic Fluid Simulation* (2005). It is **qualitative and educational**, not a quantitatively validated fluid solver.

## 1. What are we actually simulating?

Each particle $i$ has a 2D position and velocity:

$$
\mathbf x_i=\begin{bmatrix}x_i\\y_i\end{bmatrix},
\qquad
\mathbf v_i=\begin{bmatrix}v_{x,i}\\v_{y,i}\end{bmatrix}.
$$

The drop is a **cloud of interacting points**, with no explicit surface mesh. The visible blue markers give the visual impression of liquid. This model does not calculate surface tension from real material constants or solve the incompressible Navier–Stokes PDE on a mesh.

The display is $800\times600$ coordinate units. Unlike the pendulum drawings, **positive $y$ points downward**. The floor is at $y=570$, so positive gravitational acceleration increases $v_y$.

| Constant | Role | Value |
| --- | --- | --- |
| `WIDTH`, `HEIGHT` | Canvas dimensions | 800, 600 |
| `FLOOR` | Bottom collision boundary | 570 |
| `PARTICLE_SPACING` | Spacing of initial candidate points | 6.5 |
| `INTERACTION_RADIUS` | Distance over which particles interact | 18.0 |
| `GRAVITY` | Downward acceleration per step | 0.60 |
| `VISCOSITY` | Local velocity-smoothing coefficient | 0.018 |
| `STIFFNESS` | Density-dependent pressure coefficient | 0.22 |
| `NEAR_STIFFNESS` | Short-range repulsion coefficient | 0.32 |
| `SOLVER_ITERATIONS` | Position-relaxation iterations per step | 3 |

These are **simulation/display parameters** (approximately pixels and steps), not calibrated SI units like meters, seconds, pascals, or kilograms. The implicit simulation time increment is **one step**.

## 2. Create the initial circular drop

We want a disk-shaped drop composed of exactly 1,000 samples. Rather than drawing a single circle, we:

1. Generate a square lattice of candidate particle positions.
2. Retain only points inside a circle of radius $R=117$:

   $$x^2+y^2<R^2.$$

3. Randomly select **1,000** of the interior points without replacement.
4. Add small positional noise so the arrangement is not perfectly regular.
5. Translate the disk so its center is at $(400,165)$ on the canvas.

```python
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
```

The seed `7` makes the random selection **reproducible**. Each initial velocity is zero:

```python
self.pos = initial_positions()
self.vel = np.zeros_like(self.pos)
```

The particles initially form a cohesive circular cloud above the floor.

## 3. Why does each particle need neighbors?

If every particle interacted with all others, there would be roughly

$$\frac{N(N-1)}{2}$$

possible distinct pairs per update: with $N=1000$, that's **499,500 pairs**. Most are far apart and cannot affect one another under this local model.

Instead, restrict interactions to pairs with separation **less than $h=18$**, where $h$ is `INTERACTION_RADIUS`.

SciPy's `cKDTree` is a spatial search data structure that efficiently identifies nearby particles:

```python
pairs = cKDTree(self.pos).query_pairs(
    INTERACTION_RADIUS, output_type="ndarray"
)
```

Each returned pair $(i,j)$ appears once. Define the separation vector from $i$ to $j$ and its distance:

$$
\mathbf d_{ij}=\mathbf x_j-\mathbf x_i,
\qquad
d_{ij}=\|\mathbf d_{ij}\|.
$$

To avoid division by zero when particles nearly overlap, the code computes

$$d_{ij}\leftarrow\max(d_{ij},10^{-6}).$$

Then calculate the **proximity weight**:

$$
\boxed{q_{ij}=1-\frac{d_{ij}}{h},\qquad d_{ij}<h.}
$$

For particles very close together, $q$ is near 1. Near the interaction cutoff, $q$ approaches 0. Pairs with $d\geq h$ are not included.

This simple compact-support kernel lets us approximate local density and near-field interactions.

## 4. Determine the reference or rest density

Particles close together should exert pressure that prevents excessive compression. First compute a local density proxy:

$$
\boxed{\rho_i=\sum_{j\in\mathcal N(i)}q_{ij}^2.}
$$

Here $\mathcal N(i)$ denotes all neighbors of particle $i$ inside the interaction radius. These are **dimensionless density estimates**, not physical density in kg/m³.

The code accumulates each pair's contribution on *both* participating particles:

```python
rho = np.zeros(len(self.pos))
np.add.at(rho, i, q*q)
np.add.at(rho, j, q*q)
```

`np.add.at` handles repeated particle indices correctly, accumulating multiple neighbors rather than overwriting their contributions.

The **reference density** is calculated once using the initial particle arrangement:

```python
self.rest_density = float(np.percentile(rho, 60))
```

So $\rho_0$ is the **60th percentile** of the initial density distribution. This is a **heuristic** to capture a typical interior particle density without hand-tuning a constant. It is not a measured water density.

## 5. Smooth neighboring velocities (viscosity)

Real liquids resist rapid local differences in velocity. In the code, viscosity is approximated by directly smoothing neighbor velocities.

For each nearby pair $(i,j)$, define:

$$
\boxed{\Delta\mathbf v_{ij}
=\nu q_{ij}(\mathbf v_j-\mathbf v_i),}
$$

where $\nu=0.018$ is the numerical viscosity parameter. Update particle $i$ by $+\Delta\mathbf v_{ij}$ and particle $j$ by $-\Delta\mathbf v_{ij}$. Equal and opposite changes tend to bring their velocities closer together.

```python
vdelta = (self.vel[j] - self.vel[i]) * (VISCOSITY * q)[:, None]
vchange = np.zeros_like(self.vel)
np.add.at(vchange, i, vdelta)
np.add.at(vchange, j, -vdelta)
self.vel += vchange
```

This is **not** a calibrated solution of the viscous stress term from Navier–Stokes. It is a simple, symmetric local velocity correction. The `[:, None]` converts a vector of scalar pair weights into an array broadcastable over the two coordinates $(x,y)$.

## 6. Apply gravity and predict new positions

Gravity acts on the vertical velocity only. The simulation uses a unit time increment $\Delta t=1$:

$$
\boxed{v_{y,i}^{*}=v_{y,i}+g_{\mathrm{sim}},\quad g_{\mathrm{sim}}=0.60.}
$$

Then predict the next particle position:

$$
\boxed{\mathbf x_i^{*}=\mathbf x_i+\mathbf v_i^{*}.}
$$

```python
self.vel[:, 1] += GRAVITY
self.pos += self.vel
```

This is a **semi-implicit Euler-like prediction step**: update velocity first, then use the updated velocity for position. It makes the drop fall toward the bottom of the screen.

The code saves `old_pos` before changing any velocities or positions so that it can **reconstruct a velocity later** after enforcing the density and collision constraints.

## 7. Compute density and *near density* after prediction

Once positions have moved, neighbors change. Rebuild the spatial neighbor pairs, then compute two density proxies:

$$
\boxed{
\rho_i=\sum_{j\in\mathcal N(i)}q_{ij}^{2},
\qquad
\rho_i^{\mathrm{near}}=\sum_{j\in\mathcal N(i)}q_{ij}^{3}.
}
$$

`near_rho` uses $q^3$, so it emphasizes **shorter separations** more strongly than the ordinary density proxy.

```python
rho = np.zeros(len(self.pos))
near_rho = np.zeros(len(self.pos))

np.add.at(rho, i, q**2)
np.add.at(rho, j, q**2)
np.add.at(near_rho, i, q**3)
np.add.at(near_rho, j, q**3)
```

Because particles can move during a relaxation iteration, the code rebuilds the neighbor search at **each iteration**, not just once per frame.

## 8. Convert densities into pressure-like terms

The code uses two empirical pressure laws:

$$
\boxed{P_i=k(\rho_i-\rho_0),}
$$

$$
\boxed{P_i^{\mathrm{near}}=k_{\mathrm{near}}\rho_i^{\mathrm{near}}.}
$$

Here $k=0.22$ and $k_{\mathrm{near}}=0.32$ are `STIFFNESS` and `NEAR_STIFFNESS`. They control how strongly density differences influence particle separation.

**Ordinary pressure** becomes positive when a particle is too crowded ($\rho_i>\rho_0$), pushing it away from neighbors. It can become negative when local density falls below rest density, causing **weak cohesion** in this implementation.

**Near pressure** is always nonnegative and becomes relatively important at short distances, resisting clumping.

```python
pressure = STIFFNESS * (rho - self.rest_density)
near_pressure = NEAR_STIFFNESS * near_rho
```

These values are dimensionless **correction coefficients**, not pressures measured in pascals.

## 9. Relax positions using neighboring pressures

For a pair $(i,j)$, calculate a scalar displacement magnitude based on both particles' pressure estimates:

$$
\boxed{
a_{ij}=\frac12\left[
(P_i+P_j)q_{ij}
+(P_i^{\mathrm{near}}+P_j^{\mathrm{near}})q_{ij}^{2}
\right].
}
$$

Limit the magnitude for robustness:

$$\hat a_{ij}=\operatorname{clip}(a_{ij},-2,2).$$

Turn it into a vector pointing along the line between the particles:

$$
\boxed{\mathbf c_{ij}=\hat a_{ij}\frac{\mathbf d_{ij}}{d_{ij}}.}
$$

Then shift the pair in opposite directions:

$$
\boxed{
\mathbf x_i\leftarrow\mathbf x_i-\mathbf c_{ij},\qquad
\mathbf x_j\leftarrow\mathbf x_j+\mathbf c_{ij}.
}
$$

When $\hat a_{ij}>0$, the pair moves **apart**; when it is negative, the pair moves **together**. This mechanism gives the particle cloud a fluid-like degree of cohesion while preventing too much overlap.

In code:

```python
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
```

All corrections are **accumulated in `shift`** and applied together at the end of the iteration. This avoids updating pair positions sequentially based on an inconsistent mixture of old and new states.

Why use **three relaxation iterations**? One pass cannot necessarily resolve a strongly compressed cluster. Repeating pressure correction and neighbor search helps the drop maintain a more coherent distribution. The code uses:

```python
for _ in range(SOLVER_ITERATIONS):  # SOLVER_ITERATIONS = 3
    # Rebuild neighbors, densities and pressures
    # Accumulate pair corrections and update positions
```

## 10. Enforce walls and floor collisions

Without boundaries, particles could pass through the floor or leave the display area. After every pressure relaxation iteration, positions are clipped:

$$
8\leq x_i\leq 792,\qquad 8\leq y_i\leq 570.
$$

```python
self.pos[:, 0] = np.clip(self.pos[:, 0], 8, WIDTH - 8)
self.pos[:, 1] = np.clip(self.pos[:, 1], 8, FLOOR)
```

The floor is a **hard geometric constraint**, not a simulated deformable solid. Pressure corrections plus this constraint cause the particles to spread sideways on impact, producing a splash-like effect.

## 11. Reconstruct velocities after the constraints

The pressure solver and collision handling have **changed positions** after the initial gravity step. To make future motion reflect those changes, reconstruct each particle's velocity from its full-step displacement:

$$
\boxed{\mathbf v_i^{n+1}=\frac{\mathbf x_i^{n+1}-\mathbf x_i^n}{\Delta t}.}
$$

Because the simulation uses $\Delta t=1$, the code is simply:

```python
self.vel = self.pos - old_pos
```

Then handle particles touching the floor. Horizontal velocity is reduced to model friction-like energy loss, and downward motion into the floor is removed:

```python
grounded = self.pos[:, 1] >= FLOOR - 0.01
self.vel[grounded, 0] *= 0.93
self.vel[grounded, 1] = np.minimum(self.vel[grounded, 1], 0)
```

**Sign convention:** positive $y$ is down. Therefore `min(v_y, 0)` prevents a particle on the floor from retaining **positive downward velocity**. Negative $v_y$ (upward motion) remains possible.

The factor `0.93` is an empirical **tangential damping** coefficient, not a friction coefficient derived from a material law.

## 12. See the full algorithm at a glance

Each animation update calls `fluid.step()`, whose order matters:

```text
Initialize 1,000 positions and velocities
Compute initial reference density

For each simulation step:
    Save old particle positions
    Find nearby particle pairs
    Smooth velocity differences (viscosity)
    Apply gravity to velocities
    Predict new particle positions

    Repeat 3 times:
        Find nearby pairs at current positions
        Compute density and near density
        Compute pressure and near pressure
        Accumulate and apply position corrections
        Clamp positions inside the container

    Reconstruct velocities from position changes
    Dampen horizontal speed at the floor
    Remove downward velocity through the floor
    Advance the step counter
```

The resulting algorithm is an example of **position-based fluid dynamics**. It emphasizes visual plausibility and computational simplicity over exact conservation laws or experimental accuracy.

## 13. Animate and save the water drop

Matplotlib's scatter plot draws 1,000 blue particles at once:

```python
scat = ax.scatter(
    fluid.pos[:, 0], fluid.pos[:, 1],
    s=28, color="#39c8ff", edgecolor="none"
)
```

Every frame updates the scatter offsets (all particles' positions):

```python
def update(_frame):
    if not paused:
        for _ in range(STEPS_PER_FRAME):
            fluid.step()
    scat.set_offsets(fluid.pos)
    time_label.set_text(f"Simulation step: {fluid.time}   |   SPACE pause   R reset")
    return scat, time_label
```

The export is **finite** and reproducible:

```python
fps = 30
animation = FuncAnimation(
    fig, update,
    frames=180,
    interval=1000 / fps,
    blit=True,
    cache_frame_data=False
)
animation.save("videos/water_drop_simulation.gif", writer="pillow", fps=fps)
```

180 frames at 30 FPS produce approximately **6 seconds of GIF playback**, with one simulation step per animation frame (apart from the small number of extra initializer calls Matplotlib may perform). There is no calibrated mapping from simulation steps to real seconds.

**Important:** The code registers keyboard handlers for space (pause) and R (reset), but its default behavior is **rendering a GIF and closing the figure**, not opening an interactive window. Those keyboard controls are not usable on the exported GIF.

## 14. Accuracy, cost, and future improvements

- **Neighbor search:** A spatial tree avoids checking all $O(N^2)$ particle pairs, although very dense clusters can still create many interactions.
- **No physical mass/pressure calibration:** Particle density and pressure are proxies in arbitrary units.
- **No explicit Navier–Stokes solve:** Real water involves momentum conservation, pressure constraints, viscosity, a free surface, and possibly surface tension. This model approximates some visible effects with local corrections.
- **Simplified boundaries:** A hard flat floor and rectangular walls; no true fluid–solid contact model.
- **Useful experiments:** Vary `GRAVITY`, `INTERACTION_RADIUS`, `VISCOSITY`, `STIFFNESS`, `NEAR_STIFFNESS`, or `SOLVER_ITERATIONS` one at a time. Compare splash width, drop cohesion, run time, and particle penetration near boundaries.
- **Future challenge:** Compare this solver to a conventional **smoothed particle hydrodynamics (SPH)** implementation or a grid-based Navier–Stokes solver before training an ML surrogate on its outputs.

**Run from the repository root:**

```bash
python simulations/water_drop_simulation.py
```

The GIF will be written to `videos/water_drop_simulation.gif`.

## Further reading

Clavet, S., Beaudoin, P., and Poulin, P. (2005). *Particle-based Viscoelastic Fluid Simulation*. This project borrows the idea of **double-density relaxation** but implements its own simplified set of coefficients, corrections, and boundaries.
