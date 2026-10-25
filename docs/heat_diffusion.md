# Heat Diffusion: From Fourier's Law to a 2D Finite-Difference Simulation

[← Back to README](../README.md) · [View Python source](../simulations/heat_diffusion.py)

Unlike the pendulums, **temperature depends on both space and time**. Its governing equation is therefore a **partial differential equation (PDE)**. This guide explains where the equation comes from, how to discretize it, how the code updates a 2D grid, and why the time step must be small enough for numerical stability.

## 1. Define the physical problem

Imagine a thin, flat plate initially at $20^\circ$C, except for a hot square at $100^\circ$C near the center. Heat flows from hotter regions to colder ones.

The unknown is the temperature **field**:

$$T=T(x,y,t).$$

Three independent variables appear: the horizontal position $x$, the vertical position $y$, and the time $t$.

The script chooses:

| Parameter | Meaning | Value |
| --- | --- | --- |
| $N$ | Number of grid cells on each side | $50$ |
| $\alpha$ | Thermal diffusivity | $1.0$ in arbitrary units |
| $\Delta x$ | Spatial grid spacing (same in both directions) | $1.0$ |
| $\Delta t$ | Simulation time step | $0.1$ |
| $T_{\mathrm{ambient}}$ | Initial background temperature | $20^\circ$C |
| $T_{\mathrm{hot}}$ | Initial central temperature | $100^\circ$C |

This is a **2D model**. No thickness-dependent effects or heat sources are included. The chosen units are educational, not calibrated to a particular real metal or material.

## 2. Why does heat spread? Fourier's law

The heat flux vector $\mathbf q$ represents the rate of heat transfer through unit area. **Fourier's law** states:

$$
\boxed{\mathbf q=-k\nabla T,}
$$

where $k$ is thermal conductivity and the temperature gradient in 2D is:

$$
\nabla T=
\begin{bmatrix}
\partial T/\partial x\\
\partial T/\partial y
\end{bmatrix}.
$$

The minus sign means heat flows **down the temperature gradient**: from high temperature toward low temperature.

The gradient is zero in a uniform-temperature region, so there is no net conductive heat flux there.

## 3. Apply conservation of thermal energy

For a material with density $\rho$ and specific heat capacity $c$, the stored thermal energy per unit volume changes at a rate proportional to $\rho c\,\partial T/\partial t$.

With no external heat generation, local energy conservation states:

$$
\rho c\frac{\partial T}{\partial t}=-\nabla\cdot\mathbf q.
$$

Substitute Fourier's law:

$$
\rho c\frac{\partial T}{\partial t}
=\nabla\cdot(k\nabla T).
$$

Assume $k$, $\rho$, and $c$ are constant in space. Then:

$$
\rho c\frac{\partial T}{\partial t}
=k\nabla^2 T,
$$

where $\nabla^2$ is the **Laplacian**. Define thermal diffusivity:

$$\alpha=\frac{k}{\rho c}.$$

Its SI units would be $\mathrm{m^2/s}$, although the example uses arbitrary simulation units. Dividing by $\rho c$ gives:

$$
\boxed{\frac{\partial T}{\partial t}=\alpha\nabla^2T.}
$$

For our 2D grid, expand the Laplacian:

$$
\boxed{
\frac{\partial T}{\partial t}
=\alpha\left(\frac{\partial^2T}{\partial x^2}
+\frac{\partial^2T}{\partial y^2}\right).
}
$$

This is the **two-dimensional heat equation**. It is a PDE because the unknown function has derivatives with respect to both space and time.

## 4. Convert the continuous plate into a grid

A computer cannot store the temperature at every real-valued $(x,y)$. Divide the plate into discrete cells indexed by $i,j$:

$$T^n_{i,j}\approx T(x_i,y_j,t_n).$$

Here $n$ is the **time-step index**, and $i,j$ are the two **spatial indices**. On a uniform grid:

$$
 x_i=i\Delta x,\quad y_j=j\Delta x,\quad t_n=n\Delta t.
$$

In code:

```python
N = 50
T = np.full((N, N), 20.0)
T[20:30, 20:30] = 100.0
```

This creates a $50\times50$ temperature array (2,500 values). The slice `[20:30, 20:30]` selects a $10\times10$ hot square. All other cells start at 20°C.

## 5. Approximate the second spatial derivatives

Use a centered finite-difference approximation to the second derivative in the $x$ direction:

$$
\frac{\partial^2T}{\partial x^2}
\bigg|_{i,j}
\approx \frac{T_{i+1,j}-2T_{i,j}+T_{i-1,j}}{\Delta x^2}.
$$

Do the same for $y$, using the same spatial step:

$$
\frac{\partial^2T}{\partial y^2}
\bigg|_{i,j}
\approx \frac{T_{i,j+1}-2T_{i,j}+T_{i,j-1}}{\Delta x^2}.
$$

Adding them yields the **five-point Laplacian stencil**:

$$
\boxed{
\nabla^2T_{i,j}\approx
\frac{T_{i+1,j}+T_{i-1,j}+T_{i,j+1}+T_{i,j-1}-4T_{i,j}}{\Delta x^2}.
}
$$

Intuition: compare each cell's temperature to its four immediate neighbors. A hot cell surrounded by colder neighbors has a negative Laplacian, so its temperature should decrease; a cold cell surrounded by hot neighbors has a positive Laplacian, so its temperature should rise.

The second-order central difference is a reasonable approximation when the temperature varies smoothly over a few grid cells.

## 6. Approximate the time derivative

Use **forward Euler** in time:

$$
\frac{\partial T}{\partial t}
\approx\frac{T_{i,j}^{n+1}-T_{i,j}^n}{\Delta t}.
$$

Insert the spatial stencil into the PDE and solve for the new temperature:

$$
\frac{T_{i,j}^{n+1}-T_{i,j}^n}{\Delta t}
=\alpha\frac{T_{i+1,j}^n+T_{i-1,j}^n+T_{i,j+1}^n+T_{i,j-1}^n-4T_{i,j}^n}{\Delta x^2}.
$$

Rearrange:

$$
\boxed{
T_{i,j}^{n+1}=T_{i,j}^n+r\left(
T_{i+1,j}^n+T_{i-1,j}^n+T_{i,j+1}^n+T_{i,j-1}^n-4T_{i,j}^n
\right),
}
$$

with the dimensionless diffusion parameter:

$$\boxed{r=\frac{\alpha\Delta t}{\Delta x^2}.}$$

This is the entire **numerical model**: calculate the next grid from the current grid using four neighbors per interior cell. You do **not** need SciPy to execute this particular simulation; NumPy is sufficient for the update.

## 7. Understand numerical stability

An explicit method can become unstable if $\Delta t$ is too large. For the 2D five-point stencil with equal grid spacings, a standard stability condition is:

$$
\boxed{r\leq\frac14\quad\Longleftrightarrow\quad
\Delta t\leq\frac{\Delta x^2}{4\alpha}.}
$$

In this simulation:

$$
 r=\frac{1.0\times0.1}{1.0^2}=0.1\leq0.25.
$$

So the chosen time step satisfies the condition. Another way to see the update is:

$$
T_{i,j}^{n+1}=(1-4r)T_{i,j}^n
+r(T_{i+1,j}^n+T_{i-1,j}^n+T_{i,j+1}^n+T_{i,j-1}^n).
$$

When $0\leq r\leq1/4$, these weights are nonnegative and sum to 1. Therefore each new interior temperature is a **weighted average** of the current cell and its neighbors, consistent with smooth diffusion rather than explosive oscillations.

## 8. Implement the update efficiently with NumPy

The script updates all **interior cells simultaneously** using array slices:

```python
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
```

The index expressions mean:

| Slice | Neighbor |
| --- | --- |
| `T[1:-1, 1:-1]` | Current interior cell |
| `T[2:, 1:-1]` | Next row |
| `T[:-2, 1:-1]` | Previous row |
| `T[1:-1, 2:]` | Next column |
| `T[1:-1, :-2]` | Previous column |

**Why make a copy?** All values of $T^{n+1}$ must be computed from the **same old state** $T^n$. If we overwrote cells one by one and read partially updated neighbors, we would implement a different algorithm.

**What happens at the edges?** The code does not update the outermost row or column. Because they are copied unchanged, the boundary cells remain at **20°C**: an implicit **fixed-temperature (Dirichlet) boundary condition**. This means heat may effectively leave the modeled interior through boundaries maintained by an external thermal reservoir; do not assume total grid energy stays constant forever.

## 9. Visualize and export

`imshow` displays each cell as a colored square, with a hot color map and a color bar:

```python
im = ax.imshow(T, cmap="hot", vmin=20, vmax=100)
plt.colorbar(im, ax=ax, label="Temperature (°C)")

animation = FuncAnimation(
    fig, update, frames=300,
    interval=1000 / 30, blit=True
)
animation.save("videos/heat_diffusion.gif", writer="pillow", fps=30)
```

There are 300 steps of length $\Delta t=0.1$, so the numerical evolution covers approximately **30 simulation-time units**, while the GIF playback lasts **10 seconds** at 30 FPS. The axis labels in degrees Celsius describe the stored temperature values, not calibrated real-world evolution rates.

## 10. Experiments and extensions

- Increase $\alpha$: heat spreads more rapidly, but verify the stability condition.
- Decrease $\Delta x$: a finer grid can improve spatial resolution, but the maximum stable $\Delta t$ decreases **quadratically**.
- Set the boundary to a constant different from 20°C, or implement **insulated (zero-flux)** edges instead.
- Extend the equation to include a heat source $S(x,y,t)$, or use an implicit solver for larger time steps.
- Track total temperature above ambient over time; with fixed-temperature edges, it will eventually decrease when heat reaches those edges.

**Run from the repository root:**

```bash
python simulations/heat_diffusion.py
```
