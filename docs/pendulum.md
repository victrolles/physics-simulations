# Simple Pendulum: Derive, Solve, and Animate the Motion

[← Back to README](../README.md) · [View Python source](../simulations/pendulum.py)

This tutorial reproduces the physics and numerical algorithm in `simulations/pendulum.py`, starting from **gravity and an initial angle**. We will derive the angular acceleration, turn the equation into a form that SciPy can integrate, and animate the resulting trajectory.

## 1. Define the physical model

We consider a point mass (the *bob*) attached to an inextensible, massless rod of length $L$, rotating around a fixed, frictionless pivot. Gravity is uniform. There is **no air resistance or damping**.

Define:

| Symbol | Meaning | Script value |
| --- | --- | --- |
| $m$ | Mass of the bob | Cancels from the final equation |
| $g$ | Gravitational acceleration | $9.81\ \mathrm{m/s^2}$ |
| $L$ | Distance from pivot to bob | $1.0\ \mathrm m$ |
| $\theta(t)$ | Angle from the downward vertical | Unknown function of time |
| $\omega(t)=\dot\theta(t)$ | Angular velocity | Unknown function of time |
| $\ddot\theta(t)$ | Angular acceleration | Unknown function of time |

A dot means a **time derivative**: $\dot\theta=d\theta/dt$, and $\ddot\theta=d^2\theta/dt^2$. Angles passed to NumPy's trigonometric functions are in **radians**, not degrees.

The initial conditions are:

$$
\theta(0)=\theta_0=\frac{\pi}{4}\quad(45^\circ),\qquad
\omega(0)=\omega_0=0.
$$

The bob starts from rest, but that does **not** mean its acceleration is zero.

## 2. Decompose the force of gravity

Gravity acts downward:

$$\vec F_g=m\vec g.$$

The bob is constrained to move on a circle, so we decompose gravity into a **radial** component (along the rod) and a **tangential** component (along the direction of motion). The rod's tension handles the radial constraint; the tangential component changes the angle.

Project gravity onto the tangent of the circle:

$$
\boxed{F_{\mathrm{tan}}=-mg\sin\theta.}
$$

Why $\sin\theta$? Because the tangent is perpendicular to the rod, and the projection of the vertical gravitational force onto it is proportional to the sine of the deflection angle.

Why the **minus sign**? If $\theta>0$ (bob displaced right), gravity points toward decreasing $\theta$; if $\theta<0$, gravity points toward increasing $\theta$. It is a *restoring* force. At $\theta=0$, $\sin\theta=0$ and tangential gravity vanishes.

## 3. Apply Newton's second law

The distance along the circular arc is $s=L\theta$ (with $\theta$ in radians). For a fixed $L$:

$$
\frac{ds}{dt}=L\dot\theta=L\omega
\qquad\text{and}\qquad
\frac{d^2s}{dt^2}=L\ddot\theta.
$$

Newton's second law along the tangent gives:

$$
F_{\mathrm{tan}}=ma_{\mathrm{tan}}
\quad\Rightarrow\quad
-mg\sin\theta=mL\ddot\theta.
$$

Cancel $m$ and divide by $L$:

$$
\boxed{\ddot\theta=-\frac{g}{L}\sin\theta.}
$$

This is the **nonlinear second-order ordinary differential equation (ODE)** solved in the script. The bob's mass cancels: for this ideal model, it does not change the motion.

The same result follows from torque. Gravitational torque about the pivot is $\tau=-mgL\sin\theta$; moment of inertia is $I=mL^2$; therefore $I\ddot\theta=\tau$ produces the same equation.

At the initial angle:

$$
\ddot\theta(0)=-\frac{9.81}{1}\sin(\pi/4)
\approx -6.94\ \mathrm{rad/s^2}.
$$

The initial angular **velocity** is zero; the angular **acceleration** is not.

## 4. Distinguish angle, velocity, and acceleration

The angular velocity is the first time derivative of the angle:

$$\omega(t)=\frac{d\theta}{dt}.$$

The angular acceleration is the first derivative of velocity:

$$
\frac{d\omega}{dt}=\ddot\theta=-\frac{g}{L}\sin\theta.
$$

We cannot just calculate a single velocity from the starting angle. Acceleration changes **at every instant** as $\theta(t)$ changes. We therefore integrate the evolving acceleration to obtain $\omega(t)$, and integrate velocity to obtain $\theta(t)$.

## 5. Rewrite as a first-order system

SciPy's `solve_ivp` operates on equations of the form $\dot{\mathbf y}=f(t,\mathbf y)$. Introduce a two-element **state vector**:

$$
\mathbf y=\begin{bmatrix}\theta\\\omega\end{bmatrix}.
$$

Its derivative is:

$$
\boxed{
\frac{d\mathbf y}{dt}
=\begin{bmatrix}
\omega\\
-(g/L)\sin\theta
\end{bmatrix}.
}
$$

The code translates the equations literally:

```python
def pendulum(t, y):
    theta, omega = y

    dtheta = omega
    domega = -(g / L) * np.sin(theta)

    return [dtheta, domega]
```

- `dtheta` represents $d\theta/dt=\omega$, **not** a small change in the angle during one step.
- `domega` represents $d\omega/dt=-(g/L)\sin\theta$.
- The function receives `t` because SciPy requires it, even though this model has no explicit time-dependent forcing.

## 6. Solve numerically with SciPy

```python
t = np.linspace(0, 10, 500)

solution = solve_ivp(
    pendulum,                   # Derivative of the state
    [0, 10],                    # Physical time interval, in seconds
    [theta0, omega0],          # Initial angle and angular velocity
    t_eval=t                   # Requested output times
)

theta = solution.y[0]
omega = solution.y[1]
```

The solver receives the current state $[\theta,\omega]$, computes its derivative, advances the system, and repeats. By default it uses **RK45**, an adaptive Runge–Kutta integrator, which can take smaller internal time steps when more precision is needed. The 500 values in `t_eval` specify where to **record** results, not necessarily the solver's internal step size.

For reliability you may also check `solution.success` and inspect `solution.message`. The numerical solution approximates the continuous equations; it is not a symbolic closed-form expression.

## 7. Convert the simulated angle to Cartesian coordinates

A plotting library needs $(x,y)$ rather than the angular state. Set the pivot at $(0,0)$, with $y$ increasing upward:

$$
\boxed{x(t)=L\sin\theta(t),\qquad y(t)=-L\cos\theta(t).}
$$

Why the negative sign in $y$? When $\theta=0$, the bob hangs directly **below** the pivot at $(0,-L)$.

```python
x = L * np.sin(theta)
y = -L * np.cos(theta)
```

## 8. Animate the solution with Matplotlib

Create a line (the rod) and a circular marker (the **bob**):

```python
line, = ax.plot([], [], "k-", linewidth=2)
bob, = ax.plot([], [], "bo", markersize=15)
```

On each frame, `update()` places the rod between the pivot and bob and moves the bob to the calculated position:

```python
def update(frame):
    line.set_data([0, x[frame]], [0, y[frame]])
    bob.set_data([x[frame]], [y[frame]])
    return line, bob
```

Finally:

```python
animation = FuncAnimation(
    fig, update, frames=len(theta),
    interval=1000 / 30, blit=True
)
animation.save("videos/pendulum.gif", writer="pillow", fps=30)
```

There are 500 animation frames displayed at 30 frames/s: approximately **16.7 seconds of GIF playback**, although the simulated physical times cover **10 seconds**. These two time scales are not synchronized in the current script; that affects playback speed but not the computed trajectory.

## 9. Validate and extend the model

Without friction, total mechanical energy is approximately conserved:

$$
E=\frac12 mL^2\omega^2+mgL(1-\cos\theta).
$$

Check whether $E(t)$ remains nearly constant, allowing for numerical error. For very small $|\theta|$, $\sin\theta\approx\theta$ and the equation becomes $\ddot\theta+(g/L)\theta=0$, with approximate period $T\approx2\pi\sqrt{L/g}$. **The script retains the full sine** instead of making this approximation.

Try varying $L$, $\theta_0$, or $\omega_0$, and observe how the motion changes.

**Run from the repository root:**

```bash
python simulations/pendulum.py
```
