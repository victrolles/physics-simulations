# Double Pendulum: Positions, Velocities, and Lagrange's Equations

[← Back to README](../README.md) · [View Python source](../simulations/double_pendulum.py)

A double pendulum has **two degrees of freedom**: the angles of its two rods. Unlike a simple pendulum, the motion of each mass depends on the other. In this tutorial, we derive the complete equations implemented in `simulations/double_pendulum.py`, **starting from geometry**, then differentiating coordinates to get velocities, forming energies, and applying **Euler–Lagrange equations**.

## 1. Define the system and generalized coordinates

We assume two point masses $m_1,m_2$, two massless rigid rods of lengths $L_1,L_2$, frictionless joints, and uniform gravity $g$. The system moves in a vertical plane.

Both angles are measured from the **downward vertical**:

- $\theta_1$ is the absolute angle of the first rod.
- $\theta_2$ is the absolute angle of the second rod, **not** the angle relative to the first.
- $\omega_1=\dot\theta_1$ and $\omega_2=\dot\theta_2$ are angular velocities.
- $\Delta=\theta_1-\theta_2$ is the difference between angles.

The code uses $L_1=L_2=1\ \mathrm m$, $m_1=m_2=1\ \mathrm{kg}$ and $g=9.81\ \mathrm{m/s^2}$. Initial conditions are $\theta_1(0)=120^\circ$, $\theta_2(0)=-10^\circ$, and $\omega_1(0)=\omega_2(0)=0$.

We choose *generalized coordinates* $q_1=\theta_1$ and $q_2=\theta_2$. They automatically respect constant rod lengths, so we do not need to simulate rod tension explicitly.

## 2. Compute the coordinates of M1 and M2

Place the fixed top pivot at $(0,0)$, and take positive $y$ to point **upward**.

The first mass $M_1$ is exactly $L_1$ away from the pivot:

$$
\boxed{
\begin{aligned}
x_1 &= L_1\sin\theta_1,\\
y_1 &=-L_1\cos\theta_1.
\end{aligned}
}
$$

The second mass $M_2$ is located a distance $L_2$ away from $M_1$ along the second rod. **Add** that displacement to the position of $M_1$:

$$
\boxed{
\begin{aligned}
x_2 &= L_1\sin\theta_1+L_2\sin\theta_2,\\
y_2 &=-L_1\cos\theta_1-L_2\cos\theta_2.
\end{aligned}
}
$$

These are the formulas used for drawing:

```python
x1 = L1 * np.sin(theta1)
y1 = -L1 * np.cos(theta1)
x2 = x1 + L2 * np.sin(theta2)
y2 = y1 - L2 * np.cos(theta2)
```

Check: $(x_1)^2+(y_1)^2=L_1^2$ and $(x_2-x_1)^2+(y_2-y_1)^2=L_2^2$. Therefore both rods keep their lengths.

## 3. Differentiate to obtain V1 and V2

Velocity is the time derivative of position. The chain rule gives:

$$
\frac{d}{dt}\sin\theta=\cos\theta\,\dot\theta,
\qquad
\frac{d}{dt}\cos\theta=-\sin\theta\,\dot\theta.
$$

For $M_1$:

$$
\boxed{
\begin{aligned}
v_{1x}&=\dot x_1=L_1\cos\theta_1\,\omega_1,\\
v_{1y}&=\dot y_1=L_1\sin\theta_1\,\omega_1.
\end{aligned}
}
$$

Its **speed squared**, obtained by squaring and adding the components, simplifies using $\sin^2\theta+\cos^2\theta=1$:

$$
\begin{aligned}
v_1^2&=v_{1x}^2+v_{1y}^2\\
&=L_1^2\omega_1^2(\cos^2\theta_1+\sin^2\theta_1)\\
&=\boxed{L_1^2\omega_1^2}.
\end{aligned}
$$

For $M_2$, **both rods contribute** to the velocity:

$$
\boxed{
\begin{aligned}
v_{2x}&=L_1\cos\theta_1\,\omega_1+L_2\cos\theta_2\,\omega_2,\\
v_{2y}&=L_1\sin\theta_1\,\omega_1+L_2\sin\theta_2\,\omega_2.
\end{aligned}
}
$$

Expand $v_{2x}^2+v_{2y}^2$. The mixed terms combine by the identity

$$
\cos\theta_1\cos\theta_2+\sin\theta_1\sin\theta_2
=\cos(\theta_1-\theta_2)=\cos\Delta.
$$

Hence:

$$
\boxed{
v_2^2=L_1^2\omega_1^2+L_2^2\omega_2^2
+2L_1L_2\omega_1\omega_2\cos\Delta.
}
$$

The last term encodes the **mechanical coupling** between the two degrees of freedom.

## 4. Calculate kinetic energy K

The two masses have translational kinetic energies:

$$K=\frac12m_1v_1^2+\frac12m_2v_2^2.$$

Substitute the squared speeds and group like terms:

$$
\boxed{
K=\frac12(m_1+m_2)L_1^2\omega_1^2
+\frac12m_2L_2^2\omega_2^2
+m_2L_1L_2\omega_1\omega_2\cos\Delta.
}
$$

Notice that $m_2$ contributes to the first rod's kinetic term: when the first rod rotates, it moves *both* masses.

## 5. Calculate gravitational potential energy V

With the origin at the top pivot, gravitational potential energy is $V=m_1gy_1+m_2gy_2$:

$$
\begin{aligned}
V &= -m_1gL_1\cos\theta_1 \\
&\quad -m_2g(L_1\cos\theta_1+L_2\cos\theta_2).
\end{aligned}
$$

Combining terms gives:

$$
\boxed{V=-(m_1+m_2)gL_1\cos\theta_1-m_2gL_2\cos\theta_2.}
$$

We may add any constant to $V$ without changing the motion. This expression simply chooses the zero of gravitational potential at the height of the pivot.

## 6. Construct the Lagrangian

Define the **Lagrangian** $\mathcal L=K-V$:

$$
\begin{aligned}
\mathcal L={}&\frac12(m_1+m_2)L_1^2\omega_1^2
+\frac12m_2L_2^2\omega_2^2\\
&+m_2L_1L_2\omega_1\omega_2\cos\Delta\\
&+(m_1+m_2)gL_1\cos\theta_1
+m_2gL_2\cos\theta_2.
\end{aligned}
$$

For each generalized coordinate $q_i\in\{\theta_1,\theta_2\}$, apply **Euler–Lagrange**:

$$
\boxed{
\frac{d}{dt}\left(\frac{\partial\mathcal L}{\partial\dot q_i}\right)
-\frac{\partial\mathcal L}{\partial q_i}=0.
}
$$

The $\partial$ symbols are *partial derivatives* of the energy expression with respect to one argument at a time. **They do not make the pendulum equations PDEs**: the unknown angles depend only on time, so the final motion equations are **coupled nonlinear ODEs**. The 2D heat simulation, in contrast, does solve a PDE numerically.

## 7. Derive the first Euler–Lagrange equation

Take derivatives with respect to $\omega_1$ and $\theta_1$:

$$
\frac{\partial\mathcal L}{\partial\omega_1}
=(m_1+m_2)L_1^2\omega_1
+m_2L_1L_2\omega_2\cos\Delta,
$$

$$
\frac{\partial\mathcal L}{\partial\theta_1}
=-m_2L_1L_2\omega_1\omega_2\sin\Delta
-(m_1+m_2)gL_1\sin\theta_1.
$$

For the total time derivative, remember $\dot\Delta=\omega_1-\omega_2$ and use the product rule:

$$
\begin{aligned}
\frac{d}{dt}\left(\frac{\partial\mathcal L}{\partial\omega_1}\right)
={}&(m_1+m_2)L_1^2\ddot\theta_1\\
&+m_2L_1L_2\big[\ddot\theta_2\cos\Delta
-\omega_2(\omega_1-\omega_2)\sin\Delta\big].
\end{aligned}
$$

Substitute into Euler–Lagrange. Terms containing $\omega_1\omega_2\sin\Delta$ cancel. The result is:

$$
\boxed{
\begin{aligned}
&(m_1+m_2)L_1^2\ddot\theta_1
+m_2L_1L_2\cos\Delta\,\ddot\theta_2\\
&\qquad+m_2L_1L_2\sin\Delta\,\omega_2^2
+(m_1+m_2)gL_1\sin\theta_1=0.
\end{aligned}
}
$$

## 8. Derive the second Euler–Lagrange equation

Differentiate with respect to $\omega_2$ and $\theta_2$:

$$
\frac{\partial\mathcal L}{\partial\omega_2}
=m_2L_2^2\omega_2+m_2L_1L_2\omega_1\cos\Delta,
$$

$$
\frac{\partial\mathcal L}{\partial\theta_2}
=+m_2L_1L_2\omega_1\omega_2\sin\Delta
-m_2gL_2\sin\theta_2.
$$

The signs differ because $\partial\Delta/\partial\theta_2=-1$:

$$
\begin{aligned}
\frac{d}{dt}\left(\frac{\partial\mathcal L}{\partial\omega_2}\right)
={}&m_2L_2^2\ddot\theta_2\\
&+m_2L_1L_2\big[\ddot\theta_1\cos\Delta
-\omega_1(\omega_1-\omega_2)\sin\Delta\big].
\end{aligned}
$$

Substitution and cancellation of mixed terms gives:

$$
\boxed{
\begin{aligned}
&m_2L_1L_2\cos\Delta\,\ddot\theta_1
+m_2L_2^2\ddot\theta_2\\
&\qquad-m_2L_1L_2\sin\Delta\,\omega_1^2
+m_2gL_2\sin\theta_2=0.
\end{aligned}
}
$$

## 9. Rearrange into the matrix system used in Python

We have **two equations for two unknown angular accelerations**. Divide the first equation by $L_1$ and the second by $m_2L_2$ (assuming positive lengths and masses), then put known gravity/velocity terms on the right.

$$
\underbrace{\begin{bmatrix}
(m_1+m_2)L_1 & m_2L_2\cos\Delta\\
L_1\cos\Delta & L_2
\end{bmatrix}}_{A}
\underbrace{\begin{bmatrix}
\ddot\theta_1\\\ddot\theta_2
\end{bmatrix}}_{\boldsymbol\alpha}
=
\underbrace{\begin{bmatrix}
-m_2L_2\sin\Delta\,\omega_2^2-(m_1+m_2)g\sin\theta_1\\
L_1\sin\Delta\,\omega_1^2-g\sin\theta_2
\end{bmatrix}}_{B}.
$$

In compact form, $A\boldsymbol\alpha=B$. **Both $A$ and $B$ change every time the pendulum moves**, because they depend on current angles and angular velocities.

The implementation matches this derivation:

```python
delta = theta1 - theta2

A = np.array([
    [(M1 + M2) * L1, M2 * L2 * np.cos(delta)],
    [L1 * np.cos(delta), L2]
])

B = np.array([
    -M2 * L2 * np.sin(delta) * omega2**2
    - (M1 + M2) * g * np.sin(theta1),
    L1 * np.sin(delta) * omega1**2
    - g * np.sin(theta2)
])

alpha1, alpha2 = np.linalg.solve(A, B)
```

`np.linalg.solve(A, B)` computes the accelerations without explicitly forming a matrix inverse. The displayed $A$ is a **row-scaled** version of the underlying kinetic-energy mass matrix, which is why it is not necessarily symmetric.

## 10. Convert to four first-order ODEs

An IVP solver needs derivatives of the state. Define:

$$
\mathbf y=\begin{bmatrix}\theta_1\\\omega_1\\\theta_2\\\omega_2\end{bmatrix}.
$$

At each solver evaluation, calculate $\alpha_1,\alpha_2$ as above and return:

$$
\boxed{\dot{\mathbf y}=\begin{bmatrix}\omega_1\\\alpha_1\\\omega_2\\\alpha_2\end{bmatrix}.}
$$

```python
def double_pendulum(t, y):
    theta1, omega1, theta2, omega2 = y
    # Build A and B using the current angles and angular velocities.
    alpha1, alpha2 = np.linalg.solve(A, B)
    return [omega1, alpha1, omega2, alpha2]
```

The actual Python function defines `A` and `B` before calling `np.linalg.solve`.

The script uses `solve_ivp(..., rtol=1e-9, atol=1e-9)` for tight error tolerances. It requests 500 samples over **0–10 s**, although the solver's interval is **0–15 s**. Only the requested 0–10 s data are used for visualization.

## 11. Draw both rods and masses

Take angles from `solution.y[0]` and `solution.y[2]`, calculate both Cartesian trajectories with the formulas in Section 2, and update two lines and two markers:

```python
line1.set_data([0, x1[frame]], [0, y1[frame]])
line2.set_data([x1[frame], x2[frame]], [y1[frame], y2[frame]])
bob1.set_data([x1[frame]], [y1[frame]])
bob2.set_data([x2[frame]], [y2[frame]])
```

`FuncAnimation` exports 500 frames at 30 FPS to `videos/double_pendulum.gif`. Therefore GIF playback lasts ~16.7 s even though the sampled physical motion spans 10 s.

## 12. Validate and experiment

- **Constraint check:** The distances $|OM_1|$ and $|M_1M_2|$ should always equal $L_1,L_2$.
- **Energy check:** With no damping, $E=K+V$ should stay approximately constant (modulo numerical error).
- **Sensitive dependence:** Try changing $\theta_2(0)$ by $0.1^\circ$ and compare the paths. Different trajectories after some time are characteristic of chaotic motion, not randomness.
- **Mass ratio:** Vary $m_2/m_1$ to see how the coupling changes.

**Run from the repository root:**

```bash
python simulations/double_pendulum.py
```
