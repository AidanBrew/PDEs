# Simulating Wound Healing

Consider how wound healing can be described by the Porous-Fisher equation:

$$
\frac{\partial u}{\partial t}
=
\frac{\partial}{\partial x}
\left(
\alpha
\left(\frac{u}{K}\right)^m
\frac{\partial u}{\partial x}
\right)
+
\beta u
\left(1-\frac{u}{K}\right).
$$

where $u(x,t)$ denotes the density of cells and $K > 0$ is the carrying capacity. Consider $0 \leq x \leq L$ and $t > 0$, subject to the following initial and boundary conditions:

$$
u(x,0)
=
f(x)
=
\begin{cases}
1, & x \leq \frac{L}{4} \text{ or } x \geq \frac{3L}{4}, \\[6pt]
0, & \frac{L}{4} < x < \frac{3L}{4}.
\end{cases}
$$

$$
\frac{\partial u}{\partial x}(0,t)=0,
\qquad
\frac{\partial u}{\partial x}(L,t)=0,
\qquad
t>0.
$$

The initial condition represents a one-dimensional cross-section of a wound located on the interval

$$
\frac{L}{4} < x < \frac{3L}{4},
$$

where the cell population is initially zero.

## 1. Finite-volume spatial discretisation

We discretise the initial-boundary value problem in space using the finite volume method. Averaging is used to discretise the nonlinear diffusion term, and a non-uniform mesh is allowed with nodes

$$
x=x_i,
\qquad
i=1,\dots,N,
$$

where

$$
x_1=0,
\qquad
x_N=L,
\qquad
N \in \mathbb{Z}^{+}.
$$

Let $u_i$ be the numerical approximation to $u(x_i,t)$. The spatial discretisation is expressed in the form

$$
\frac{d\vec{u}}{dt}
=
\vec{G}(\vec{u}),
\qquad
\vec{u}(0)
=
\vec{u}^{(0)}.
$$

Here,

$$
\vec{u}
=
(u_1,\dots,u_N)^T,
\qquad
\vec{u}^{(0)}
=
(f(x_1),\dots,f(x_N))^T,
$$

and

$$
\vec{G}(\vec{u})
=
(g_1(\vec{u}),\dots,g_N(\vec{u}))^T.
$$

## Discretisation final result

Define the nonlinear diffusivity and reaction term as

$$
D(u_i)
=
\alpha
\left(
\frac{u_i}{K}
\right)^m,
$$

$$
R(u_i)
=
\beta u_i
\left(
1-\frac{u_i}{K}
\right).
$$

For a non-uniform mesh, define the node spacings

$$
h_i
=
x_{i+1}-x_i,
\qquad
i=1,\dots,N-1.
$$

The control volume lengths are

$$
V_1
=
\frac{h_1}{2},
\qquad
V_i
=
\frac{h_{i-1}+h_i}{2},
\qquad
V_N
=
\frac{h_{N-1}}{2}.
$$

Using East and West notation, define the averaged diffusivities at the west and east faces as

$$
D_{W,i}
=
\frac{D(u_{i-1})+D(u_i)}{2},
\qquad
D_{E,i}
=
\frac{D(u_i)+D(u_{i+1})}{2}.
$$

The spatial derivatives at the west and east faces are approximated by

$$
\frac{\partial u}{\partial x}(w_i,t)
\approx
\frac{u_i-u_{i-1}}{h_{i-1}},
\qquad
\frac{\partial u}{\partial x}(e_i,t)
\approx
\frac{u_{i+1}-u_i}{h_i}.
$$

Then,

$$
g_i(\vec{u})
=
\begin{cases}
\dfrac{1}{V_1}
\left[
- D_{E,1}
\dfrac{\partial u}{\partial x}(e_1,t)
\right]
+
R(u_1),
& i=1,
\\[12pt]
\dfrac{1}{V_i}
\left[
D_{W,i}
\dfrac{\partial u}{\partial x}(w_i,t)
-
D_{E,i}
\dfrac{\partial u}{\partial x}(e_i,t)
\right]
+
R(u_i),
& i=2,\dots,N-1,
\\[12pt]
\dfrac{1}{V_N}
\left[
D_{W,N}
\dfrac{\partial u}{\partial x}(w_N,t)
\right]
+
R(u_N),
& i=N.
\end{cases}
$$
