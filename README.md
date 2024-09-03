# Crank_Nicolson


The Crank-Nicolson method is a finite difference method used for numerically solving partial differential equations (PDEs), especially parabolic equations like the heat equation. It is known for being unconditionally stable and second-order accurate in both time and space.


Consider the one-dimensional heat equation:
∂u/∂t = α * ∂²u/∂x²
where:
u(x, t) is the temperature distribution as a function of space x and time t.
α (alpha) is the thermal diffusivity constant.

The spatial domain [0, L] is discretized into M points with a spatial step size Δx = L / (M - 1). The time domain [0, T] is discretized into N time steps with a time step size Δt = T / N.

Let u_i^n represent the approximation of u(x_i, t_n) at spatial position x_i and time t_n. The Crank-Nicolson method is a time-stepping method that is implicit, meaning it requires solving a system of equations at each time step.


The Crank-Nicolson method takes the average of the explicit (forward in time) and implicit (backward in time) schemes. The finite difference approximation is given by:
(u_i^(n+1) - u_i^n) / Δt = (α / 2) * [ (u_(i+1)^(n+1) - 2u_i^(n+1) + u_(i-1)^(n+1)) / Δx² + (u_(i+1)^n - 2u_i^n + u_(i-1)^n) / Δx² ]

Rearranging terms:
(-r/2) * u_(i-1)^(n+1) + (1 + r) * u_i^(n+1) - (r/2) * u_(i+1)^(n+1) = (r/2) * u_(i-1)^n + (1 - r) * u_i^n + (r/2) * u_(i+1)^n
where r = (α * Δt) / Δx² is a dimensionless parameter.

This can be written in matrix form as:
A * u^(n+1) = B * u^n
where:
A and B are tridiagonal matrices.
u^(n+1) is the vector of unknowns at time step n+1.
u^n is the known vector from the previous time step.

The tridiagonal matrices A and B have the following form:
For matrix A:
A = [ 1+r   -r/2    0    ...    0  ]
    [-r/2   1+r   -r/2   ...    0  ]
    [  0   -r/2   1+r   ...    0  ]
    [ ...                 ...     ]
    [  0    ...   -r/2   1+r  -r/2]
    [  0    ...    0    -r/2  1+r ]
    
For matrix B:
B = [ 1-r   r/2     0    ...    0  ]
    [ r/2   1-r    r/2   ...    0  ]
    [  0    r/2    1-r   ...    0  ]
    [ ...                 ...     ]
    [  0    ...    r/2   1-r   r/2]
    [  0    ...     0    r/2   1-r ]

    
Time-Stepping Algorithm
Initialization: Set the initial condition u(x, 0) = u0(x).

Time Loop:
For each time step n:
Compute the right-hand side vector B * u^n.
Solve the linear system A * u^(n+1) = B * u^n for u^(n+1).
Boundary Conditions: Apply boundary conditions (e.g., Dirichlet or Neumann) at each time step.
