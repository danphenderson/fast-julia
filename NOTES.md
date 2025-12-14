# References:

1. [This is a common code pattern from high-level languages like MATLAB, SciPy, or R's deSolve. However, the issue with this form is that it allocates a vector, [dx,dy,dz], at each step. Let's benchmark the solution process with this choice of function"](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/faster_ode_example/#:~:text=R%27s%20deSolve,with%20this%20choice%20of%20function)


# Questions

1. Do techniques scale to "small-systems" instead
of the "tiny" Rössler system?

2. Do results change on a "stiff" problem (e.g,
the Lorenz System)? 
- Or, are they consistent with non-stiff and the only difference
is that we must use an adaptive step size, making it hard to compare
various implementations.

