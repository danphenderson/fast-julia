![Rössler System Animation](./rossler.gif)




Benchmarks are compared via result object:

## Benchmark matrix (what variants you will compare)
C variants (same algorithm, different compilation / data layout)
  - C1: baseline: -O3 (no fast-math)
  - C2: fast-math: -O3 -ffast-math
  - C3: aggressive (optional): -Ofast -march=native (often implies fast-math; document it explicitly)

Implementation choices inside C:
  - State layout: double y[3] (stack) vs struct {double x,y,z;} (optional)
  - Use static inline RHS + restrict pointers to help optimization

Julia variants (match your poster narrative)
  - J1: naïve out-of-place Vector RHS (allocates each call)
  - J2: in-place Vector RHS (no per-call alloc)
  - J3: StaticArrays SVector RHS (stack/register-friendly)
  - J4: each of the above with @fastmath (separate benchmarks)
  - (Optional) J5: “production grade” version: @inbounds, @inline, - - preallocated temporaries, no globals

Important: Use a Julia-written RK4 loop for apples-to-apples with C or with python. (Comparing to DifferentialEquations.jl is valuable, but that’s a separate experiment because the algorithm differs unless you force fixed-step RK4.)


## Julia Implementations
  1. out-of-place
  2. in-place
  3. static
  4. annotated to support AD
  7. @inbounds and @inline


```json
{
  "lang":"C",
  "variant":"O3_ffast-math",
  "scenario":"final_state",
  "steps":20000,
  "time_ns_median":123456789,
  "ns_per_step":6172.8,
  "notes":"clang-17, -march=native"
}
```

```mermaid
mindmap
  root(Solving Small System of ODEs)
    IVP Specification
      \newlines["`State & parameters
        dimension n
        parameters p
        units/scales -->
      `"]
      \newlines["`Time domain
        t0, tf
        output times (saveat)`"]
      \newlines["`IC x0 (and x0' if DAE)`"]
      \newlines["`Model properties
        smoothness of f
        invariants/constraints (mass, positivity)
        events/discontinuities`"]
    Well-posedness checks
      \newlines["`Existence & uniqueness (Picard–Lindelöf)
        f locally Lipschitz in x
      `"]
      \newlines["`Domain validity
        singularities / blow-up
        parameter ranges
      `"]
      \newlines["`Scaling / nondimensionalization
        reduce condition issues
        balance magnitudes
      `"]
    Numerical method selection
      \newlines["`Stiffness?
        nonstiff: explicit RK (e.g., RK4 / Tsit5)
        mildly stiff: stabilized explicit / IMEX
        stiff: implicit (BDF, Rosenbrock, Radau)
      `"]
      \newlines["`Accuracy target
        low/medium: fixed step
        high: adaptive step + error control
      `"]
      \newlines["`Structure to exploit
        Hamiltonian / symplectic: symplectic integrator
        oscillatory: methods tuned for oscillations
        large sparse: Jacobian sparsity + Krylov
      `"]
          \newlines["`Discretization details
      Step size / tolerance
        dt (fixed) OR reltol/abstol (adaptive)
        maxiters, dtmin/dtmax
      Error estimation
        embedded pair
        local truncation error
      Jacobian handling (implicit)
        analytic vs AD vs finite diff
        factorization / preconditioner
    `"]   
    \newlines["`Implementation pipeline
      Define f(t,x,p)
        allocate-free / in-place
        type stability
      Choose solver + options
        saveat / dense output
        callbacks (events)
      Run solve
        check return code
    `"]   
    \newlines["`Validation & verification
      Sanity checks
        conservation / invariants
        positivity / bounds
      Convergence
        refine dt or tighten tolerances
        compare against reference solution
      Sensitivity
        vary parameters p
        check qualitative stability
    `"]   
    \newlines["`Diagnostics & reporting
      Plots
        components vs t
        phase portraits
        norms / energy
      Performance
        allocations / profiling
        stiffness indicators
      Reproducibility
        fixed seeds (if noise)
        record solver/tols/version
    `"]
```















```mermaid
flowchart TD
    markdown["`Solve ODE(S) IVP`"]
    newLines["`Initial State
    ODE System
    Time Span
    Error Tolerance`"]
    markdown --> newLines
```

```mermaid
flowchart LR
  %% ------------------------------------------------------------
  %% ODE Pipeline (Model → Implement → Solve → Analyze → Deliver)
  %% ------------------------------------------------------------

  
```
