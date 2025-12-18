[Rössler System Animation](./rossler.gif)

*Figure 1 — Performance ladder for Rössler RHS variants (L0–L5).*
```mermaid
flowchart TB
  %% Performance / “production-grade” ladder for Rössler RHS variants (increasing maturity/perf left→right)

  A["L0: Baseline (naïve, out-of-place)<br/>rossler_naive(u,p,t)<br/>• allocates new Vector each call<br/>• bounds checks<br/>• simplest reference"]
    --> B["L1: Bounds-check-free (out-of-place)<br/>rossler(u,p,t)<br/>• still allocates Vector<br/>• @inbounds / @inline<br/>• faster RHS, same API"]
    --> C["L2: Allocation-free (in-place)<br/>rossler!(du,u,p,t)<br/>• no per-call heap allocation<br/>• preallocated du<br/>• standard SciML performance step"]
    --> D["L3: Concrete eltype (type-stable OOP)<br/>rossler_type_stable(u,p,t)<br/>• allocates, but concrete eltype<br/>• robust under mixed types / AD<br/>• better compiler specialization"]
    --> E["L4: Stack-resident tiny system (StaticArrays)<br/>rossler_static(u::SVector,p,t)<br/>• allocation-free<br/>• SVector in/out<br/>• best for very small state (3D)"]
    --> F["L5: AD-ready + stack-resident<br/>rossler_ad(u::SVector,p::SVector,t)<br/>• allocation-free<br/>• explicit promote_type<br/>• forward-mode AD friendly"]

  %% Side notes
  B -.-> N1["Note: Out-of-place variants allocate by construction<br/>(useful for clarity / correctness baselines)"]
  C -.-> N2["Note: In-place is the typical “production” RHS<br/>(minimizes GC pressure)"]
  E -.-> N3["Note: StaticArrays often wins for tiny systems<br/>(stack/register temporaries; low overhead)"] --> 
  F["L5: AD-ready + stack-resident\nrossler_ad(u::SVector,p::SVector,t)\n• allocation-free\n• explicit promote_type\n• forward-mode AD friendly"]

  %% Optional side notes for the story
  B --- N1["Note: Out-of-place variants allocate by construction\n(useful for clarity / correctness baselines)"]
  C --- N2["Note: In-place is the typical 'production' RHS\n(minimizes GC pressure)"]
  E --- N3["Note: StaticArrays often wins for very small systems\n(less overhead, stack/register temporaries)"]
```

*Figure 2 — Mindmap of a practical small-ODE workflow (spec → method → validation → reporting).*
```mermaid
%%{init: {'theme': 'base', 'themeVariables': { 'primaryColor': '#3ec0c5', 'primaryBorderColor': '#0a1d1eff', 'primaryTextColor': '#0f2f30', 'lineColor': '#0f2f30', 'fontFamily': 'Helvetica Neue, Arial, sans-serif'}} }%%
mindmap
  root((Fast Julia: Small ODE Workflow))
    IVP setup
      "State and parameters: small n, scaled units"
      "Time window: t0 → tf, saveat grid"
      "Model checks: smooth f, invariants, events"
    Method choice
      "Stiffness guide: explicit RK4 / stabilized / BDF"
      "Accuracy target: fixed step vs adaptive tolerances"
      "Structure to exploit: symplectic flows, sparsity, Jacobian plan"
    Implementation
      "Allocate-free RHS: in-place `f!(du,u,p,t)`"
      "Type-stable inputs; promote eltypes"
      "Tiny systems: `StaticArrays` to stay on the stack"
    Validation
      "Sanity: bounds, invariants, positivity"
      "Convergence: refine dt or tolerances"
      "Sensitivity: sweep parameters, watch qualitative stability"
    Reporting
      "Plots: components, phase portrait, norms"
      "Performance: runtime plus allocations/GC"
      "Reproducibility: solver/tols, seeds, versions"
```

## Experiment Specification

*Figure 3 — Benchmark pipeline used in Experiment 1 (warmup → RHS microbench → solve sweep).*
```mermaid
flowchart TB
  %% Matches run_experiment1(): warmup compile, RHS microbench once/variant, solve bench per (variant, dt)

  S0["Experiment1Spec<br/>u0=[1,1,1], p=[0.1,0.1,14], tspan=(0,200)<br/>dt0=1e-2, halvings=6, warmup_steps=10"]
    --> S1["dt_schedule(dt0, halvings)<br/>dts = [dt0, dt0/2, dt0/4, ...]"]

  S1 --> S2["build_systems_for_experiment1(spec)<br/>Vector u0 for vector variants<br/>SVector u0 for static variants<br/>Optional rossler_ad if defined"]

  S2 --> W["Warm-up compilation<br/>tw=(t0, t0 + warmup_steps*dts[1])<br/>solve_fixed_rk4(prob_warm, dts[1])"]

  W --> R["RHS microbenchmark (dt-independent)<br/>bench_rhs(sys)<br/>• in-place: alloc du=similar(u)<br/>• out-of-place: call f(u,p,t)<br/>@benchmark ... samples=200 evals=1"]

  W --> P["Full solve benchmarks (dt sweep)<br/>prob = make_prob(sys) on full tspan<br/>for dt in dts: bench_solve(prob; dt)<br/>@benchmark ... samples=20 seconds=1 evals=1"]

  P --> K["RK4 fixed-step solve wrapper<br/>solve(prob, RK4(); dt=dt, ode_benchmark_solve_kwargs()...)<br/>minimal saving, dense=false, maxiters=1e12"]
```

*Figure 4 — RHS variant families and signatures used by the benchmark driver (Vector vs StaticArrays).*
```mermaid
flowchart LR
  %% Matches build_systems_for_experiment1() variant list and impl.jl signatures.

  subgraph V["Vector-based state (u::AbstractVector, p::AbstractVector)"]
    A["L0: OOP baseline (alloc + bounds)<br/>rossler_naive(u,p,t)<br/>returns Vector"] 
      --> B["L1: OOP + @inbounds/@inline (alloc)<br/>rossler(u,p,t)<br/>returns Vector"] 
    B --> C0["L2a: In-place baseline (no heap alloc; may have bounds)<br/>rossler_naive!(du,u,p,t)<br/>writes du"] 
      --> C1["L2b: In-place + @inbounds/@inline (no heap alloc)<br/>rossler!(du,u,p,t)<br/>writes du"]
    B --> D["L3: OOP type-stable (alloc, promoted eltype)<br/>rossler_type_stable(u,p,t)<br/>returns Vector{promote_type(Tu,Tp)}"]
  end

  subgraph S["Static state (u::SVector{3}, params vary by variant)"]
    E0["L4a: Static baseline (SVector output)<br/>rossler_static_naive(u,p,t)<br/>returns @SVector [...]"]
      --> E1["L4b: Static + @inbounds/@inline (allocation-free)<br/>rossler_static(u::SVector{3}, p, t)<br/>returns SVector{3,T}"]
    E1 --> F["L5: AD-ready static (allocation-free)<br/>rossler_ad(u::SVector{3}, p::SVector{3}, t)<br/>returns SVector{3,promote_type(Tu,Tp)}"]
  end

  %% Relationship between families
  B -.-> E1["Switch u0 to SVector in Experiment 1"]
  F -.-> N["Included only if rossler_ad is defined (conditional add)"]
```
# TODO.md

## Goal
Finish the poster with a single clear narrative: **small, production-grade code changes produce measurable ODE performance improvements in SciML** (Rössler as the case study), supported by a **reproducible benchmark + figure pipeline** and **clean, readable visuals**.

---

## 0) Resolve headline metric + claim consistency (blocker)
Your current poster text claims a “tenfold speedup,” but the current Midpoint fixed-step results you’ve been working with (as reflected in your table/fig placeholders) do **not** obviously support 10× solve-time speedup. Decide which metric is headline and make the claim match computed numbers everywhere.

## 1) Lock scope and narrative (do first)
- [x] Intro to Rössler; then illustrate out-of-place vs in-place updating in pass-by-reference language.
- [x] Add 2–3 poster-distance bullets explaining stack allocation (`StaticArrays`) and why allocations matter (GC pressure / throughput).
- [x] Add 2–3 poster-distance bullets explaining `@inbounds` / `@inline` as “remove bounds checks / encourage inlining,” and why that can matter in tight loops.
  - Keep it factual; avoid over-claiming compiler behavior.
- [x] Confirmed: **performance ladder variants are already implemented in code** (benchmark driver constructs 8 variants).
- [x] Choose which **6–8 variants** you will actually show on the poster (recommend 6 for clarity):
  - Suggested core ladder (6):
    1) `rossler_naive` (allocating)
    2) `rossler` (allocating + `@inbounds/@inline`)
    3) `rossler_naive!` (in-place “plain”)
    4) `rossler!` (in-place + `@inbounds/@inline`)
    5) `rossler_static_naive` (static “plain”)
    6) `rossler_static` (static + `@inbounds/@inline`)
  - Optional add-ons (only if they strengthen the story without clutter):
    - `rossler_type_stable`
    - `rossler_ad`
- [x] Define baseline and speedup formula (single choice, used everywhere):
  - Baseline: `rossler_naive`
  - Speedup: `median_time(baseline) / median_time(variant)`

**Acceptance criteria**
- Variant names in the poster match repository function names (no drift, no renaming).
- The narrative ladder is visible directly in the figures/table (same ordering).

---

## 2) Make the numerical experiment authoritative (single source of truth)
- [ ] Ensure the benchmark driver produces consistent results for:
  - RHS timing
  - Solve timing
  - allocations/bytes
- [ ] Ensure the experiment uses consistent settings and that `poster.tex` matches them exactly:
  - fixed-step: `tspan`, `dt`, solver (RK4), output disabling policy
  - (Avoid adaptive tolerances unless you explicitly add an adaptive experiment panel)
- [ ] Document benchmark policy in one place (script header or poster “Experiment Spec” box):
  - warmup policy
  - `BenchmarkTools` configuration
  - whether output is disabled (e.g., `save_on=false`, `dense=false`, etc.)
- [ ] Make the parameter story consistent in text:
  - “canonical parameters” used for the attractor visualization vs. the benchmark parameters (currently these differ in `poster.tex`)

**Acceptance criteria**
- Running the benchmark script twice yields stable median ordering (minor variance acceptable).
- Captions state: metric, baseline, units, solver spec, and output policy.

---

## 2.Stretch) Add `C` and `Python` apples-to-apples RK4 fixed comparisons
Only do this if it does not jeopardize the poster completion. Treat as optional appendix / “next step” unless already nearly done.

- [ ] Define apples-to-apples rules:
  - Same RK4 algorithm
  - Same dt, tspan, initial condition, parameters
  - Same output policy (prefer no saving / minimal saving)
- [ ] C variants (same algorithm, different compilation / flags):
  - C1: `-O3` (no fast-math)
  - C2: `-O3 -ffast-math`
  - C3: optional: `-Ofast -march=native` (document implied fast-math)
- [ ] C implementation choices:
  - `double y[3]` (stack) vs `struct` layout (optional)
  - `static inline` RHS + `restrict` pointers
- [ ] Julia variants to match poster ladder:
  - J1: naïve out-of-place Vector RHS (allocates each call)
  - J2: in-place Vector RHS (no per-call alloc)
  - J3: StaticArrays SVector RHS (stack/register-friendly)
  - J4: each above with `@fastmath` (separate benchmarks; label clearly)
- [ ] Use a Julia-written RK4 loop for apples-to-apples with C/Python.
  - (Comparing to `DifferentialEquations.jl` is valuable but is a *separate* experiment unless algorithm is identical.)

**Acceptance criteria**
- All implementations run the same spec and report comparable metrics.

---

## 3) Build a one-command figure pipeline (`./poster/figures/`)
Create/finish:
- `poster/make_poster_figures.jl`

This script must:
- [ ] Run the canonical study entry point (e.g., `run_studies()` or a single-case runner).
- [ ] Flatten results into one DataFrame:
  - columns: `solver`, `variant`, `metric`, `median_time_ns`, `median_time_s`, `memory`, `allocs`
- [ ] Write all poster assets to: `./poster/figures/`
- [ ] Produce **exactly the figures referenced by `poster.tex`** (no extra, no missing).

### Minimum figure set (updated to match your current poster intent)
- [ ] **Figure A — Solve-time speedup (log scale)**
  - Currently referenced as: `rk4_time_normalized_log10.png`
  - Keep speedup definition consistent with captions.
- [ ] **Figure B — Solve-time bar chart (normalized)**
  - Create a second, distinct asset (do not reuse the same file twice).
  - Update `poster.tex` to reference the new filename.
- [ ] **Figure C — RHS + allocations**
  - RHS median time + allocations/bytes per RHS call (combined or two panels).
- [ ] **Figure D — Attractor / solution sanity check**
  - Use this to satisfy the “Add Figure: Rössler attractor” placeholder in the Introduction.
- [ ] Generate LaTeX table from data (avoid drift):
  - `./poster/figures/rk4_table.tex` (or a naming scheme you standardize)

**Acceptance criteria**
- One command regenerates all assets:
  - `julia --project -e 'include("poster/make_poster_figures.jl")'`
- `poster.tex` compiles without manual figure/table edits.
- `poster.tex` uses `\input{./figures/<table>.tex}` rather than a hand-typed table.


## 6) Content polish: make claims consistent and readable at poster distance
- [ ] Replace placeholders (“XX-fold speedup”) with computed values.
- [ ] Captions must specify:
  - what is measured (solve vs RHS),
  - baseline,
  - units,
  - output policy (saving disabled / saveat behavior),
  - solver spec (Midpoint fixed dt, RK4 fixed dt, etc.)
- [ ] Trim prose density:
  - 3–5 bullets per section max
  - prioritize numeric results and what they imply

**Acceptance criteria**
- A reader at ~1–2 meters can follow the story without reading everything.

---

## 7) Write a “production checklist” conclusion (right column)
Your current conclusion veers into symbolic computation / ModelingToolkit. Unless you add results supporting that, keep it as “next steps” only.

- [ ] Replace conclusion with a field-guide checklist:
  - 3 key takeaways (allocations, type stability, data layout/macros)
- [ ] Add reproducibility block:
  - repo path
  - commands to generate figures
  - command to build poster
- [ ] Limits / next steps (1–2 bullets):
  - larger systems
  - ensemble runs
  - stiff vs nonstiff
  - AD/Jacobian pipeline
  - (Optional) ModelingToolkit as next step only if not demonstrated
- [ ] Optional: note compiler/benchmark nondeterminism carefully (variance, warmup, pinning).

**Acceptance criteria**
- Conclusion reads like actionable guidance, not a generic summary.

---

## 8) Definition of Done (DoD)
You are done when all are true:
- [ ] `latexmk -pdf poster.tex` succeeds cleanly.
- [ ] All figures/tables in `./poster/figures/` are generated from a single Julia script run.
- [ ] The headline speedup number appears in **exactly two places** (Abstract + Introduction) and matches plotted data.
- [ ] The poster references **no missing files**, and no “manual edits” are required after figure generation.
- [ ] Total visuals per column are kept reasonable (readability > completeness).

---

## Commands (recommended)
- Regenerate poster assets:
  - `julia --project -e 'include("poster/make_poster_figures.jl")'`
- Build poster:
  - `latexmk -pdf poster.tex`
