# TODO.md

## Goal
Finish the poster with a single clear narrative: **small, production-grade code changes produce measurable ODE performance improvements in SciML** (Rössler as the case study), supported by a **reproducible benchmark + figure pipeline** and **clean, readable visuals**.

---
## 1) Lock scope and narrative (do first)
- [x] Intro to Rössler; then illustrate out-of-place vs in-place
updating in pass-by-reference language.
- Disscus stack allocated "static/fixed" arrays; display decrease
in memory allocations. (look at old commits -- this writeup was in
more detail at one time)
- Disscuss "in-line" and "in-bounds" functions. Explain these optimizations leverage the llvm compiler framework.
- [ ] Write/confirm the **performance ladder** variants you will show (6–8 max):
  1) naive out-of-place (allocates)
  2) in-place (alloc-free per RHS call)
  3) static small-system (SVector / stack allocation)
  4) optional: type-stable out-of-place vector
  5) optional: AD-ready tiny system
- [ ] Choose one **headline experiment spec** (single source of truth):
  - Recommended: fixed-step method (e.g., RK4 fixed dt or Midpoint fixed dt) with consistent tolerances/output settings.
- [ ] Define the **baseline** and **speedup formula**:
  - Baseline: `rossler_naive` (or whichever you decide—must be consistent everywhere)
  - Speedup: `median_time(baseline) / median_time(variant)`

**Acceptance criteria**
- Poster claims (e.g., “10×”) match the computed headline number for the chosen experiment.
---

## 2) Make the numerical experiment authoritative (single source of truth)
- [ ] Ensure the benchmark driver produces consistent results for:
  - RHS timing
  - Solve timing
  - allocations/bytes
- [ ] Ensure the experiment uses consistent:
  - `tspan`, `dt` (if fixed-step), tolerances (if adaptive), `saveat`/output disabling
  - warmup policy and `BenchmarkTools` setup

## 2.Strech) Add `C` and `Python` Apples-to-Apples RK4 fixed comparisons. (ideally, ``DifferentialEquations.jl`` wins and we emphasize it's bad to roll our own)
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
**Acceptance criteria**
- Running the benchmark script twice yields stable median ordering (minor variance acceptable).

## 3) Build a one-command figure pipeline (`./poster/figures/`)
Create/finish a script (recommended path):
- `poster/make_poster_figures.jl`

This script must:
- [ ] Run the canonical study entry point (e.g., `run_studies()`).
- [ ] Flatten results into a single DataFrame:
  - columns: `solver`, `variant`, `metric`, `median_time_ns`, `median_time_s`, `memory`, `allocs`
- [ ] Write all poster assets to: `./poster/figures/`
- [ ] Produce **exactly the figures referenced by `poster.tex`**.

### Minimum figure set (recommended)
- [ ] **Figure A — Solve time bar chart**
  - `midpoint_solve_times.png` (or `rk4_solve_times.png`)
  - Show 6–8 variants; consider log y-axis if dynamic range is large.
- [ ] **Figure B — RHS + allocations**
  - Either one combined figure or two separate figures:
    - median RHS time
    - allocations/bytes per RHS call (high impact “production lesson”)
- [ ] **Figure C — Solution sanity check**
  - A single 2D projection or 3D trajectory plot for one chosen variant (same ICs/params).
- [ ] Generate the LaTeX table from data (avoid drift):
  - `./poster/figures/midpoint_table.tex`

**Acceptance criteria**
- One command regenerates all assets:
  - `julia --project -e 'include("poster/make_poster_figures.jl")'`
- `poster.tex` compiles without manual figure/table edits.

---

## 4) Replace fragile TikZ diagram with a stable included graphic
- [ ] Remove the large TikZ pipeline diagram block from `poster.tex`.
- [ ] Replace with `\includegraphics{./figures/pipeline_overview.(pdf|png)}`.
- [ ] Generate `pipeline_overview` via:
  - Mermaid → SVG/PDF, or
  - a minimal standalone TikZ file compiled separately, or
  - a simple vector diagram tool of choice

**Acceptance criteria**
- Poster compiles reliably with fewer TikZ dependencies and faster iteration.

---

## 5) Make poster code listings match the repository (no copy/paste drift)
- [ ] Replace hand-copied listings with `\lstinputlisting` slices from source files:
  - `\lstinputlisting[firstline=..., lastline=...]{../rossler/impl.jl}`
- [ ] Show only what supports the ladder (keep snippets short).

**Acceptance criteria**
- Code shown on the poster is identical to repo code (or directly included from it).

---

## 6) Content polish: make claims consistent and readable at poster distance
- [ ] Replace any placeholders (e.g., “XX-fold speedup”) with computed values.
- [ ] Ensure captions specify:
  - what is measured (solve vs RHS),
  - baseline,
  - units,
  - whether output is disabled / saveat policy.
- [ ] Trim prose to poster-friendly density:
  - 3–5 bullets per section max
  - prioritize numeric results and what they imply

**Acceptance criteria**
- A reader standing ~1–2 meters away can follow the story without reading everything.

---

## 7) Write a “production checklist” conclusion (right column)
Replace generic conclusions with actionable takeaways:
- [ ] 3 key takeaways (allocations, type stability, data layout)
- [ ] Reproducibility block:
  - repo path
  - commands to generate figures
  - command to build poster
- [ ] Limits / next steps (1–2 bullets):
  - e.g., behavior on larger systems; ensemble runs; threading; stiff vs nonstiff; AD/Jacobian pipeline
- [ ] Possible plug how the compiler is nondeterministic. 

**Acceptance criteria**
- The conclusion reads like a field guide, not a summary.

---

## 8) Definition of Done (DoD)
You are done when all are true:
- [ ] `latexmk -pdf poster.tex` succeeds cleanly.
- [ ] All figures/tables in `./poster/figures/` are generated from a single Julia script run.
- [ ] The headline speedup number appears in **exactly two places** (Abstract + Introduction) and matches the plotted data.
- [ ] The poster references **no missing files**, and no “manual edits” are required after figure generation.
- [ ] Total visuals per column are kept reasonable (readability > completeness).

---

## Commands (recommended)
- Regenerate poster assets:
  - `julia --project -e 'include("poster/make_poster_figures.jl")'`
- Build poster:
  - `latexmk -pdf poster.tex`

---
