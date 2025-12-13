module PosterFigures

using Statistics
using DataFrames
using Plots

# Poster-friendly defaults (override in your own code if desired)
default(
    dpi = 300,
    titlefontsize = 16,
    guidefontsize = 14,
    tickfontsize = 12,
    legendfontsize = 10,
)

# -----------------------------------------------------------------------------
# Small utilities
# -----------------------------------------------------------------------------

"""Coerce a value to a Symbol, accepting Symbols or Strings."""
as_symbol(x) = x isa Symbol ? x : Symbol(x)

"""A conservative, file-name-safe slug for solver titles."""
function solver_slug(s::AbstractString)
    t = lowercase(String(s))
    t = replace(t, r"[^a-z0-9]+" => "_")
    t = replace(t, r"^_+|_+$" => "")
    return t
end

"""Make variant labels a bit more human-friendly for poster figures."""
function pretty_variant(v::AbstractString)
    v == "rossler_naive"        && return "naïve (alloc)"
    v == "rossler"              && return "oop (tuned)"
    v == "rossler_naive!"       && return "in-place (naïve)"
    v == "rossler!"             && return "in-place (tuned)"
    v == "rossler_static_naive" && return "SVector (naïve)"
    v == "rossler_static"       && return "SVector (tuned)"
    v == "rossler_type_stable"  && return "type-stable"
    v == "rossler_ad"           && return "AD-ready"
    return v
end

"""Preferred ordering for your case-study variants."""
function preferred_variant_order(df::DataFrame)
    preferred = [
        "rossler_naive",
        "rossler",
        "rossler_naive!",
        "rossler!",
        "rossler_static_naive",
        "rossler_static",
        "rossler_type_stable",
        "rossler_ad",
    ]
    present = Set(String.(df.variant))
    ordered = String[]
    for v in preferred
        v in present && push!(ordered, v)
    end
    # Append any unseen variants (deterministically)
    for v in sort(collect(present))
        v in ordered || push!(ordered, v)
    end
    return ordered
end

"""Ensure required columns exist; add derived time columns in-place."""
function normalize_schema!(df::DataFrame)
    for col in (:solver, :variant, :metric, :median_time_ns)
        hasproperty(df, col) || error("DataFrame missing required column $(col)")
    end
    # Ensure metric is Symbol.
    df.metric = as_symbol.(df.metric)

    # Derived columns (idempotent).
    if !hasproperty(df, :time_ms)
        df.time_ms = df.median_time_ns .* 1e-6
    end
    if !hasproperty(df, :time_s)
        df.time_s = df.median_time_ns .* 1e-9
    end
    return df
end

"""Filter a DataFrame by solver + metric."""
function subset(df::DataFrame; solver::AbstractString, metric::Symbol)
    return filter(r -> (String(r.solver) == solver) && (r.metric == metric), df)
end

# -----------------------------------------------------------------------------
# Plot builders
# -----------------------------------------------------------------------------

"""
    plot_solver_times(df; solver, metric, variants, xscale=:identity)

Horizontal bar chart of median runtimes (milliseconds) for a single solver.
"""
function plot_solver_times(
    df::DataFrame;
    solver::AbstractString,
    metric::Symbol,
    variants::Vector{String},
    xscale::Symbol = :identity,
    title_suffix::AbstractString = "",
)
    d = subset(df; solver=solver, metric=metric)

    # Align into the requested ordering.
    times = Float64[]
    labels = String[]
    for v in variants
        rows = d[d.variant .== v, :]
        nrow(rows) == 0 && continue
        push!(times, rows.time_ms[1])
        push!(labels, pretty_variant(v))
    end

    p = bar(
        reverse(labels),
        reverse(times);
        orientation = :h,
        xlabel = "median time (ms)",
        ylabel = "",
        title = "$(solver): $(String(metric)) times $(title_suffix)",
        legend = false,
        xscale = xscale,
        size = (1100, 650),
        left_margin = 8Plots.mm,
        bottom_margin = 8Plots.mm,
    )
    return p
end

"""
    plot_solver_speedup(df; solver, metric, baseline, variants, xscale=:log10)

Horizontal bar chart of speedups (baseline_time / time). Values > 1 mean faster
than baseline.
"""
function plot_solver_speedup(
    df::DataFrame;
    solver::AbstractString,
    metric::Symbol,
    baseline::AbstractString,
    variants::Vector{String},
    xscale::Symbol = :log10,
)
    d = subset(df; solver=solver, metric=metric)

    base_rows = d[d.variant .== baseline, :]
    if nrow(base_rows) == 0
        # Fall back to the first available variant for this solver+metric.
        baseline = String(first(d.variant))
        base_rows = d[d.variant .== baseline, :]
    end
    base = base_rows.time_ms[1]

    speedup = Float64[]
    labels = String[]
    for v in variants
        rows = d[d.variant .== v, :]
        nrow(rows) == 0 && continue
        push!(speedup, base / rows.time_ms[1])
        push!(labels, pretty_variant(v))
    end

    p = bar(
        reverse(labels),
        reverse(speedup);
        orientation = :h,
        xlabel = "speedup vs $(pretty_variant(String(baseline))) (×)",
        ylabel = "",
        title = "$(solver): speedup (metric = $(String(metric)))",
        legend = false,
        xscale = xscale,
        size = (1100, 650),
        left_margin = 8Plots.mm,
        bottom_margin = 8Plots.mm,
    )
    return p
end

"""Bar chart of allocations per call for a given solver+metric."""
function plot_solver_allocs(
    df::DataFrame;
    solver::AbstractString,
    metric::Symbol,
    variants::Vector{String},
    xscale::Symbol = :identity,
)
    hasproperty(df, :allocs) || error("DataFrame missing required column :allocs")
    d = subset(df; solver=solver, metric=metric)

    allocs = Float64[]
    labels = String[]
    for v in variants
        rows = d[d.variant .== v, :]
        nrow(rows) == 0 && continue
        push!(allocs, float(rows.allocs[1]))
        push!(labels, pretty_variant(v))
    end

    p = bar(
        reverse(labels),
        reverse(allocs);
        orientation = :h,
        xlabel = "allocations per call",
        ylabel = "",
        title = "$(solver): allocations (metric = $(String(metric)))",
        legend = false,
        xscale = xscale,
        size = (1100, 650),
        left_margin = 8Plots.mm,
        bottom_margin = 8Plots.mm,
    )
    return p
end

"""
    plot_pareto(df; metric=:solve)

Scatter plot of (allocations, time_ms), grouped by solver.
This provides a compact, poster-friendly visual: performance vs allocations.
"""
function plot_pareto(df::DataFrame; metric::Symbol = :solve)
    hasproperty(df, :allocs) || error("DataFrame missing required column :allocs")
    d = filter(r -> r.metric == metric, df)

    p = scatter(
        d.allocs,
        d.time_ms;
        group = d.solver,
        xlabel = "allocations per call",
        ylabel = "median time (ms)",
        title = "Pareto view: time vs allocations (metric = $(String(metric)))",
        legend = :topright,
        xscale = :log10,
        yscale = :log10,
        markersize = 6,
        size = (1100, 650),
        left_margin = 8Plots.mm,
        bottom_margin = 8Plots.mm,
    )
    return p
end

# -----------------------------------------------------------------------------
# Main routine
# -----------------------------------------------------------------------------

"""
    make_all_figures(df; outdir="./poster/figures", baseline="rossler_naive")

Create a standard figure suite for the poster and write PNGs into `outdir`.

Outputs (for each solver, with slug `TAG`):
  - `TAG_speedup_solve_log10.png`
  - `TAG_solve_times_ms.png`
  - `TAG_rhs_times_ms.png`
  - `TAG_allocs_rhs.png`
  - `TAG_allocs_solve.png`

Plus overall summaries:
  - `pareto_solve.png`

Additionally, for compatibility with your current `poster.tex`, this function
writes:
  - `rk4_time_normalized_log10.png`
  - `rk4_allocs_solve.png`
when the solver name contains "RK4".
"""
function make_all_figures(
    df::DataFrame;
    outdir::AbstractString = "./poster/figures",
    baseline::AbstractString = "rossler_naive",
)
    normalize_schema!(df)
    mkpath(outdir)

    variants = preferred_variant_order(df)
    solvers = sort(unique(String.(df.solver)))

    # Per-solver figures
    for solver in solvers
        tag = solver_slug(solver)

        # Speedup (solve)
        p_speed = plot_solver_speedup(df; solver=solver, metric=:solve, baseline=baseline, variants=variants, xscale=:log10)
        savefig(p_speed, joinpath(outdir, "$(tag)_speedup_solve_log10.png"))

        # Solve times (ms)
        p_solve = plot_solver_times(df; solver=solver, metric=:solve, variants=variants, xscale=:identity)
        savefig(p_solve, joinpath(outdir, "$(tag)_solve_times_ms.png"))

        # RHS times (ms)
        p_rhs = plot_solver_times(df; solver=solver, metric=:rhs, variants=variants, xscale=:identity)
        savefig(p_rhs, joinpath(outdir, "$(tag)_rhs_times_ms.png"))

        # Allocation charts
        if hasproperty(df, :allocs)
            p_alloc_rhs = plot_solver_allocs(df; solver=solver, metric=:rhs, variants=variants, xscale=:identity)
            savefig(p_alloc_rhs, joinpath(outdir, "$(tag)_allocs_rhs.png"))

            p_alloc_solve = plot_solver_allocs(df; solver=solver, metric=:solve, variants=variants, xscale=:identity)
            savefig(p_alloc_solve, joinpath(outdir, "$(tag)_allocs_solve.png"))
        end

        # Backward-compatible file names referenced by your current poster.
        if occursin(r"(?i)rk4", solver)
            # Reuse the speedup plot as the canonical RK4 figure.
            savefig(p_speed, joinpath(outdir, "rk4_time_normalized_log10.png"))
            if hasproperty(df, :allocs)
                savefig(p_alloc_solve, joinpath(outdir, "rk4_allocs_solve.png"))
            end
        end
    end

    # Global overview
    if hasproperty(df, :allocs)
        p_pareto = plot_pareto(df; metric=:solve)
        savefig(p_pareto, joinpath(outdir, "pareto_solve.png"))
    end

    return outdir
end

end # module PosterFigures

# -----------------------------------------------------------------------------
# Optional CLI wrapper
# -----------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    @info "Running PosterFigures CLI" ARGS

    # We keep the CLI optional to avoid forcing CSV.jl on your project.
    # Usage example (from repo root):
    #   julia --project=. make_poster_figures.jl results.csv ./poster/figures

    if length(ARGS) < 1
        println("Usage: julia make_poster_figures.jl <results.csv> [outdir] [baseline_variant]")
        exit(2)
    end

    csv_path = ARGS[1]
    outdir = length(ARGS) >= 2 ? ARGS[2] : "./poster/figures"
    baseline = length(ARGS) >= 3 ? ARGS[3] : "rossler_naive"

    try
        @eval using CSV
    catch
        error("CSV.jl not available. Either add CSV.jl, or call PosterFigures.make_all_figures(df) from a Julia session.")
    end

    df = CSV.read(csv_path, DataFrame)
    PosterFigures.make_all_figures(df; outdir=outdir, baseline=baseline)
    println("Wrote figures into: $(outdir)")
end
