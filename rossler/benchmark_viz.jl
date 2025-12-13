module BenchmarkViz

using Statistics
using BenchmarkTools
using DataFrames
using Plots
using StatsPlots

"""
    aggregate_median(flat::AbstractVector{<:NamedTuple}) -> DataFrame
]
Compute summary statistics for each combination of solver, variant and metric.
The input `flat` is expected to be the flattened benchmark output, e.g. as
returned by `BenchmarkFlatten.flatten_run_studies(results)`.  The function
groups by `:solver`, `:variant` and `:metric` and returns a `DataFrame` with
columns:

* `solver`, `variant`, `metric` – grouping keys.
* `median_time_ns` – the median runtime in nanoseconds.
* `median_time_ms` – median runtime in milliseconds.
* `mean_time_ns` – the mean runtime in nanoseconds.
* `q25_ns`, `q75_ns` – the 25th and 75th percentiles of the runtime.
* `memory` – bytes allocated per call.
* `allocs` – allocations per call.

These statistics can be used to create bar plots with error bars representing
interquartile ranges.
"""
function aggregate_median(flat)
    # Convert to DataFrame for grouping
    df = DataFrame(flat)
    # Compute additional columns: median ms
    df.median_time_ms = df.median_time_ns .* 1e-6
    # Group and aggregate
    g = groupby(df, [:solver, :variant, :metric])
    agg = combine(g, [:median_time_ns] => median => :median_time_ns,
                     [:median_time_ns] => mean => :mean_time_ns,
                     [:median_time_ns] => x -> quantile(x, 0.25) => :q25_ns,
                     [:median_time_ns] => x -> quantile(x, 0.75) => :q75_ns,
                     [:memory] => first => :memory,
                     [:allocs] => first => :allocs)
    # Compute ms columns
    agg.median_time_ms = agg.median_time_ns .* 1e-6
    return agg
end

"""
    plot_solve_times(agg::DataFrame; bar_width=0.6, palette=nothing)

Create a grouped bar chart of median solve times for each solver and variant.
The input `agg` should be a `DataFrame` produced by `aggregate_median` and
filtered to rows where `metric == :solve`.  The `bar_width` controls the
fractional width of the bars, and `palette` may be a list of colors or a
ColorScheme; if omitted, the default palette is used.  The function
returns a `Plot` object that can be further customised or saved with
`savefig`.
"""
function plot_solve_times(agg; bar_width=0.6, palette=nothing)
    df = filter(row -> row.metric == :solve, agg)
    # Sort variants for consistent ordering
    unique_variants = sort(unique(df.variant))
    # Build a matrix of values: rows = variants, cols = solvers
    solvers = sort(unique(df.solver))
    data = [first(filter(r -> r.solver == s && r.variant == v, df)).median_time_ms for v in unique_variants, s in solvers]
    p = groupedbar(unique_variants, data;
                   group = solvers,
                   bar_position = :dodge,
                   bar_width = bar_width,
                   xlabel = "Variant",
                   ylabel = "Median solve time (ms)",
                   legend = :topright,
                   title = "Median solve times by solver and variant")
    if palette !== nothing
        p.seriesargs[:color] = palette
    end
    return p
end

"""
    compute_speedups(agg::DataFrame; baseline::String) -> DataFrame

Compute speedup factors relative to a baseline variant for each solver.  The
input `agg` must be filtered to `metric == :solve`.  For each solver, the
median runtime of the `baseline` variant is used as the reference.  The
returned `DataFrame` contains `solver`, `variant`, `median_time_ns` and
`speedup` where `speedup = baseline_time / median_time_ns`.
"""
function compute_speedups(agg; baseline::String)
    df = filter(row -> row.metric == :solve, agg)
    # Build a dictionary of baseline times per solver
    base_times = Dict(row.solver => row.median_time_ns for row in eachrow(df) if row.variant == baseline)
    rows = NamedTuple[]
    for row in eachrow(df)
        ref = get(base_times, row.solver, missing)
        if ref === missing
            continue
        end
        push!(rows, (solver = row.solver, variant = row.variant,
                     median_time_ns = row.median_time_ns,
                     speedup = ref / row.median_time_ns))
    end
    return DataFrame(rows)
end

"""
    plot_speedups(speed::DataFrame; bar_width=0.6, palette=nothing)

Plot bar charts of speedup factors computed with `compute_speedups`.  The
`speed` DataFrame should have columns `:solver`, `:variant` and `:speedup`.
Bars are grouped by solver and displayed for each variant.  A legend and
axes labels are added automatically.
"""
function plot_speedups(speed; bar_width=0.6, palette=nothing)
    # Determine variants and solvers
    variants = sort(unique(speed.variant))
    solvers = sort(unique(speed.solver))
    # Build matrix of speedup values
    values = [first(filter(r -> r.solver == s && r.variant == v, speed)).speedup for v in variants, s in solvers]
    p = groupedbar(variants, values;
                   group = solvers,
                   bar_position = :dodge,
                   bar_width = bar_width,
                   xlabel = "Variant",
                   ylabel = "Speedup (× baseline)",
                   legend = :topright,
                   title = "Speedup relative to baseline variant")
    if palette !== nothing
        p.seriesargs[:color] = palette
    end
    return p
end

"""
    plot_rhs_violin(results::Dict)

Create violin plots of right–hand-side micro-benchmark timing distributions for
each variant and solver.  The `results` argument should be the dictionary
returned by `run_studies()`.  This function extracts the `times` vector from
each `BenchmarkTools.Trial` in `rhs_trials` and plots a violin plot grouped
by `(solver, variant)`.  It returns the plot; call `savefig` to write
the figure to disk.
"""
function plot_rhs_violin(results)
    # Build list of entries: solver, variant, times (ms)
    solver_list = String[]
    variant_list = String[]
    times_list = Vector{Float64}[]
    for (solver, study) in results
        for (fn, trial) in study.rhs_trials
            push!(solver_list, solver)
            push!(variant_list, fn)
            # Convert to ms
            push!(times_list, trial.times .* 1e-6)
        end
    end
    # Create a categorical label combining solver and variant for grouping
    labels = string.(solver_list, " - ", variant_list)
    p = violin(labels, times_list;
               xlabel = "Solver and variant",
               ylabel = "RHS call time (ms)",
               title = "Distribution of RHS call times",
               side = :both,
               linewidth = 0)
    return p

end
function plot_each_figure(results)
    for (solver, study) in results
        for (fn, trial) in study.rhs_trials
            p = plot(trial.times .* 1e-6;
                     xlabel = "Time (ms)",
                     ylabel = "Frequency",
                     title = "RHS call time distribution for $solver - $fn",
                     legend = false)
            savefig(p, "rhs_$solver-$fn.png")
        end
    end
end
function plot_rhs_histograms(results)
    for (solver, study) in results
        for (fn, trial) in study.rhs_trials
            p = histogram(trial.times .* 1e-6;
                          xlabel = "Time (ms)",
                          ylabel = "Frequency",
                          title = "RHS call time distribution for $solver - $fn",
                          legend = false)
            savefig(p, "rhs_$solver-$fn.png")
        end
    end
end

end # module BenchmarkViz