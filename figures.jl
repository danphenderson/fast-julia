# figures_julia.jl
#
# Poster visual pipeline for Fast Julia / Rössler benchmarks.
#
# Reads the *flattened* benchmark CSV (from `BenchmarkFlatten.dataframe(...)` in utils.jl)
# and generates poster-ready PNGs plus optional LaTeX (booktabs) table snippets.
#
# Example:
#   julia --project=. figures_julia.jl \
#     --csv poster/results.csv \
#     --outdir poster/figures \
#     --tables poster/tables \
#     --poster-compat \
#     --dt-strategy max \
#     --dt-sweep

module PosterFigures

using DataFrames
using Statistics
using Printf
using Dates
using Plots

# CSV is optional unless you pass --csv
# We only do `using CSV` inside main when needed.

# -----------------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------------

const DEFAULT_DPI = 300
const DEFAULT_FIGSIZE = (1200, 820)

# Poster-friendly plot defaults
default(
    dpi = DEFAULT_DPI,
    # Poster readability: prioritize large, legible text.
    titlefontsize = 26,
    guidefontsize = 22,
    tickfontsize = 18,
    legendfontsize = 16,
)

# -----------------------------------------------------------------------------
# Labels / ordering
# -----------------------------------------------------------------------------

as_symbol(x) = x isa Symbol ? x : Symbol(replace(String(x), ":" => ""))

function solver_slug(s::AbstractString)
    t = lowercase(String(s))
    t = replace(t, r"[^a-z0-9]+" => "_")
    t = replace(t, r"^_+|_+$" => "")
    while occursin("__", t)
        t = replace(t, "__" => "_")
    end
    return t
end

function pretty_variant(v::AbstractString)
    v == "rossler_naive"        && return "naïve (alloc)"
    v == "rossler"              && return "oop (tuned)"
    v == "rossler_naive!"       && return "in-place (naïve)"
    v == "rossler!"             && return "in-place (tuned)"
    v == "rossler_static_naive" && return "SVector (naïve)"
    v == "rossler_static"       && return "SVector (tuned)"
    v == "rossler_type_stable"  && return "type-stable"
    v == "rossler_ad"           && return "AD-ready"
    return String(v)
end

function preferred_variant_order(present)::Vector{String}
    preferred = String[
        "rossler_naive",
        "rossler",
        "rossler_naive!",
        "rossler!",
        "rossler_static_naive",
        "rossler_static",
        "rossler_type_stable",
        "rossler_ad",
    ]
    present_set = Set(String.(present))
    ordered = String[v for v in preferred if v in present_set]
    for v in sort(collect(present_set))
        v in ordered || push!(ordered, v)
    end
    return ordered
end

# -----------------------------------------------------------------------------
# Schema normalization
# -----------------------------------------------------------------------------

function _ensure_col(df::DataFrame, col::Symbol)
    hasproperty(df, col) || error("CSV missing required column :$(col). Columns: $(names(df))")
end

"""Normalize expected schema in-place.

Required columns:
  :solver, :variant, :metric, :median_time_ns
Optional columns:
  :dt, :nsteps, :memory, :allocs

Adds:
  :time_ms, :time_s

Also normalizes metric to Symbol and strips a leading ':' if present.
"""
function normalize_schema!(df::DataFrame)
    for col in (:solver, :variant, :metric, :median_time_ns)
        _ensure_col(df, col)
    end

    df.metric = as_symbol.(df.metric)

    # Ensure optional columns exist
    for col in (:dt, :nsteps, :memory, :allocs)
        if !hasproperty(df, col)
            df[!, col] = Vector{Union{Missing, Float64}}(missing, nrow(df))
        end
    end

    # Derived columns
    if !hasproperty(df, :time_ms)
        df.time_ms = Float64.(df.median_time_ns) .* 1e-6
    end
    if !hasproperty(df, :time_s)
        df.time_s = Float64.(df.median_time_ns) .* 1e-9
    end

    return df
end

# -----------------------------------------------------------------------------
# dt selection for solve rows (one row per variant)
# -----------------------------------------------------------------------------

"""Select a deterministic dt slice for solve metric.

The experiment produces multiple dt rows per (solver, variant). For bar charts and
LaTeX tables we must choose a single row per variant.

Strategies:
  - :max     => choose the largest dt per variant (typically the first dt in your dt_schedule)
  - :min     => choose the smallest dt per variant
  - :nearest => choose dt closest to dt_target per variant

If dt is missing for all rows, we fall back to choosing the row with minimum median_time_ns per variant.
"""
function select_solve_slice(
    df::DataFrame;
    solver::AbstractString,
    dt_strategy::Symbol = :max,
    dt_target::Union{Nothing,Float64} = nothing,
)
    d = filter(r -> (String(r.solver) == solver) && (r.metric == :solve), df)
    nrow(d) == 0 && return d

    # If dt is entirely missing, pick best-effort single row per variant
    all_dt_missing = all(ismissing, d.dt)
    if all_dt_missing
        g = groupby(d, :variant)
        rows = DataFrame()
        idxs = Int[]
        for sub in g
            # pick the minimum-time row deterministically
            i = argmin(sub.median_time_ns)
            push!(idxs, parentindices(sub)[1][i])
        end
        return d[idxs, :]
    end

    if dt_strategy == :nearest
        dt_target === nothing && error("dt_strategy=:nearest requires a dt_target")
    end

    g = groupby(d, :variant)
    idxs = Int[]
    for sub in g
        # sub.dt is Union{Missing, T}. Filter missing.
        dtvals = collect(skipmissing(sub.dt))
        if isempty(dtvals)
            i = argmin(sub.median_time_ns)
            push!(idxs, parentindices(sub)[1][i])
            continue
        end

        # choose index within `sub`
        chosen_i = 1
        if dt_strategy == :max
            chosen_i = argmax(collect(coalesce.(sub.dt, -Inf)))
        elseif dt_strategy == :min
            chosen_i = argmin(collect(coalesce.(sub.dt, Inf)))
        elseif dt_strategy == :nearest
            # compute abs(dt - target)
            diffs = map(x -> ismissing(x) ? Inf : abs(Float64(x) - dt_target), sub.dt)
            chosen_i = argmin(diffs)
        else
            error("Unknown dt_strategy=$(dt_strategy). Use :max, :min, or :nearest")
        end

        push!(idxs, parentindices(sub)[1][chosen_i])
    end

    return d[idxs, :]
end

# -----------------------------------------------------------------------------
# Plot helpers
# -----------------------------------------------------------------------------

function _hbar(labels::Vector{String}, values::Vector{<:Real};
    xlabel::AbstractString,
    title::AbstractString,
    xscale::Symbol = :identity,
)
    y = collect(1:length(labels))
    bar(
        y,
        Float64.(values);
        orientation = :h,
        xlabel = xlabel,
        ylabel = "",
        yticks = (y, labels),
        yflip = true,
        title = title,
        legend = false,
        xscale = xscale,
        size = DEFAULT_FIGSIZE,
        left_margin = 20Plots.mm,
        bottom_margin = 12Plots.mm,
        grid = :x,
    )
end

function plot_solver_times_ms(
    df_slice::DataFrame;
    solver::AbstractString,
    variants::Vector{String},
    metric_label::AbstractString,
    title_suffix::AbstractString = "",
    xscale::Symbol = :identity,
)
    labels = String[]
    times = Float64[]
    for v in variants
        rows = df_slice[df_slice.variant .== v, :]
        nrow(rows) == 0 && continue
        push!(labels, pretty_variant(v))
        push!(times, Float64(rows.time_ms[1]))
    end

    # poster ordering: reverse so last is at top
    labels = reverse(labels)
    times  = reverse(times)

    return _hbar(labels, times;
        xlabel = "median time (ms)",
        title = "$(solver): $(metric_label) times $(title_suffix)" |> strip,
        xscale = xscale,
    )
end

function plot_solver_speedup(
    df_slice::DataFrame;
    solver::AbstractString,
    variants::Vector{String},
    baseline::AbstractString,
    metric_label::AbstractString,
    xscale::Symbol = :identity,
)
    # baseline preference: rossler_naive
    base_rows = df_slice[df_slice.variant .== baseline, :]
    if nrow(base_rows) == 0 && baseline != "rossler_naive"
        base_rows = df_slice[df_slice.variant .== "rossler_naive", :]
        baseline = "rossler_naive"
    end
    nrow(base_rows) == 0 && error("Baseline $(baseline) missing for solver=$(solver)")

    base = Float64(base_rows.median_time_ns[1])

    labels = String[]
    speed  = Float64[]
    for v in variants
        rows = df_slice[df_slice.variant .== v, :]
        nrow(rows) == 0 && continue
        push!(labels, pretty_variant(v))
        push!(speed, base / Float64(rows.median_time_ns[1]))
    end

    labels = reverse(labels)
    speed  = reverse(speed)

    return _hbar(labels, speed;
        xlabel = "speedup vs baseline ($(baseline)) (×)",
        title = "$(solver): speedup ($(metric_label))",
        xscale = xscale,
    )
end

function plot_solver_allocs(
    df_slice::DataFrame;
    solver::AbstractString,
    variants::Vector{String},
    metric_label::AbstractString,
    xscale::Symbol = :identity,
)
    hasproperty(df_slice, :allocs) || error("DataFrame has no :allocs column")

    labels = String[]
    allocs = Float64[]
    for v in variants
        rows = df_slice[df_slice.variant .== v, :]
        nrow(rows) == 0 && continue
        push!(labels, pretty_variant(v))
        push!(allocs, Float64(rows.allocs[1]))
    end

    labels = reverse(labels)
    allocs = reverse(allocs)

    unit = metric_label == "solve" ? "solve" : "call"

    xlabel = "allocations per $(unit)"
    values = allocs
    xscale_plot = :identity
    if xscale == :log10
        # Plots.jl horizontal bar + xscale=:log10 can render poorly; instead plot log10-values explicitly.
        values = log10.(allocs)
        xlabel = "log10 allocations per $(unit)"
    else
        xscale_plot = xscale
    end

    return _hbar(labels, values;
        xlabel = xlabel,
        title = "$(solver): allocations ($(metric_label))",
        xscale = xscale_plot,
    )
end

function plot_pareto(df::DataFrame; metric::Symbol = :solve)
    hasproperty(df, :allocs) || error("DataFrame has no :allocs column")

    d = filter(r -> r.metric == metric, df)
    d = filter(r -> !(ismissing(r.allocs) || ismissing(r.time_ms)), d)
    nrow(d) == 0 && error("No rows available for pareto plot (metric=$(metric))")

    scatter(
        Float64.(d.allocs),
        Float64.(d.time_ms);
        group = String.(d.solver),
        xlabel = "allocations per call",
        ylabel = "median time (ms)",
        title = "Pareto view: time vs allocations (metric=$(String(metric)))",
        xscale = :log10,
        yscale = :log10,
        legend = :topright,
        markersize = 6,
        size = DEFAULT_FIGSIZE,
        left_margin = 8Plots.mm,
        bottom_margin = 8Plots.mm,
    )
end

function plot_dt_sweep(
    df::DataFrame;
    solver::AbstractString,
    variants::Vector{String},
    x::Symbol = :nsteps,
    y::Symbol = :time_ms,
    loglog::Bool = true,
)
    d = filter(r -> (String(r.solver) == solver) && (r.metric == :solve), df)
    nrow(d) == 0 && return plot(title="$(solver): solve dt sweep (no data)")

    # drop missing y/dt
    d = filter(r -> !(ismissing(r.dt) || ismissing(getproperty(r, y))), d)
    nrow(d) == 0 && return plot(title="$(solver): solve dt sweep (no dt/y)")

    if x == :nsteps
        # If nsteps missing for all rows, fallback to dt
        if !hasproperty(d, :nsteps) || all(ismissing, d.nsteps)
            x = :dt
        end
    end

    p = plot(
        xlabel = String(x),
        ylabel = (y == :time_ms ? "median time (ms)" : String(y)),
        title = "$(solver): solve cost vs $(String(x))",
        legend = :best,
        size = DEFAULT_FIGSIZE,
        left_margin = 8Plots.mm,
        bottom_margin = 8Plots.mm,
    )

    for v in variants
        g = d[d.variant .== v, :]
        nrow(g) == 0 && continue
        # Sort by dt descending (matches Python default)
        if hasproperty(g, :dt)
            sort!(g, :dt, rev=true)
        end
        xs = Float64.(coalesce.(g[!, x], NaN))
        ys = Float64.(g[!, y])
        plot!(p, xs, ys; marker=:circle, label=pretty_variant(v))
    end

    if loglog
        # Plots.jl sets axis scaling via keywords (no universal `xscale!`/`yscale!`).
        plot!(p; xscale=:log10, yscale=:log10)
    end

    return p
end

# -----------------------------------------------------------------------------
# LaTeX tables
# -----------------------------------------------------------------------------

function latex_escape(s::AbstractString)
    # Minimal escaping sufficient for these variant labels.
    s = replace(String(s), "\\" => "\\textbackslash{}")
    s = replace(s, "_" => "\\_")
    s = replace(s, "%" => "\\%")
    s = replace(s, "&" => "\\&")
    s = replace(s, "#" => "\\#")
    return s
end

function _fmt_float(x; sig=3)
    if x isa Missing
        return "--"
    end
    if x isa AbstractFloat
        isfinite(x) || return "--"
        return @sprintf("%.*g", sig, Float64(x))
    end
    # ints etc
    return string(x)
end

"""Create a booktabs tabular snippet (no table float wrapper)."""
function make_latex_table(
    df_slice::DataFrame;
    solver::AbstractString,
    metric::Symbol,
    variants::Vector{String},
)
    cols = ["Variant"]
    if metric == :rhs
        push!(cols, "Median (ns)")
    else
        push!(cols, "Median (ms)")
    end

    has_allocs = hasproperty(df_slice, :allocs) && !all(ismissing, df_slice.allocs)
    has_mem    = hasproperty(df_slice, :memory) && !all(ismissing, df_slice.memory)
    has_dt     = (metric == :solve) && hasproperty(df_slice, :dt) && !all(ismissing, df_slice.dt)
    has_steps  = (metric == :solve) && hasproperty(df_slice, :nsteps) && !all(ismissing, df_slice.nsteps)

    if has_dt
        insert!(cols, 2, "dt")
    end
    if has_steps
        insert!(cols, has_dt ? 3 : 2, "nsteps")
    end
    if has_allocs
        push!(cols, metric == :rhs ? "Allocs/call" : "Allocs/solve")
        push!(cols, metric == :rhs ? "Bytes/call" : "Bytes/solve")
    elseif has_mem
        push!(cols, metric == :rhs ? "Bytes/call" : "Bytes/solve")
    end

    io = IOBuffer()

    # column spec: l then r...
    print(io, "\\begin{tabular}{l")
    for _ in 2:length(cols)
        print(io, "r")
    end
    println(io, "}")
    println(io, "\\toprule")
    println(io, join(cols, " & ") * " \\\\")
    println(io, "\\midrule")

    for v in variants
        rows = df_slice[df_slice.variant .== v, :]
        nrow(rows) == 0 && continue
        r0 = rows[1, :]

        parts = String[]
        push!(parts, latex_escape(pretty_variant(v)))

        if has_dt
            push!(parts, _fmt_float(r0.dt))
        end
        if has_steps
            ns = r0.nsteps
            push!(parts, ns isa Missing ? "--" : string(Int(ns)))
        end

        if metric == :rhs
            push!(parts, _fmt_float(r0.median_time_ns))
        else
            push!(parts, _fmt_float(r0.time_ms))
        end

        if has_allocs
            a = r0.allocs
            m = r0.memory
            push!(parts, a isa Missing ? "--" : string(Int(a)))
            push!(parts, m isa Missing ? "--" : string(Int(m)))
        elseif has_mem
            m = r0.memory
            push!(parts, m isa Missing ? "--" : string(Int(m)))
        end

        println(io, join(parts, " & ") * " \\\\")
    end

    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")

    return String(take!(io))
end

# -----------------------------------------------------------------------------
# Orchestration
# -----------------------------------------------------------------------------

function make_all_figures(
    df::DataFrame;
    outdir::AbstractString = "./poster/figures",
    tables_dir::Union{Nothing,AbstractString} = nothing,
    baseline::AbstractString = "rossler_naive",
    dt_strategy::Symbol = :max,
    dt_target::Union{Nothing,Float64} = nothing,
    poster_compat::Bool = false,
    include_dt_sweep::Bool = false,
)
    normalize_schema!(df)
    mkpath(outdir)
    if tables_dir !== nothing
        mkpath(tables_dir)
    end

    solvers = sort(unique(String.(df.solver)))

    for solver in solvers
        tag = solver_slug(solver)
        present_variants = unique(String.(df[df.solver .== solver, :variant]))
        variants = preferred_variant_order(present_variants)

        rhs   = filter(r -> (String(r.solver) == solver) && (r.metric == :rhs), df)
        solve = select_solve_slice(df; solver=solver, dt_strategy=dt_strategy, dt_target=dt_target)

        # Figures
        if nrow(solve) > 0
            p_speed = plot_solver_speedup(solve; solver=solver, variants=variants, baseline=baseline, metric_label="solve", xscale=:identity)
            savefig(p_speed, joinpath(outdir, "$(tag)_speedup_solve.png"))

            p_solve = plot_solver_times_ms(solve; solver=solver, variants=variants, metric_label="solve")
            savefig(p_solve, joinpath(outdir, "$(tag)_solve_times_ms.png"))

            if hasproperty(solve, :allocs) && !all(ismissing, solve.allocs)
                # log scale makes low-allocation variants visible next to 10^5-scale baselines
                p_alloc_solve = plot_solver_allocs(solve; solver=solver, variants=variants, metric_label="solve", xscale=:log10)
                savefig(p_alloc_solve, joinpath(outdir, "$(tag)_allocs_solve.png"))
            end
        end

        if nrow(rhs) > 0
            p_rhs = plot_solver_times_ms(rhs; solver=solver, variants=variants, metric_label="rhs")
            savefig(p_rhs, joinpath(outdir, "$(tag)_rhs_times_ms.png"))

            if hasproperty(rhs, :allocs) && !all(ismissing, rhs.allocs)
                p_alloc_rhs = plot_solver_allocs(rhs; solver=solver, variants=variants, metric_label="rhs")
                savefig(p_alloc_rhs, joinpath(outdir, "$(tag)_allocs_rhs.png"))
            end
        end

        if include_dt_sweep
            # dt sweep uses the full solve rows (not the dt slice)
            dsolve = filter(r -> (String(r.solver) == solver) && (r.metric == :solve), df)
            if nrow(dsolve) > 0 && hasproperty(dsolve, :dt) && !all(ismissing, dsolve.dt)
                p_sweep = plot_dt_sweep(df; solver=solver, variants=variants, x=:nsteps, y=:time_ms, loglog=true)
                savefig(p_sweep, joinpath(outdir, "$(tag)_solve_dt_sweep_loglog.png"))
            end
        end

        # Tables
        if tables_dir !== nothing
            if nrow(rhs) > 0
                tex = make_latex_table(rhs; solver=solver, metric=:rhs, variants=variants)
                open(joinpath(tables_dir, "$(tag)_rhs_table.tex"), "w") do io
                    write(io, tex)
                end
            end
            if nrow(solve) > 0
                tex = make_latex_table(solve; solver=solver, metric=:solve, variants=variants)
                open(joinpath(tables_dir, "$(tag)_solve_table.tex"), "w") do io
                    write(io, tex)
                end
            end
        end

        # Poster compatibility filenames
        if poster_compat && occursin(r"(?i)rk4", solver)
            src = joinpath(outdir, "$(tag)_speedup_solve.png")
            isfile(src) && cp(src, joinpath(outdir, "rk4_time_normalized_log10.png"); force=true)

            src2 = joinpath(outdir, "$(tag)_allocs_solve.png")
            isfile(src2) && cp(src2, joinpath(outdir, "rk4_allocs_solve.png"); force=true)
        end
    end

    # Global Pareto (uses all solve rows)
    if hasproperty(df, :allocs) && !all(ismissing, df.allocs)
        p_pareto = plot_pareto(df; metric=:solve)
        savefig(p_pareto, joinpath(outdir, "pareto_solve.png"))
    end

    return outdir
end

# -----------------------------------------------------------------------------
# Minimal CLI parsing
# -----------------------------------------------------------------------------

struct CLIArgs
    csv::String
    outdir::String
    tables::String
    baseline::String
    dt_strategy::Symbol
    dt_target::Union{Nothing,Float64}
    poster_compat::Bool
    dt_sweep::Bool
end

function _print_help()
    println("figures_julia.jl - poster visual pipeline")
    println("\nOptions:")
    println("  --csv PATH           (required) flattened CSV")
    println("  --outdir DIR         output directory for PNGs (default: ./poster/figures)")
    println("  --tables DIR         if set, write LaTeX table snippets to this directory")
    println("  --baseline VAR       baseline variant for speedup (default: rossler_naive)")
    println("  --dt-strategy STR    max|min|nearest (default: max)")
    println("  --dt-target FLOAT    required when dt-strategy=nearest")
    println("  --poster-compat      also write poster-compat filenames")
    println("  --dt-sweep           also generate dt-sweep plots")
    println("  --help               show this help")
end

function parse_cli(args::Vector{String})::CLIArgs
    csv = ""
    outdir = "./poster/figures"
    tables = ""
    baseline = "rossler_naive"
    dt_strategy = :max
    dt_target::Union{Nothing,Float64} = nothing
    poster_compat = false
    dt_sweep = false

    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--help" || a == "-h"
            _print_help()
            exit(0)
        elseif a == "--csv"
            i += 1; i <= length(args) || error("--csv requires a value")
            csv = args[i]
        elseif a == "--outdir"
            i += 1; i <= length(args) || error("--outdir requires a value")
            outdir = args[i]
        elseif a == "--tables"
            i += 1; i <= length(args) || error("--tables requires a value")
            tables = args[i]
        elseif a == "--baseline"
            i += 1; i <= length(args) || error("--baseline requires a value")
            baseline = args[i]
        elseif a == "--dt-strategy"
            i += 1; i <= length(args) || error("--dt-strategy requires a value")
            s = lowercase(args[i])
            s in ("max", "min", "nearest") || error("--dt-strategy must be max|min|nearest")
            dt_strategy = Symbol(s)
        elseif a == "--dt-target"
            i += 1; i <= length(args) || error("--dt-target requires a value")
            dt_target = parse(Float64, args[i])
        elseif a == "--poster-compat"
            poster_compat = true
        elseif a == "--dt-sweep"
            dt_sweep = true
        else
            error("Unknown argument: $(a). Use --help.")
        end
        i += 1
    end

    isempty(csv) && error("--csv is required (path to flattened CSV). Use --help.")
    if dt_strategy == :nearest && dt_target === nothing
        error("--dt-strategy nearest requires --dt-target")
    end

    return CLIArgs(csv, outdir, tables, baseline, dt_strategy, dt_target, poster_compat, dt_sweep)
end

end # module PosterFigures

# -----------------------------------------------------------------------------
# Script entrypoint
# -----------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    args = PosterFigures.parse_cli(copy(ARGS))

    # Load CSV lazily
    @eval using CSV
    @eval using DataFrames

    df = CSV.read(args.csv, DataFrames.DataFrame)
    PosterFigures.normalize_schema!(df)

    tables_dir = isempty(args.tables) ? nothing : args.tables

    PosterFigures.make_all_figures(
        df;
        outdir = args.outdir,
        tables_dir = tables_dir,
        baseline = args.baseline,
        dt_strategy = args.dt_strategy,
        dt_target = args.dt_target,
        poster_compat = args.poster_compat,
        include_dt_sweep = args.dt_sweep,
    )

    println("Wrote figures into: $(args.outdir)")
    if tables_dir !== nothing
        println("Wrote LaTeX tables into: $(tables_dir)")
    end
end
