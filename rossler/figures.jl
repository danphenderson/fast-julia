# rossler/make_poster_figures.jl

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

as_symbol(x) = x isa Symbol ? x : Symbol(x)

function solver_slug(s::AbstractString)
    t = lowercase(String(s))
    t = replace(t, r"[^a-z0-9]+" => "_")
    t = replace(t, r"^_+|_+$" => "")
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
    return v
end

function preferred_variant_order(df::DataFrame)
    preferred = [
        "rossler_naive",
        "rossler",
        "rossler_naive!",
        "rossler!",
        "rossler_static_naive",
        "rossler_static",
    ]
    present = Set(String.(df.variant))
    ordered = String[]
    for v in preferred
        v in present && push!(ordered, v)
    end
    for v in sort(collect(present))
        v in ordered || push!(ordered, v)
    end
    return ordered
end

function normalize_schema!(df::DataFrame)
    for col in (:solver, :variant, :metric, :median_time_ns)
        hasproperty(df, col) || error("DataFrame missing required column $(col)")
    end
    df.metric = as_symbol.(df.metric)

    if !hasproperty(df, :time_ms)
        df.time_ms = df.median_time_ns .* 1e-6
    end
    if !hasproperty(df, :time_s)
        df.time_s = df.median_time_ns .* 1e-9
    end
    return df
end

function subset(df::DataFrame; solver::AbstractString, metric::Symbol)
    return filter(r -> (String(r.solver) == solver) && (r.metric == metric), df)
end

# -----------------------------------------------------------------------------
# Plot builders
# -----------------------------------------------------------------------------

function plot_solver_times(
    df::DataFrame;
    solver::AbstractString,
    metric::Symbol,
    variants::Vector{String},
    xscale::Symbol = :identity,
    title_suffix::AbstractString = "",
)
    d = subset(df; solver=solver, metric=metric)

    times = Float64[]
    labels = String[]
    for v in variants
        rows = d[d.variant .== v, :]
        nrow(rows) == 0 && continue
        push!(times, rows.time_ms[1])
        push!(labels, pretty_variant(v))
    end

    return bar(
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
end

function plot_solver_speedup(
    df::DataFrame;
    solver::AbstractString,
    metric::Symbol,
    baseline::AbstractString,
    variants::Vector{String},
    xscale::Symbol = :log10,
)
    d = subset(df; solver=solver, metric=metric)

    base_rows = d[d.variant .== "rossler_naive", :]
    if nrow(base_rows) == 0
        base_rows = d[d.variant .== baseline, :]
    end
    nrow(base_rows) == 0 && error("Baseline variant $(baseline) missing for solver $(solver)")
    base = base_rows.median_time_ns[1]

    speedup = Float64[]
    labels = String[]
    for v in variants
        rows = d[d.variant .== v, :]
        nrow(rows) == 0 && continue
        push!(speedup, base / rows.median_time_ns[1])
        push!(labels, pretty_variant(v))
    end

    return bar(
        reverse(labels),
        reverse(speedup);
        orientation = :h,
        xlabel = "speedup vs naïve baseline (rossler_naive) (×)",
        ylabel = "",
        title = "$(solver): speedup (metric = $(String(metric)))",
        legend = false,
        xscale = xscale,
        size = (1100, 650),
        left_margin = 8Plots.mm,
        bottom_margin = 8Plots.mm,
    )
end

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

    return bar(
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
end

function plot_pareto(df::DataFrame; metric::Symbol = :solve)
    hasproperty(df, :allocs) || error("DataFrame missing required column :allocs")
    d = filter(r -> r.metric == metric, df)

    return scatter(
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
end

# -----------------------------------------------------------------------------
# Main routine (figure suite)
# -----------------------------------------------------------------------------

function make_all_figures(
    df::DataFrame;
    outdir::AbstractString = "./poster/figures",
    baseline::AbstractString = "rossler_naive",
    variants::Vector{String} = preferred_variant_order(df),
)
    normalize_schema!(df)
    mkpath(outdir)
    solvers = sort(unique(String.(df.solver)))

    for solver in solvers
        tag = solver_slug(solver)

        p_speed = plot_solver_speedup(df; solver=solver, metric=:solve, baseline=baseline, variants=variants, xscale=:log10)
        savefig(p_speed, joinpath(outdir, "$(tag)_speedup_solve_log10.png"))

        p_solve = plot_solver_times(df; solver=solver, metric=:solve, variants=variants)
        savefig(p_solve, joinpath(outdir, "$(tag)_solve_times_ms.png"))

        p_rhs = plot_solver_times(df; solver=solver, metric=:rhs, variants=variants)
        savefig(p_rhs, joinpath(outdir, "$(tag)_rhs_times_ms.png"))

        if hasproperty(df, :allocs)
            p_alloc_rhs = plot_solver_allocs(df; solver=solver, metric=:rhs, variants=variants)
            savefig(p_alloc_rhs, joinpath(outdir, "$(tag)_allocs_rhs.png"))

            p_alloc_solve = plot_solver_allocs(df; solver=solver, metric=:solve, variants=variants)
            savefig(p_alloc_solve, joinpath(outdir, "$(tag)_allocs_solve.png"))
        end

        # Poster compatibility filenames
        if occursin(r"(?i)rk4", solver)
            savefig(p_speed, joinpath(outdir, "rk4_time_normalized_log10.png"))
            if hasproperty(df, :allocs)
                savefig(p_alloc_solve, joinpath(outdir, "rk4_allocs_solve.png"))
            end
        end
    end

    if hasproperty(df, :allocs)
        p_pareto = plot_pareto(df; metric=:solve)
        savefig(p_pareto, joinpath(outdir, "pareto_solve.png"))
    end

    return outdir
end

end # module PosterFigures

# -----------------------------------------------------------------------------
# Convenience entrypoint for REPL usage: include(...); main()
# -----------------------------------------------------------------------------
function main(; outdir::AbstractString = "./poster/figures", baseline::AbstractString = "rossler_naive")
    rossler_dir = @__DIR__

    if !isdefined(Main, :run_studies)
        include(joinpath(rossler_dir, "benchmark.jl"))
    end
    if !isdefined(Main, :BenchmarkFlatten)
        include(joinpath(rossler_dir, "flatten_results.jl"))
    end

    results = run_studies()
    df = BenchmarkFlatten.dataframe(results)
    return PosterFigures.make_all_figures(df; outdir=outdir, baseline=baseline)
end

# -----------------------------------------------------------------------------
# Script entrypoint:
#   julia --project=. rossler/make_poster_figures.jl [outdir] [baseline] [results.csv]
# -----------------------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    outdir   = length(ARGS) >= 1 ? ARGS[1] : "./poster/figures"
    baseline = length(ARGS) >= 2 ? ARGS[2] : "rossler_naive"
    csv_path = length(ARGS) >= 3 ? ARGS[3] : ""

    if !isempty(csv_path)
        try
            @eval using CSV
        catch
            error("CSV.jl not available. Either add CSV.jl or run without a CSV path.")
        end
        df = CSV.read(csv_path, DataFrames.DataFrame)
        PosterFigures.make_all_figures(df; outdir=outdir, baseline=baseline)
    else
        main(; outdir=outdir, baseline=baseline)
    end

    println("Wrote figures into: $(outdir)")
end
