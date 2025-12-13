#!/usr/bin/env julia
#= 
Generate poster figures and tables from benchmark data.

This script expects to run inside the project environment with the
following packages available: JSON3, DataFrames, Statistics, Plots,
StatsPlots, StaticArrays, DifferentialEquations, and OrdinaryDiffEq.
It reads benchmark metrics from `poster/data.json`, produces summary
figures under `poster/figures/`, and writes a LaTeX table of solve times.
Run with:
    julia --project -e 'include("poster/make_poster_figures.jl")'
=#

using JSON3
using DataFrames
using Statistics
using Plots
using StatsPlots
using StaticArrays
using DifferentialEquations
using OrdinaryDiffEq

const DATA_PATH = joinpath(@__DIR__, "data.json")
const FIG_DIR = joinpath(@__DIR__, "figures")
const DEFAULT_SOLVER = "Midpoint Fixed"  # switch to "RK4 Fixed" if desired
const TRAJECTORY_VARIANT = "rossler!"

include(joinpath(@__DIR__, "..", "rossler", "impl.jl"))

mkpath(FIG_DIR)

function load_benchmark_data(path::AbstractString)
    raw = JSON3.read(read(path, String))
    return DataFrame(raw)
end

solver_slug(solver::AbstractString) = lowercase(first(split(solver)))

function plot_solve_times(df; solver::AbstractString = DEFAULT_SOLVER, logscale::Bool = false)
    filtered = filter(row -> row.solver == solver && row.metric == "solve", df)
    isempty(filtered) && error("No solve-time rows found for solver '$solver'")
    sort!(filtered, :median_time_ns)
    plt = bar(
        filtered.variant,
        filtered.median_time_ns,
        legend=false,
        xlabel="Variant",
        ylabel="Median solve time (ns)",
        yscale = logscale ? :log10 : :identity,
        xtickfont=font(7),
        rotation=25,
        title="$(solver) solve times",
    )
    return plt
end

function aggregate_rhs(df)
    rhs_rows = filter(:metric => ==("rhs"), df)
    isempty(rhs_rows) && error("No RHS rows found in benchmark data")
    return combine(groupby(rhs_rows, :variant),
        :median_time_ns => median => :median_time_ns,
        :memory => median => :memory,
        :allocs => median => :allocs)
end

function plot_rhs_metrics(df)
    rhs_df = aggregate_rhs(df)
    sort!(rhs_df, :median_time_ns)

    time_plot = bar(
        rhs_df.variant,
        rhs_df.median_time_ns,
        legend=false,
        xlabel="Variant",
        ylabel="Median RHS time (ns)",
        xtickfont=font(7),
        rotation=25,
        title="RHS performance",
    )

    alloc_plot = groupedbar(
        rhs_df.variant,
        [rhs_df.memory rhs_df.allocs],
        bar_position=:dodge,
        labels=["Bytes" "Allocations"],
        xlabel="Variant",
        ylabel="Per-call cost",
        xtickfont=font(7),
        rotation=25,
        title="RHS allocations",
    )

    return time_plot, alloc_plot
end

function rhs_bytes_plot(df)
    rhs_df = aggregate_rhs(df)
    sort!(rhs_df, :memory)
    plt = bar(
        rhs_df.variant,
        rhs_df.memory,
        legend=false,
        xlabel="Variant",
        ylabel="Bytes per RHS call",
        xtickfont=font(7),
        rotation=25,
        title="RHS memory footprint",
    )
    return plt
end

function trajectory_plot()
    u0 = [1.0, 0.0, 0.0]
    p = [0.2, 0.2, 5.7]
    tspan = (0.0, 100.0)
    prob = ODEProblem(rossler!, u0, tspan, p)
    sol = solve(prob, Tsit5(); saveat=0.1)
    plt = plot(sol, vars=(1,2), legend=false, xlabel="x", ylabel="y", title="$(TRAJECTORY_VARIANT) phase portrait")
    return plt
end

function write_solve_table(df; solver::AbstractString = DEFAULT_SOLVER, outfile::AbstractString)
    filtered = filter(row -> row.solver == solver && row.metric == "solve", df)
    isempty(filtered) && error("No solve rows available for solver '$solver'")
    sort!(filtered, :median_time_ns)
    open(outfile, "w") do io
        println(io, "\\begin{tabular}{lrrr}")
        println(io, "\\toprule")
        println(io, "Variant & Median time (ns) & Bytes & Allocs \\\")
        println(io, "\\midrule")
        for row in eachrow(filtered)
            variant = replace(row.variant, "_" => "\\_")
            println(io, "\\texttt{$variant} & $(round(row.median_time_ns; digits=2)) & $(row.memory) & $(row.allocs) \\\")
        end
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
    end
end

function main()
    df = load_benchmark_data(DATA_PATH)

    solve_plot = plot_solve_times(df)
    solve_out = joinpath(FIG_DIR, string(solver_slug(DEFAULT_SOLVER), "_solve_times.png"))
    savefig(solve_plot, solve_out)

    rhs_time_plot, rhs_alloc_plot = plot_rhs_metrics(df)
    savefig(rhs_time_plot, joinpath(FIG_DIR, "rhs_median_times.png"))
    savefig(rhs_alloc_plot, joinpath(FIG_DIR, "rhs_allocations.png"))

    savefig(rhs_bytes_plot(df), joinpath(FIG_DIR, "rhs_memory.png"))

    savefig(trajectory_plot(), joinpath(FIG_DIR, "trajectory_phase.png"))

    write_solve_table(df; solver=DEFAULT_SOLVER, outfile=joinpath(FIG_DIR, string(solver_slug(DEFAULT_SOLVER), "_table.tex")))

    println("Figures and tables written to $(FIG_DIR)")
end

main()
