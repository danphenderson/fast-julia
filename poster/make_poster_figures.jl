using DataFrames
using JSON3

include(joinpath(@__DIR__, "..", "rossler", "figures.jl"))

const LADDER_VARIANTS = [
    "rossler_naive",
    "rossler",
    "rossler_naive!",
    "rossler!",
    "rossler_static_naive",
    "rossler_static",
]

function load_dataframe(path::AbstractString)
    raw = JSON3.read(read(path, String))
    df = DataFrame(raw)
    df.metric = Symbol.(df.metric)
    return df
end

function ladder_variants(df::DataFrame)
    present = Set(String.(df.variant))
    ordered = String[]
    for v in LADDER_VARIANTS
        v in present && push!(ordered, v)
    end
    for v in sort(present)
        v in ordered && continue
        push!(ordered, v)
    end
    return ordered
end

function main(; data_path::AbstractString = joinpath(@__DIR__, "data.json"), outdir::AbstractString = joinpath(@__DIR__, "figures"))
    df = load_dataframe(data_path)
    PosterFigures.make_all_figures(df; outdir=outdir, baseline="rossler_naive", variants=ladder_variants(df))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
    println("Wrote figures into $(joinpath(@__DIR__, "figures"))")
end
