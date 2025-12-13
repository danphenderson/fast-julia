using JSON3
using DataFrames

const REQUIRED_COLUMNS = (
    :solver,
    :variant,
    :metric,
    :median_time_ns,
    :median_time_s,
    :memory,
    :allocs,
)

const REQUIRED_METRICS = Set(["rhs"])

"""
    load_poster_data(; path=joinpath(@__DIR__, "data.json"))

Read benchmarking results from `path` (defaults to `poster/data.json`) and return a
validated `DataFrame` suitable for plotting. The helper enforces several sanity
checks to keep poster figures and tables consistent:

* All required columns (`solver`, `variant`, `metric`, `median_time_ns`,
  `median_time_s`, `memory`, `allocs`) must be present and stored as
  vector-like fields of equal, nonzero length.
* Expected metrics (currently `"rhs"`, counting right-hand-side evaluations)
  must be present in the dataset.
* Timing columns are interpreted as median wall-clock runtimes (`median_time_s`
  in seconds; `median_time_ns` in nanoseconds) and must be numeric. Memory is
  assumed to be measured in bytes and `allocs` counts heap allocations; both
  must also be numeric.
* String columns (`solver`, `variant`, `metric`) are normalized to `String`.

The returned DataFrame contains `String` columns for identifiers and concrete
numeric types (`Int64` for counts and `Float64` for seconds) so downstream
plotting code can rely on consistent dtypes.
"""
function load_poster_data(; path=joinpath(@__DIR__, "data.json"))
    isfile(path) || error("data file not found at $(abspath(path))")

    raw = JSON3.read(read(path, String))

    missing = filter(col -> !(col in keys(raw)), REQUIRED_COLUMNS)
    isempty(missing) || error("missing required columns in $(path): $(collect(missing))")

    nonvectors = filter(col -> !(raw[col] isa AbstractVector), REQUIRED_COLUMNS)
    isempty(nonvectors) || error("columns are not vector-like: $(collect(nonvectors))")

    lengths = map(col -> length(raw[col]), REQUIRED_COLUMNS)
    isempty(lengths) && error("dataset in $(path) is empty")
    (length(unique(lengths)) == 1 && first(lengths) > 0) || error(
        "columns must have equal, nonzero length (observed lengths: $(lengths))",
    )

    metrics = raw[:metric]
    required_missing = setdiff(REQUIRED_METRICS, Set(metrics))
    isempty(required_missing) || error(
        "missing required metrics $(collect(required_missing)) in $(path)",
    )

    numeric_columns = [:median_time_ns, :median_time_s, :memory, :allocs]
    for col in numeric_columns
        all(v -> v isa Real, raw[col]) || error("column $(col) must be numeric")
    end

    df = DataFrame(
        solver = String.(raw[:solver]),
        variant = String.(raw[:variant]),
        metric = String.(metrics),
        median_time_ns = Int.(raw[:median_time_ns]),
        median_time_s = Float64.(raw[:median_time_s]),
        memory = Int.(raw[:memory]),
        allocs = Int.(raw[:allocs]),
    )

    return df
end

"""
    split_by_metric(df::DataFrame)

Convenience helper that groups the validated poster data by `:metric` and
returns a `Dict{String,DataFrame}` for quick plotting. Use this to access
metric-specific slices without re-validating the input.
"""
function split_by_metric(df::DataFrame)
    :metric in names(df) || error(":metric column missing from DataFrame")
    return Dict(String(key[1]) => copy(subdf) for (key, subdf) in pairs(groupby(df, :metric)))
end
