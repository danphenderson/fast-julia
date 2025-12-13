Generated poster assets live here. The figure pipeline should write:
- midpoint_solve_times.png
- rhs_allocations.png
- trajectory.png
- midpoint_table.tex

Run `julia --project -e 'include("poster/make_poster_figures.jl")'` to regenerate the files before compiling `poster/poster.tex`.
