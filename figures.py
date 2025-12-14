"""figures.py

Poster visual pipeline for Fast Julia / Rössler benchmarks.

This script consumes the *flattened* benchmark results produced by the Julia
pipeline (see `BenchmarkFlatten.dataframe(...)` in `utils.jl`) and generates:

  - poster-ready PNGs (bar charts, speedups, Pareto)
  - optional LaTeX table snippets (booktabs)

Key design goals:
  - deterministic selection of `dt` for solve metrics (no accidental "first row")
  - filenames compatible with typical `poster.tex` includes

Expected CSV schema (minimum):
  solver, variant, metric, median_time_ns
Optional (recommended):
  dt, nsteps, memory, allocs

Usage examples:

  # Generate figures from an existing CSV
  python figures.py --csv poster/results.csv --outdir ./poster/figures --tables ./poster/tables

  # Choose which dt slice to use for solve bar charts
  python figures.py --csv poster/results.csv --dt-strategy max

Notes:
  - This file does *not* run Julia; it assumes you already exported results.
  - All charts use matplotlib defaults (no forced color palettes).
"""

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


DTStrategy = Literal["max", "min", "nearest"]


# -----------------------------------------------------------------------------
# Labels / ordering (mirrors the intent in `rossler/make_poster_figures.jl`)
# -----------------------------------------------------------------------------

def solver_slug(s: str) -> str:
    t = s.lower()
    t = "".join(ch if ch.isalnum() else "_" for ch in t)
    t = t.strip("_")
    while "__" in t:
        t = t.replace("__", "_")
    return t


def pretty_variant(v: str) -> str:
    mapping = {
        "rossler_naive": "naïve (alloc)",
        "rossler": "oop (tuned)",
        "rossler_naive!": "in-place (naïve)",
        "rossler!": "in-place (tuned)",
        "rossler_static_naive": "SVector (naïve)",
        "rossler_static": "SVector (tuned)",
        "rossler_type_stable": "type-stable",
        "rossler_ad": "AD-ready",
    }
    return mapping.get(v, v)


def preferred_variant_order(present: Iterable[str]) -> list[str]:
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
    present_set = set(map(str, present))
    ordered: list[str] = [v for v in preferred if v in present_set]
    # add any other variants deterministically
    for v in sorted(present_set):
        if v not in ordered:
            ordered.append(v)
    return ordered


# -----------------------------------------------------------------------------
# Schema normalization
# -----------------------------------------------------------------------------

def _ensure_col(df: pd.DataFrame, col: str) -> None:
    if col not in df.columns:
        raise ValueError(f"CSV missing required column '{col}'. Columns: {list(df.columns)}")


def normalize_schema(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for col in ("solver", "variant", "metric", "median_time_ns"):
        _ensure_col(df, col)

    # Julia side uses Symbol for metric; CSV may carry ":solve" or "solve".
    df["metric"] = df["metric"].astype(str).str.replace(":", "", regex=False)

    # Optional columns
    for opt in ("dt", "nsteps", "memory", "allocs"):
        if opt not in df.columns:
            df[opt] = np.nan

    df["median_time_ns"] = pd.to_numeric(df["median_time_ns"], errors="coerce")
    df["time_ms"] = df["median_time_ns"] * 1e-6
    df["time_s"] = df["median_time_ns"] * 1e-9

    # coerce optional numeric columns
    for c in ("dt", "nsteps", "memory", "allocs"):
        df[c] = pd.to_numeric(df[c], errors="coerce")

    return df


# -----------------------------------------------------------------------------
# dt selection for solve rows
# -----------------------------------------------------------------------------

def select_solve_slice(
    df: pd.DataFrame,
    *,
    solver: str,
    dt_strategy: DTStrategy,
    dt_target: float | None,
) -> pd.DataFrame:
    """Return a solve-only DataFrame with exactly one row per variant.

    The Julia experiment produces multiple dt values per (solver, variant).
    For bar charts and tables we must pick a single dt slice deterministically.
    """

    d = df[(df["solver"].astype(str) == solver) & (df["metric"] == "solve")].copy()
    if d.empty:
        return d

    if d["dt"].isna().all():
        # no dt column present; best-effort: pick first per variant
        return d.sort_values(["variant", "median_time_ns"]).groupby("variant", as_index=False).head(1)

    if dt_strategy == "nearest":
        if dt_target is None or not math.isfinite(dt_target):
            raise ValueError("--dt-strategy nearest requires --dt-target")

        def _pick(g: pd.DataFrame) -> pd.DataFrame:
            idx = (g["dt"] - dt_target).abs().idxmin()
            return g.loc[[idx]]

        return d.groupby("variant", group_keys=False).apply(_pick)

    if dt_strategy == "max":
        # dt0 (largest dt): max dt per variant
        idx = d.groupby("variant")["dt"].idxmax()
        return d.loc[idx].copy()

    if dt_strategy == "min":
        idx = d.groupby("variant")["dt"].idxmin()
        return d.loc[idx].copy()

    raise ValueError(f"Unknown dt_strategy: {dt_strategy}")


# -----------------------------------------------------------------------------
# Plotting helpers
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class PlotSpec:
    dpi: int = 300
    figsize: tuple[float, float] = (11.0, 6.5)


def _hbar(
    labels: list[str],
    values: list[float],
    *,
    xlabel: str,
    title: str,
    xscale: str = "linear",
    spec: PlotSpec = PlotSpec(),
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=spec.figsize, dpi=spec.dpi)
    y = np.arange(len(labels))
    ax.barh(y, values)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.set_xscale(xscale)
    ax.grid(True, axis="x", which="both", linestyle=":", linewidth=0.7)
    fig.tight_layout()
    return fig


def plot_solver_times_ms(
    df_slice: pd.DataFrame,
    *,
    solver: str,
    variants: list[str],
    metric_label: str,
    title_suffix: str = "",
    xscale: str = "linear",
    spec: PlotSpec = PlotSpec(),
) -> plt.Figure:
    labels: list[str] = []
    times: list[float] = []
    for v in variants:
        rows = df_slice[df_slice["variant"].astype(str) == v]
        if rows.empty:
            continue
        labels.append(pretty_variant(v))
        times.append(float(rows.iloc[0]["time_ms"]))
    # match poster ordering: top is most interesting (reverse)
    labels = list(reversed(labels))
    times = list(reversed(times))
    return _hbar(
        labels,
        times,
        xlabel="median time (ms)",
        title=f"{solver}: {metric_label} times {title_suffix}".strip(),
        xscale=xscale,
        spec=spec,
    )


def plot_solver_speedup(
    df_slice: pd.DataFrame,
    *,
    solver: str,
    variants: list[str],
    baseline: str,
    metric_label: str,
    xscale: str = "log",
    spec: PlotSpec = PlotSpec(),
) -> plt.Figure:
    base_rows = df_slice[df_slice["variant"].astype(str) == baseline]
    if base_rows.empty and baseline != "rossler_naive":
        base_rows = df_slice[df_slice["variant"].astype(str) == "rossler_naive"]
        baseline = "rossler_naive"
    if base_rows.empty:
        raise ValueError(f"Baseline '{baseline}' missing for solver '{solver}'.")

    base = float(base_rows.iloc[0]["median_time_ns"])
    labels: list[str] = []
    speed: list[float] = []
    for v in variants:
        rows = df_slice[df_slice["variant"].astype(str) == v]
        if rows.empty:
            continue
        labels.append(pretty_variant(v))
        speed.append(base / float(rows.iloc[0]["median_time_ns"]))
    labels = list(reversed(labels))
    speed = list(reversed(speed))
    return _hbar(
        labels,
        speed,
        xlabel=f"speedup vs baseline ({baseline}) (×)",
        title=f"{solver}: speedup ({metric_label})",
        xscale=xscale,
        spec=spec,
    )


def plot_solver_allocs(
    df_slice: pd.DataFrame,
    *,
    solver: str,
    variants: list[str],
    metric_label: str,
    xscale: str = "linear",
    spec: PlotSpec = PlotSpec(),
) -> plt.Figure:
    if "allocs" not in df_slice.columns:
        raise ValueError("DataFrame has no 'allocs' column")
    labels: list[str] = []
    allocs: list[float] = []
    for v in variants:
        rows = df_slice[df_slice["variant"].astype(str) == v]
        if rows.empty:
            continue
        labels.append(pretty_variant(v))
        allocs.append(float(rows.iloc[0]["allocs"]))
    labels = list(reversed(labels))
    allocs = list(reversed(allocs))
    return _hbar(
        labels,
        allocs,
        xlabel="allocations per call",
        title=f"{solver}: allocations ({metric_label})",
        xscale=xscale,
        spec=spec,
    )


def plot_pareto(
    df: pd.DataFrame,
    *,
    metric: str,
    spec: PlotSpec = PlotSpec(),
) -> plt.Figure:
    d = df[df["metric"] == metric].copy()
    d = d.dropna(subset=["allocs", "time_ms"])
    fig, ax = plt.subplots(figsize=spec.figsize, dpi=spec.dpi)
    for solver, g in d.groupby("solver"):
        ax.scatter(g["allocs"], g["time_ms"], label=str(solver), s=40)
    ax.set_xlabel("allocations per call")
    ax.set_ylabel("median time (ms)")
    ax.set_title(f"Pareto view: time vs allocations (metric = {metric})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="both", linestyle=":", linewidth=0.7)
    ax.legend(loc="best")
    fig.tight_layout()
    return fig


def plot_dt_sweep(
    df: pd.DataFrame,
    *,
    solver: str,
    variants: list[str],
    x: Literal["dt", "nsteps"] = "nsteps",
    y: Literal["time_ms", "median_time_ns"] = "time_ms",
    loglog: bool = True,
    spec: PlotSpec = PlotSpec(),
) -> plt.Figure:
    d = df[(df["solver"].astype(str) == solver) & (df["metric"] == "solve")].copy()
    d = d.dropna(subset=["dt", y])
    if x == "nsteps" and d["nsteps"].isna().all():
        x = "dt"  # fallback

    fig, ax = plt.subplots(figsize=spec.figsize, dpi=spec.dpi)
    for v in variants:
        g = d[d["variant"].astype(str) == v].sort_values("dt", ascending=False)
        if g.empty:
            continue
        ax.plot(g[x], g[y], marker="o", label=pretty_variant(v))
    ax.set_xlabel(x)
    ax.set_ylabel("median time (ms)" if y == "time_ms" else y)
    ax.set_title(f"{solver}: solve cost vs {x}")
    ax.grid(True, which="both", linestyle=":", linewidth=0.7)
    if loglog:
        ax.set_xscale("log")
        ax.set_yscale("log")
    ax.legend(loc="best")
    fig.tight_layout()
    return fig


# -----------------------------------------------------------------------------
# LaTeX tables
# -----------------------------------------------------------------------------


def latex_escape(s: str) -> str:
    return (
        s.replace("\\", r"\textbackslash{}")
        .replace("_", r"\_")
        .replace("%", r"\%")
        .replace("&", r"\&")
        .replace("#", r"\#")
    )


def make_latex_table(
    df_slice: pd.DataFrame,
    *,
    solver: str,
    metric: str,
    variants: list[str],
    caption: str | None = None,
    label: str | None = None,
) -> str:
    """Return a booktabs LaTeX table snippet (tabular only)."""

    cols = ["Variant", "Median (ns)", "Median (ms)"]
    has_allocs = "allocs" in df_slice.columns and not df_slice["allocs"].isna().all()
    has_mem = "memory" in df_slice.columns and not df_slice["memory"].isna().all()
    has_dt = metric == "solve" and "dt" in df_slice.columns and not df_slice["dt"].isna().all()
    has_steps = metric == "solve" and "nsteps" in df_slice.columns and not df_slice["nsteps"].isna().all()

    if has_dt:
        cols.insert(1, "dt")
    if has_steps:
        cols.insert(2 if has_dt else 1, "nsteps")
    if has_allocs:
        cols += ["Allocs", "Bytes"]
    elif has_mem:
        cols += ["Bytes"]

    lines: list[str] = []
    lines.append(r"\begin{tabular}{l" + "r" * (len(cols) - 1) + r"}")
    lines.append(r"\toprule")
    lines.append(" ".join([c + (" &" if i < len(cols) - 1 else r" \\") for i, c in enumerate(cols)]))
    lines.append(r"\midrule")

    for v in variants:
        rows = df_slice[df_slice["variant"].astype(str) == v]
        if rows.empty:
            continue
        r0 = rows.iloc[0]
        parts: list[str] = [latex_escape(pretty_variant(v))]
        if has_dt:
            parts.append(f"{float(r0['dt']):.3g}")
        if has_steps:
            ns = r0.get("nsteps", np.nan)
            parts.append("--" if pd.isna(ns) else f"{int(ns)}")
        parts.append(f"{float(r0['median_time_ns']):.3g}")
        parts.append(f"{float(r0['time_ms']):.3g}")
        if has_allocs:
            a = r0.get("allocs", np.nan)
            m = r0.get("memory", np.nan)
            parts.append("--" if pd.isna(a) else f"{int(a)}")
            parts.append("--" if pd.isna(m) else f"{int(m)}")
        elif has_mem:
            m = r0.get("memory", np.nan)
            parts.append("--" if pd.isna(m) else f"{int(m)}")
        lines.append(" & ".join(parts) + r" \\")

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")

    # If caller wants a complete table environment, they can wrap it; here we keep it snippet-friendly.
    if caption or label:
        # Provide commented helpers for convenience.
        if caption:
            lines.insert(0, f"% caption: {caption}")
        if label:
            lines.insert(0, f"% label: {label}")

    return "\n".join(lines) + "\n"


# -----------------------------------------------------------------------------
# Main pipeline
# -----------------------------------------------------------------------------


def make_all_figures(
    df: pd.DataFrame,
    *,
    outdir: Path,
    tables_dir: Path | None,
    baseline: str,
    dt_strategy: DTStrategy,
    dt_target: float | None,
    poster_compat: bool,
    include_dt_sweep: bool,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    if tables_dir is not None:
        tables_dir.mkdir(parents=True, exist_ok=True)

    solvers = sorted(df["solver"].astype(str).unique())
    for solver in solvers:
        tag = solver_slug(solver)
        present_variants = df[df["solver"].astype(str) == solver]["variant"].astype(str).unique()
        variants = preferred_variant_order(present_variants)

        # RHS slice (one row per variant already)
        rhs = df[(df["solver"].astype(str) == solver) & (df["metric"] == "rhs")].copy()

        # Solve slice: pick a dt per variant deterministically.
        solve = select_solve_slice(df, solver=solver, dt_strategy=dt_strategy, dt_target=dt_target)

        # Build and save figures
        if not solve.empty:
            fig = plot_solver_speedup(solve, solver=solver, variants=variants, baseline=baseline, metric_label="solve", xscale="log")
            fig.savefig(outdir / f"{tag}_speedup_solve_log10.png")
            plt.close(fig)

            fig = plot_solver_times_ms(solve, solver=solver, variants=variants, metric_label="solve")
            fig.savefig(outdir / f"{tag}_solve_times_ms.png")
            plt.close(fig)

            if "allocs" in solve.columns and not solve["allocs"].isna().all():
                fig = plot_solver_allocs(solve, solver=solver, variants=variants, metric_label="solve")
                fig.savefig(outdir / f"{tag}_allocs_solve.png")
                plt.close(fig)

        if not rhs.empty:
            fig = plot_solver_times_ms(rhs, solver=solver, variants=variants, metric_label="rhs")
            fig.savefig(outdir / f"{tag}_rhs_times_ms.png")
            plt.close(fig)

            if "allocs" in rhs.columns and not rhs["allocs"].isna().all():
                fig = plot_solver_allocs(rhs, solver=solver, variants=variants, metric_label="rhs")
                fig.savefig(outdir / f"{tag}_allocs_rhs.png")
                plt.close(fig)

        if include_dt_sweep:
            dsolve = df[(df["solver"].astype(str) == solver) & (df["metric"] == "solve")].copy()
            if not dsolve.empty and not dsolve["dt"].isna().all():
                fig = plot_dt_sweep(df, solver=solver, variants=variants, x="nsteps", y="time_ms", loglog=True)
                fig.savefig(outdir / f"{tag}_solve_dt_sweep_loglog.png")
                plt.close(fig)

        # Tables
        if tables_dir is not None:
            if not rhs.empty:
                tex = make_latex_table(rhs, solver=solver, metric="rhs", variants=variants)
                (tables_dir / f"{tag}_rhs_table.tex").write_text(tex, encoding="utf-8")
            if not solve.empty:
                tex = make_latex_table(solve, solver=solver, metric="solve", variants=variants)
                (tables_dir / f"{tag}_solve_table.tex").write_text(tex, encoding="utf-8")

        # Poster-compat filenames (keeps your LaTeX stable)
        if poster_compat and ("rk4" in solver.lower()):
            # your Julia-side convention: rk4_time_normalized_log10.png
            # here we reuse the speedup plot.
            src = outdir / f"{tag}_speedup_solve_log10.png"
            if src.exists():
                (outdir / "rk4_time_normalized_log10.png").write_bytes(src.read_bytes())

            src = outdir / f"{tag}_allocs_solve.png"
            if src.exists():
                (outdir / "rk4_allocs_solve.png").write_bytes(src.read_bytes())

    # Global Pareto
    if "allocs" in df.columns and not df["allocs"].isna().all():
        fig = plot_pareto(df, metric="solve")
        fig.savefig(outdir / "pareto_solve.png")
        plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate poster figures from flattened benchmark CSV")
    p.add_argument("--csv", required=True, help="Path to flattened CSV (from Julia BenchmarkFlatten.dataframe)")
    p.add_argument("--outdir", default="./poster/figures", help="Output directory for PNGs")
    p.add_argument("--tables", default="", help="If set, write LaTeX table snippets into this directory")
    p.add_argument("--baseline", default="rossler_naive", help="Baseline variant for speedup charts")
    p.add_argument(
        "--dt-strategy",
        choices=["max", "min", "nearest"],
        default="max",
        help="How to choose dt for solve bar charts when multiple dt values exist",
    )
    p.add_argument(
        "--dt-target",
        type=float,
        default=None,
        help="Target dt used when --dt-strategy nearest",
    )
    p.add_argument(
        "--poster-compat",
        action="store_true",
        help="Also write poster-compat filenames (e.g., rk4_time_normalized_log10.png)",
    )
    p.add_argument(
        "--dt-sweep",
        action="store_true",
        help="Also generate a dt-sweep (log-log) plot for each solver",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    csv_path = Path(args.csv)
    if not csv_path.exists():
        raise FileNotFoundError(csv_path)

    df = pd.read_csv(csv_path)
    df = normalize_schema(df)

    outdir = Path(args.outdir)
    tables_dir = Path(args.tables) if args.tables.strip() else None

    make_all_figures(
        df,
        outdir=outdir,
        tables_dir=tables_dir,
        baseline=str(args.baseline),
        dt_strategy=str(args.dt_strategy),
        dt_target=args.dt_target,
        poster_compat=bool(args.poster_compat),
        include_dt_sweep=bool(args.dt_sweep),
    )

    print(f"Wrote figures into: {outdir}")
    if tables_dir is not None:
        print(f"Wrote LaTeX tables into: {tables_dir}")


if __name__ == "__main__":
    main()
