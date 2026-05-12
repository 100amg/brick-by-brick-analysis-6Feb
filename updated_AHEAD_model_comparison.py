#!/usr/bin/env python3
"""
Per-Brick Model Comparison for AHEAD Methylation ASCII Logs
============================================================
Parses ASCII log files from the modkit/B2A pipeline and compares per-brick
accuracy across methylation calling models for all AHEAD samples.

Log files must follow the naming convention:
    <anything>_<model>.log
    e.g.  ASCII_Log_AHEAD_10_11_25_3_0_sorted_bam_160226_212638_dorado_fast.log

Recognised model suffixes (edit KNOWN_MODELS to extend):
    dorado_fast, dorado_hac, dorado_sup, guppy_hac, deepmod2

Outputs:
    <prefix>_per_brick.csv
    <prefix>_summary.csv
    <prefix>_figures.pdf          (all 4 plots in a single multi-page PDF)

Usage:
    python per_brick_model_comparison.py <log_directory> [output_prefix]
"""

import re
import sys
import csv
from pathlib import Path
from collections import defaultdict

from gate_analysis import run_gate_analysis

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

# ── Global plot style ─────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         10,
    "axes.titlesize":    12,
    "axes.titleweight":  "bold",
    "axes.labelsize":    10,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.grid":         True,
    "grid.linestyle":    "--",
    "grid.alpha":        0.45,
    "axes.axisbelow":    True,
    "legend.framealpha": 0.85,
    "figure.dpi":        150,
})

# ── CpG count per brick position ─────────────────────────────────────────────
# Derived from the construct design table (brick number → CG position → n CpGs).
# Key = CG site position (bp), Value = number of CpGs in that brick.
# Edit or extend this dict if your construct changes.
CPG_COUNT_BY_SITE = {
     60: 2,   84: 2,  108: 3,  132: 2,  156: 1,  180: 2,
    204: 2,  228: 1,  252: 1,  276: 2,  300: 1,  324: 2,
    348: 1,  372: 1,  396: 3,  420: 2,  444: 1,  468: 2,
    492: 1,  516: 2,  540: 2,  564: 2,  588: 3,  612: 3,
    636: 1,  660: 1,  684: 2,  708: 2,  732: 2,  756: 2,
    780: 2,  804: 2,  828: 1,  852: 2,  876: 2,  900: 1,
}

# ── Model registry ────────────────────────────────────────────────────────────
KNOWN_MODELS = [
    "dorado_fast",
    "dorado_hac",
    "dorado_sup",
    "guppy_hac",
    "deepmod2",
    "minknow",
]

# Consistent colour per model across every plot
MODEL_COLOURS = {
    "dorado_fast": "#4C72B0",
    "dorado_hac":  "#DD8452",
    "dorado_sup":  "#55A868",
    "guppy_hac":   "#C44E52",
    "deepmod2":    "#8172B2",
    "minknow" : '#40E0D0',
}
DEFAULT_COLOUR = "#888888"

# Error-type colours (shared across plots 3 & 4)
COLOUR_CORRECT = "#4CAF50"   # green
COLOUR_FN      = "#FFC107"   # amber  — 1→0 missed methylation
COLOUR_FP      = "#F44336"   # red    — 0→1 false positive


def model_colour(m):
    return MODEL_COLOURS.get(m, DEFAULT_COLOUR)


# ── Filename helpers ───────────────────────────────────────────────────────────

def detect_model(filename):
    stem = Path(filename).stem.lower()
    for model in KNOWN_MODELS:
        if stem.endswith("_" + model):
            return model
    return None


def sample_name_from_file(filename, model):
    stem = Path(filename).stem
    return re.sub(re.escape("_" + model), "", stem, flags=re.IGNORECASE)


# ── Parsing helpers ────────────────────────────────────────────────────────────

def parse_int_list(line):
    m = re.search(r"\[([^\]]*)\]", line)
    if not m or not m.group(1).strip():
        return []
    return [int(x.strip()) for x in m.group(1).split(",") if x.strip()]


def parse_float_list(line):
    m = re.search(r"\[([^\]]*)\]", line)
    if not m or not m.group(1).strip():
        return []
    return [float(x.strip()) for x in m.group(1).split(",") if x.strip()]


def parse_log(filepath):
    """
    Parse one ASCII log file.
    Returns a dict or None if the essential CG-sites/status fields are absent.
    """
    with open(filepath, "r", errors="replace") as fh:
        lines = fh.readlines()

    r = dict(
        cg_sites=[], actual_status=[], desired_sites=[],
        flip_1to0_bricks=[], flip_0to1_bricks=[],
        total_positions=None, matches=None,
        flips_1to0=None, flips_0to1=None, error_pct=None,
        ascii_chars=[], decoded_word="",
    )
    in_1to0 = False
    in_0to1 = False

    for line in lines:
        s = line.strip()

        if s.startswith("CG sites:"):
            r["cg_sites"] = parse_int_list(s)
        elif s.startswith("Methylation status:"):
            r["actual_status"] = [int(v) for v in parse_float_list(s)]
        elif s.startswith("ASCII characters:"):
            m = re.search(r"\[([^\]]*)\]", s)
            if m:
                chars = re.findall(r"'([^']*)'", m.group(1))
                r["ascii_chars"] = chars
                r["decoded_word"] = "".join(chars)
        elif s.startswith("Desired methylated CG sites are:"):
            r["desired_sites"] = parse_int_list(s)
        elif s.startswith("Total positions compared:"):
            m = re.search(r"(\d+)", s)
            if m:
                r["total_positions"] = int(m.group(1))
        elif s.startswith("Matches (desired == actual):"):
            m = re.search(r"(\d+)", s)
            if m:
                r["matches"] = int(m.group(1))
        elif re.match(r"1 -> 0 flips", s):
            m = re.search(r":\s*(\d+)", s)
            if m:
                r["flips_1to0"] = int(m.group(1))
            in_1to0 = True; in_0to1 = False
        elif re.match(r"0 -> 1 flips", s):
            m = re.search(r":\s*(\d+)", s)
            if m:
                r["flips_0to1"] = int(m.group(1))
            in_0to1 = True; in_1to0 = False
        elif s.startswith("Error Percentage:"):
            m = re.search(r"([\d.]+)", s)
            if m:
                r["error_pct"] = float(m.group(1))
            in_1to0 = False; in_0to1 = False
        elif "Bit/Brick number=" in s:
            m = re.search(r"Bit/Brick number=\s*(\d+)", s)
            if m:
                idx = int(m.group(1)) - 1   # convert to 0-based
                if in_1to0:
                    r["flip_1to0_bricks"].append(idx)
                elif in_0to1:
                    r["flip_0to1_bricks"].append(idx)

    if not r["cg_sites"] or not r["actual_status"]:
        return None
    return r


def desired_vector(cg_sites, desired_sites):
    ds = set(desired_sites)
    return [1 if s in ds else 0 for s in cg_sites]


def compute_stats(actual, desired):
    total   = len(actual)
    correct = sum(a == d for a, d in zip(actual, desired))
    tp = sum(a == 1 and d == 1 for a, d in zip(actual, desired))
    fn = sum(a == 0 and d == 1 for a, d in zip(actual, desired))
    tn = sum(a == 0 and d == 0 for a, d in zip(actual, desired))
    fp = sum(a == 1 and d == 0 for a, d in zip(actual, desired))
    acc  = 100 * correct / total if total else 0
    sens = 100 * tp / (tp + fn) if (tp + fn) else None
    spec = 100 * tn / (tn + fp) if (tn + fp) else None
    # Matthews Correlation Coefficient — honest single metric for imbalanced classes.
    # Ranges from -1 (perfectly wrong) through 0 (chance) to +1 (perfect).
    # Unlike accuracy or F1, it uses all four confusion matrix cells so it
    # cannot be inflated by predicting the majority class (desired=0 here).
    denom = ((tp+fp) * (tp+fn) * (tn+fp) * (tn+fn)) ** 0.5
    mcc   = (tp*tn - fp*fn) / denom if denom > 0 else 0.0
    return dict(total=total, correct=correct,
                tp=tp, fn=fn, tn=tn, fp=fp,
                acc=acc, sens=sens, spec=spec, mcc=mcc)


# ── Plot 1 — Overall accuracy / sensitivity / specificity ─────────────────────

def plot_overall(ax, model_stats, all_models):
    metrics = [
        ("acc",  "Accuracy"),
        ("sens", "Sensitivity (1→1)"),
        ("spec", "Specificity (0→0)"),
        ("mcc",  "MCC ×100"),          # scaled to 0–100 for same axis
    ]
    alphas  = [1.0, 0.70, 0.40, 0.60]
    bar_w   = 0.18
    x       = np.arange(len(all_models))

    for i, ((key, label), alpha) in enumerate(zip(metrics, alphas)):
        vals = []
        for m in all_models:
            v = model_stats[m].get(key)
            # MCC is −1 to +1; scale to same axis as percentages
            if key == "mcc":
                vals.append((v or 0) * 100)
            else:
                vals.append(v if v is not None else 0)
        offset = (i - 1.5) * bar_w
        bars = ax.bar(x + offset, vals, bar_w, label=label,
                      color=[model_colour(m) for m in all_models],
                      alpha=alpha, edgecolor="white", linewidth=0.6)
        for bar, val in zip(bars, vals):
            if abs(val) > 2:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.8,
                        f"{val:.1f}", ha="center", va="bottom",
                        fontsize=6.5, rotation=40)

    ax.set_xticks(x)
    ax.set_xticklabels([m.replace("_", "\n") for m in all_models], fontsize=9)
    ax.set_ylabel("Percentage / Score (%)")
    ax.set_ylim(-20, 130)
    ax.axhline(0, color="black", linewidth=0.6, linestyle="--", alpha=0.4)
    ax.set_title("Plot 1 — Overall Model Performance\n"
                 "(MCC ×100 shown on same axis; MCC=0 → chance, MCC=100 → perfect)")

    metric_patches = [mpatches.Patch(facecolor="grey", alpha=a, label=l)
                      for a, (_, l) in zip(alphas, metrics)]
    ax.legend(handles=metric_patches, title="Metric", fontsize=8, loc="upper right")


# ── Plot 2 — Per-sample accuracy ──────────────────────────────────────────────

def plot_per_sample(ax, data, all_samples, all_models):
    n_models = len(all_models)
    bar_w    = 0.75 / n_models
    x        = np.arange(len(all_samples))

    for i, model in enumerate(all_models):
        accs = []
        for sample in all_samples:
            if model not in data[sample]:
                accs.append(0)
                continue
            p   = data[sample][model]
            des = desired_vector(p["cg_sites"], p["desired_sites"])
            st  = compute_stats(p["actual_status"], des)
            accs.append(st["acc"])
        offset = (i - (n_models - 1) / 2) * bar_w
        bars = ax.bar(x + offset, accs, bar_w,
                      label=model.replace("_", " "),
                      color=model_colour(model),
                      edgecolor="white", linewidth=0.5)
        for bar, val in zip(bars, accs):
            if val > 2:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.5,
                        f"{val:.0f}%", ha="center", va="bottom",
                        fontsize=6.5, rotation=40)

    short_samples = []
    for s in all_samples:
        m = re.search(r"(AHEAD[_.\d]+)", s, re.IGNORECASE)
        short_samples.append(m.group(1) if m else s[-35:])

    ax.set_xticks(x)
    ax.set_xticklabels(short_samples, rotation=18, ha="right", fontsize=8)
    ax.set_ylabel("Accuracy (%)")
    ax.set_ylim(0, 120)
    ax.set_title("Plot 2 — Per-Sample Accuracy by Model")
    ax.legend(title="Model", fontsize=8,
              bbox_to_anchor=(1.01, 1), loc="upper left")


# ── Plot 3 — Flip type breakdown ──────────────────────────────────────────────

def plot_flip_breakdown(ax, data, all_samples, all_models):
    totals = {m: dict(correct=0, fn=0, fp=0) for m in all_models}

    for sample in all_samples:
        for model in all_models:
            if model not in data[sample]:
                continue
            p   = data[sample][model]
            des = desired_vector(p["cg_sites"], p["desired_sites"])
            st  = compute_stats(p["actual_status"], des)
            totals[model]["correct"] += st["correct"]
            totals[model]["fn"]      += st["fn"]
            totals[model]["fp"]      += st["fp"]

    x     = np.arange(len(all_models))
    bar_w = 0.5

    corrects = [totals[m]["correct"] for m in all_models]
    fns      = [totals[m]["fn"]      for m in all_models]
    fps      = [totals[m]["fp"]      for m in all_models]

    ax.bar(x, corrects, bar_w, label="Correct",
           color=[model_colour(m) for m in all_models],
           alpha=0.88, edgecolor="white", linewidth=0.5)
    ax.bar(x, fns, bar_w, bottom=corrects,
           label="1→0  missed methylation",
           color=COLOUR_FN, alpha=0.9, edgecolor="white", linewidth=0.5)
    ax.bar(x, fps, bar_w,
           bottom=[c + f for c, f in zip(corrects, fns)],
           label="0→1  false positive",
           color=COLOUR_FP, alpha=0.9, edgecolor="white", linewidth=0.5)

    for i, (m, c, fn_v, fp_v) in enumerate(zip(all_models, corrects, fns, fps)):
        total = c + fn_v + fp_v
        if total == 0:
            continue
        ax.text(i, c / 2,
                f"{100*c/total:.0f}%",
                ha="center", va="center",
                fontsize=8.5, fontweight="bold", color="white")
        if fn_v > 0:
            ax.text(i, c + fn_v / 2,
                    f"{100*fn_v/total:.0f}%",
                    ha="center", va="center", fontsize=8)
        if fp_v > 0:
            ax.text(i, c + fn_v + fp_v / 2,
                    f"{100*fp_v/total:.0f}%",
                    ha="center", va="center", fontsize=8, color="white")

    ax.set_xticks(x)
    ax.set_xticklabels([m.replace("_", "\n") for m in all_models], fontsize=9)
    ax.set_ylabel("Number of bricks (all samples)")
    ax.set_title("Plot 3 — Error Type Breakdown per Model (all samples combined)")
    ax.legend(fontsize=8, loc="upper right")


# ── Plot 5 — MCC scatter plot ─────────────────────────────────────────────────

def plot_mcc_scatter(ax, data, all_samples, all_models):
    """
    Scatter plot: x = model (jittered), y = MCC score per sample.
    Each sample is a dot; a horizontal bar shows the mean across samples.
    A dashed line at MCC=0 marks the 'no better than chance' baseline.

    Why MCC and not accuracy:
      MCC uses all four confusion-matrix cells (TP, TN, FP, FN) and is not
      inflated when one class (desired=0) greatly outnumbers the other.
      A model that calls everything 0 would score high accuracy but MCC ≈ 0.
    """
    np.random.seed(42)   # reproducible jitter
    x_positions = {m: i for i, m in enumerate(all_models)}

    for mi, model in enumerate(all_models):
        mccs = []
        for sample in all_samples:
            if model not in data[sample]:
                continue
            p   = data[sample][model]
            des = desired_vector(p["cg_sites"], p["desired_sites"])
            st  = compute_stats(p["actual_status"], des)
            mccs.append(st["mcc"])

        if not mccs:
            continue

        jitter = np.random.uniform(-0.12, 0.12, size=len(mccs))
        xs     = [mi + j for j in jitter]

        # Individual sample points
        ax.scatter(xs, mccs,
                   color=model_colour(model), s=60, zorder=3,
                   alpha=0.85, edgecolors="white", linewidths=0.6)

        # Label each dot with sample short-name
        for x, mcc_val, sample in zip(xs, mccs, all_samples):
            m = re.search(r"AHEAD[_.](\d+[_.]\d+[_.]\d+)", sample, re.IGNORECASE)
            lbl = m.group(1).replace("_", "/") if m else ""
            ax.annotate(lbl, (x, mcc_val),
                        textcoords="offset points", xytext=(0, 5),
                        ha="center", fontsize=6, color="grey")

        # Mean bar
        mean_mcc = np.mean(mccs)
        ax.plot([mi - 0.25, mi + 0.25], [mean_mcc, mean_mcc],
                color=model_colour(model), linewidth=2.5, zorder=4,
                solid_capstyle="round")
        ax.annotate(f"{mean_mcc:.3f}",
                    (mi, mean_mcc),
                    textcoords="offset points", xytext=(0, 7),
                    va="center", fontsize=8, fontweight="bold",
                    color=model_colour(model))

    # Reference line at MCC = 0
    ax.axhline(0, color="black", linewidth=1.0, linestyle="--", alpha=0.5,
               label="MCC = 0  (chance)")
    ax.axhline(1, color="green", linewidth=0.6, linestyle=":", alpha=0.4,
               label="MCC = 1  (perfect)")

    ax.set_xticks(range(len(all_models)))
    ax.set_xticklabels([m.replace("_", "\n") for m in all_models], fontsize=9)
    ax.set_ylabel("MCC  (Matthews Correlation Coefficient)")
    ax.set_ylim(-1.05, 1.15)
    ax.set_xlim(-0.6, len(all_models) - 0.4)
    ax.set_title(
        "Plot 5 — MCC per Model per Sample\n"
        "Dots = individual replicates  |  Bar = mean across replicates\n"
        "MCC = 0 → no better than chance  |  MCC = 1 → perfect",
        fontsize=10,
    )
    ax.legend(fontsize=8, loc="lower right")
    ax.yaxis.grid(True, linestyle="--", alpha=0.4)
    ax.set_axisbelow(True)


# ── Plot 4 — Per-brick heatmap with CpG count track ──────────────────────────

def plot_brick_heatmap(fig, all_samples, all_models, data):
    """
    For each sample: a heatmap (models × bricks) plus a CpG-count bar track
    directly beneath it, sharing the same x-axis so columns align perfectly.

    Heatmap colours:
      Green  = correct call
      Amber  = 1→0 missed methylation
      Red    = 0→1 false positive
      White  = model absent for this sample

    CpG bar track:
      Steel-blue bars showing the number of CpGs in each brick.
      The number is printed inside each bar.
    """
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.gridspec import GridSpecFromSubplotSpec

    n_samples  = len(all_samples)
    cmap       = ListedColormap([COLOUR_FP, COLOUR_FN, COLOUR_CORRECT])
    norm       = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    cpg_colour = "#5B8DB8"   # steel blue for CpG bars

    # Outer grid: one row per sample, each split into heatmap + bar track
    # height_ratios keeps the heatmap ~4x taller than the bar track
    outer_gs = fig.add_gridspec(
        n_samples, 1,
        hspace=0.70,
        top=0.90, bottom=0.08,
        left=0.11, right=0.97,
    )

    for row_idx, sample in enumerate(all_samples):

        # Inner grid for this sample: row 0 = heatmap, row 1 = CpG bars
        inner_gs = GridSpecFromSubplotSpec(
            2, 1,
            subplot_spec=outer_gs[row_idx],
            height_ratios=[4, 1],
            hspace=0.08,
        )
        ax_heat = fig.add_subplot(inner_gs[0])
        ax_cpg  = fig.add_subplot(inner_gs[1])

        # ── Collect data for this sample ──────────────────────────────────────
        n_bricks = max(
            (len(data[sample][m]["cg_sites"])
             for m in all_models if m in data[sample]),
            default=0,
        )
        if n_bricks == 0:
            ax_heat.set_visible(False)
            ax_cpg.set_visible(False)
            continue

        matrix   = np.full((len(all_models), n_bricks), np.nan)
        cg_sites = None

        for mi, model in enumerate(all_models):
            if model not in data[sample]:
                continue
            p        = data[sample][model]
            cg_sites = p["cg_sites"]
            des      = desired_vector(p["cg_sites"], p["desired_sites"])
            act      = p["actual_status"]
            for bi, (a, d) in enumerate(zip(act, des)):
                if a == d:
                    matrix[mi, bi] = 2
                elif d == 0 and a == 1:
                    matrix[mi, bi] = 0
                else:
                    matrix[mi, bi] = 1

        # ── Heatmap ───────────────────────────────────────────────────────────
        masked = np.ma.masked_invalid(matrix)
        ax_heat.imshow(masked, cmap=cmap, norm=norm,
                       aspect="auto", interpolation="nearest")

        ax_heat.set_yticks(range(len(all_models)))
        ax_heat.set_yticklabels([m.replace("_", "\n") for m in all_models],
                                fontsize=7)
        ax_heat.set_xticks([])          # x-axis shared with CpG track below
        ax_heat.tick_params(bottom=False)

        m_title = re.search(r"(AHEAD[_.\d]+)", sample, re.IGNORECASE)
        short   = m_title.group(1) if m_title else sample[-45:]
        ax_heat.set_title(f"Sample: {short}", fontsize=9, pad=3)

        # ── CpG count bar track ───────────────────────────────────────────────
        if cg_sites:
            cpg_counts = [CPG_COUNT_BY_SITE.get(site, 0) for site in cg_sites]
            bar_x      = np.arange(n_bricks)

            bars = ax_cpg.bar(bar_x, cpg_counts, width=0.85,
                              color=cpg_colour, alpha=0.80,
                              linewidth=0)

            # Print the count number centred inside / above each bar
            max_cpg = max(cpg_counts) if cpg_counts else 1
            for bx, val in zip(bar_x, cpg_counts):
                if val == 0:
                    continue
                # Place text at 50% bar height; if bar is too short, place above
                y_pos  = val / 2 if val >= 1 else val + 0.05
                colour = "white" if val >= 2 else "black"
                ax_cpg.text(bx, y_pos, str(val),
                            ha="center", va="center",
                            fontsize=5.5, fontweight="bold", color=colour)

            ax_cpg.set_xlim(ax_heat.get_xlim())
            ax_cpg.set_ylim(0, max_cpg + 0.6)
            ax_cpg.set_yticks([1, max_cpg] if max_cpg > 1 else [1])
            ax_cpg.set_yticklabels(
                [str(v) for v in ([1, max_cpg] if max_cpg > 1 else [1])],
                fontsize=6,
            )
            ax_cpg.set_ylabel("CpGs", fontsize=6.5, labelpad=2)

            # X-axis tick labels: CG site numbers every Nth brick
            step = max(1, n_bricks // 24)
            ax_cpg.set_xticks(range(0, n_bricks, step))
            ax_cpg.set_xticklabels(
                [str(cg_sites[i]) for i in range(0, n_bricks, step)],
                rotation=45, ha="right", fontsize=6.5,
            )
            ax_cpg.set_xlabel("CG site / brick position", fontsize=7.5, labelpad=2)

            # Clean up spines on the CpG track
            ax_cpg.spines["top"].set_visible(False)
            ax_cpg.spines["right"].set_visible(False)
            ax_cpg.grid(axis="x", linestyle="--", alpha=0.3)
            ax_cpg.set_axisbelow(True)

    # ── Shared legends ────────────────────────────────────────────────────────
    heatmap_patches = [
        mpatches.Patch(color=COLOUR_FP,      label="0→1  false positive"),
        mpatches.Patch(color=COLOUR_FN,      label="1→0  missed methylation"),
        mpatches.Patch(color=COLOUR_CORRECT, label="Correct"),
    ]
    cpg_patch = mpatches.Patch(color=cpg_colour, alpha=0.80, label="No. of CpGs per brick")
    all_patches = heatmap_patches + [cpg_patch]

    fig.legend(handles=all_patches, loc="lower center", ncol=4,
               fontsize=8, bbox_to_anchor=(0.5, 0.01))

    fig.suptitle(
        "Plot 4 — Per-Brick Correctness Heatmap  +  CpG Count per Brick\n"
        "(heatmap rows = models, columns = bricks; bar track = CpG count)",
        fontsize=11, fontweight="bold",
    )


# ── CSV writers ────────────────────────────────────────────────────────────────

def write_csvs(data, all_samples, all_models, prefix):
    brick_rows   = []
    summary_rows = []

    for sample in all_samples:
        for model in all_models:
            if model not in data[sample]:
                summary_rows.append(dict(
                    sample=sample, model=model,
                    total_bricks="", matches="",
                    flips_1to0="", flips_0to1="",
                    accuracy_pct="", error_pct="",
                    sensitivity_pct="", specificity_pct="",
                    decoded_word="", note="log missing",
                ))
                continue

            p   = data[sample][model]
            cg  = p["cg_sites"]
            act = p["actual_status"]
            des = desired_vector(cg, p["desired_sites"])
            st  = compute_stats(act, des)

            for idx, (site, a, d) in enumerate(zip(cg, act, des)):
                flip = ("correct" if a == d
                        else ("1to0" if d == 1 else "0to1"))
                brick_rows.append(dict(
                    sample=sample, model=model,
                    brick_index=idx, brick_number=idx + 1,
                    cg_site=site, desired=d, actual=a,
                    correct=int(a == d), flip_type=flip,
                ))

            def fmt(v):
                return round(v, 4) if v is not None else "N/A"

            summary_rows.append(dict(
                sample=sample, model=model,
                total_bricks=st["total"], matches=st["correct"],
                flips_1to0=st["fn"], flips_0to1=st["fp"],
                accuracy_pct=round(st["acc"], 4),
                error_pct=round(100 - st["acc"], 4),
                sensitivity_pct=fmt(st["sens"]),
                specificity_pct=fmt(st["spec"]),
                mcc=round(st["mcc"], 4),
                decoded_word=p["decoded_word"],
                note="",
            ))

    brick_csv = f"{prefix}_per_brick.csv"
    with open(brick_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "sample", "model", "brick_index", "brick_number",
            "cg_site", "desired", "actual", "correct", "flip_type"])
        w.writeheader(); w.writerows(brick_rows)

    summary_csv = f"{prefix}_summary.csv"
    with open(summary_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "sample", "model", "total_bricks", "matches",
            "flips_1to0", "flips_0to1", "accuracy_pct", "error_pct",
            "sensitivity_pct", "specificity_pct", "mcc", "decoded_word", "note"])
        w.writeheader(); w.writerows(summary_rows)

    print(f"  Per-brick CSV : {brick_csv}")
    print(f"  Summary CSV   : {summary_csv}")


# ── Main ───────────────────────────────────────────────────────────────────────

def run(log_dir, prefix):
    log_dir   = Path(log_dir)
    log_files = sorted(log_dir.glob("*.log"))

    if not log_files:
        print(f"ERROR: No .log files found in '{log_dir}'"); sys.exit(1)

    data = defaultdict(dict)
    print(f"\nParsing {len(log_files)} log file(s) in '{log_dir}'...")
    for lf in log_files:
        model = detect_model(lf.name)
        if model is None:
            print(f"  SKIP (no recognised model suffix) : {lf.name}"); continue
        sample = sample_name_from_file(lf.name, model)
        parsed = parse_log(str(lf))
        if parsed is None:
            print(f"  SKIP (CG sites / status missing)  : {lf.name}"); continue
        data[sample][model] = parsed
        print(f"  OK  [{model:<12}] {sample}  "
              f"({len(parsed['cg_sites'])} bricks, "
              f"decoded='{parsed['decoded_word'] or '—'}')")

    if not data:
        print("ERROR: No logs could be parsed."); sys.exit(1)

    all_models  = sorted({m for s in data.values() for m in s})
    all_samples = sorted(data.keys())
    print(f"\nSamples  : {len(all_samples)}")
    print(f"Models   : {all_models}")

    # Aggregate model totals across all samples
    model_stats = {m: dict(total=0, correct=0, tp=0, fn=0, tn=0, fp=0,
                            acc=0, sens=None, spec=None, mcc=0.0)
                   for m in all_models}
    for sample in all_samples:
        for model in all_models:
            if model not in data[sample]:
                continue
            p   = data[sample][model]
            des = desired_vector(p["cg_sites"], p["desired_sites"])
            st  = compute_stats(p["actual_status"], des)
            for k in ("total", "correct", "tp", "fn", "tn", "fp"):
                model_stats[model][k] += st[k]
    for model in all_models:
        st = model_stats[model]
        st["acc"]  = 100 * st["correct"] / st["total"] if st["total"] else 0
        st["sens"] = (100 * st["tp"] / (st["tp"] + st["fn"])
                      if (st["tp"] + st["fn"]) else None)
        st["spec"] = (100 * st["tn"] / (st["tn"] + st["fp"])
                      if (st["tn"] + st["fp"]) else None)
        tp, fp, fn, tn = st["tp"], st["fp"], st["fn"], st["tn"]
        denom = ((tp+fp) * (tp+fn) * (tn+fp) * (tn+fn)) ** 0.5
        st["mcc"] = (tp*tn - fp*fn) / denom if denom > 0 else 0.0

    # ── Write CSVs ─────────────────────────────────────────────────────────────
    print("\nWriting CSVs...")
    write_csvs(data, all_samples, all_models, prefix)

    # ── Write PDF ──────────────────────────────────────────────────────────────
    pdf_path = f"{prefix}_figures.pdf"
    print(f"\nGenerating PDF: {pdf_path}")

    with PdfPages(pdf_path) as pdf:

        # Page 1 — overall accuracy
        fig1, ax1 = plt.subplots(figsize=(10, 6))
        plot_overall(ax1, model_stats, all_models)
        fig1.tight_layout()
        pdf.savefig(fig1, bbox_inches="tight")
        plt.close(fig1)
        print("  Page 1: overall accuracy")

        # Page 2 — per-sample accuracy
        fig2, ax2 = plt.subplots(figsize=(max(10, len(all_samples) * 2.2), 6))
        plot_per_sample(ax2, data, all_samples, all_models)
        fig2.tight_layout()
        pdf.savefig(fig2, bbox_inches="tight")
        plt.close(fig2)
        print("  Page 2: per-sample accuracy")

        # Page 3 — flip breakdown
        fig3, ax3 = plt.subplots(figsize=(max(8, len(all_models) * 1.8), 6))
        plot_flip_breakdown(ax3, data, all_samples, all_models)
        fig3.tight_layout()
        pdf.savefig(fig3, bbox_inches="tight")
        plt.close(fig3)
        print("  Page 3: flip breakdown")

        # Page 4 — per-brick heatmap + CpG track
        # Height scales with number of samples; each panel needs ~4.5 inches
        hmap_h = max(6, len(all_samples) * 4.5)
        fig4   = plt.figure(figsize=(16, hmap_h))
        plot_brick_heatmap(fig4, all_samples, all_models, data)
        pdf.savefig(fig4, bbox_inches="tight")
        plt.close(fig4)
        print("  Page 4: per-brick heatmap")

        # Page 5 — MCC scatter plot
        fig5, ax5 = plt.subplots(figsize=(max(8, len(all_models) * 1.6), 6))
        plot_mcc_scatter(ax5, data, all_samples, all_models)
        fig5.tight_layout()
        pdf.savefig(fig5, bbox_inches="tight")
        plt.close(fig5)
        print("  Page 5: MCC scatter plot")

        # PDF metadata
        d = pdf.infodict()
        d["Title"]   = "AHEAD Per-Brick Model Comparison"
        d["Subject"] = "Methylation calling accuracy across models and replicates"

    print(f"  PDF saved: {pdf_path}")

    # ── Console summary ────────────────────────────────────────────────────────
    ranked = sorted(all_models,
                    key=lambda m: model_stats[m]["mcc"], reverse=True)
    W = 84
    print("\n" + "=" * W)
    print("OVERALL MODEL RANKING  (sorted by MCC — most honest metric for imbalanced classes)")
    print("=" * W)
    print(f"{'Model':<16} {'Correct/Total':>15} {'Accuracy':>10} "
          f"{'Sensitivity':>13} {'Specificity':>13} {'MCC':>8}")
    print("-" * W)
    for i, model in enumerate(ranked):
        st = model_stats[model]
        if st["total"] == 0:
            continue
        tag    = "  <- Best"  if i == 0             else \
                 "  <- Worst" if i == len(ranked)-1  else ""
        sens_s = f"{st['sens']:>12.2f}%" if st["sens"] is not None else "         N/A"
        spec_s = f"{st['spec']:>12.2f}%" if st["spec"] is not None else "         N/A"
        frac   = f"{st['correct']}/{st['total']}"
        print(f"{model:<16} {frac:>15} {st['acc']:>9.2f}%{sens_s}{spec_s} "
              f"{st['mcc']:>7.4f}{tag}")
    print("=" * W + "\n")

    # ── Gate analysis (AND / OR combinations) ──────────────────────────────────
    run_gate_analysis(data, all_samples, prefix, model_stats, all_models)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        print("Usage: python per_brick_model_comparison.py <log_directory> [output_prefix]")
        sys.exit(1)

    log_directory = sys.argv[1]
    prefix        = sys.argv[2] if len(sys.argv) > 2 else "ahead_comparison"
    print(f"Log directory : {log_directory}")
    print(f"Output prefix : {prefix}")
    run(log_directory, prefix)
