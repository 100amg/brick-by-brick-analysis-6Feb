"""
gate_analysis.py  —  AND/OR gate add-on for per_brick_model_comparison.py
==========================================================================
Drop this file in the same directory as per_brick_model_comparison.py.

It is called automatically at the end of run() when you add ONE LINE there
(see the "HOW TO WIRE IN" section below).

What it does
------------
For every combination of 2, 3, and 4 models drawn from GATE_MODELS
(dorado_fast, dorado_hac, dorado_sup, minknow), and for both AND and OR
gate logic, it:

  1. Synthesises a virtual combined methylation call per brick per sample:
       AND gate → brick = 1  only if ALL models in the combo call it 1
       OR  gate → brick = 1  if  ANY model in the combo calls it 1

  2. Computes the same accuracy / sensitivity / specificity metrics as the
     main analysis.

  3. Produces one PDF per gate×combo-size tier (6 PDFs total):
       <prefix>_gate_OR_2way.pdf
       <prefix>_gate_AND_2way.pdf
       <prefix>_gate_OR_3way.pdf
       <prefix>_gate_AND_3way.pdf
       <prefix>_gate_OR_4way.pdf
       <prefix>_gate_AND_4way.pdf

     Each PDF contains:
       Page 1  — Overall accuracy / sensitivity / specificity bar chart
                 (one group of bars per model-combination)
       Page 2  — Per-sample accuracy bar chart
       Page 3  — Flip-type (error breakdown) stacked bar chart
       Page 4  — Per-brick heatmap + CpG track
                 (one panel per sample, rows = model combinations)

  4. Writes two summary CSVs:
       <prefix>_gate_OR_summary.csv
       <prefix>_gate_AND_summary.csv

HOW TO WIRE IN
--------------
Add these two lines at the top of per_brick_model_comparison.py:

    from gate_analysis import run_gate_analysis

And add this ONE line at the very end of the run() function, just before the
closing print:

    run_gate_analysis(data, all_samples, prefix)

That's it.

Adjusting which models are included
------------------------------------
Edit GATE_MODELS below.  minknow must match exactly the suffix used in your
log filenames (e.g. files ending _minknow.log).
"""

from itertools import combinations
import re
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpecFromSubplotSpec
import numpy as np
from collections import defaultdict

# ── Models to include in gate analysis ────────────────────────────────────────
GATE_MODELS = ["dorado_fast", "dorado_hac", "dorado_sup", "minknow"]

# ── Combo-size labels ──────────────────────────────────────────────────────────
SIZE_LABEL = {2: "2way", 3: "3way", 4: "4way"}

# ── Colour palette for combinations ───────────────────────────────────────────
# Automatically cycles if there are more combos than colours.
COMBO_PALETTE = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2",
    "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD",
    "#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#3B1F2B",
]

COLOUR_CORRECT = "#4CAF50"
COLOUR_FN      = "#FFC107"
COLOUR_FP      = "#F44336"

CPG_COUNT_BY_SITE = {
     60: 2,   84: 2,  108: 3,  132: 2,  156: 1,  180: 2,
    204: 2,  228: 1,  252: 1,  276: 2,  300: 1,  324: 2,
    348: 1,  372: 1,  396: 3,  420: 2,  444: 1,  468: 2,
    492: 1,  516: 2,  540: 2,  564: 2,  588: 3,  612: 3,
    636: 1,  660: 1,  684: 2,  708: 2,  732: 2,  756: 2,
    780: 2,  804: 2,  828: 1,  852: 2,  876: 2,  900: 1,
}


# ── Helpers ───────────────────────────────────────────────────────────────────

def _desired_vector(cg_sites, desired_sites):
    ds = set(desired_sites)
    return [1 if s in ds else 0 for s in cg_sites]


def _compute_stats(actual, desired):
    total   = len(actual)
    correct = sum(a == d for a, d in zip(actual, desired))
    tp = sum(a == 1 and d == 1 for a, d in zip(actual, desired))
    fn = sum(a == 0 and d == 1 for a, d in zip(actual, desired))
    tn = sum(a == 0 and d == 0 for a, d in zip(actual, desired))
    fp = sum(a == 1 and d == 0 for a, d in zip(actual, desired))
    acc  = 100 * correct / total if total else 0
    sens = 100 * tp / (tp + fn) if (tp + fn) else None
    spec = 100 * tn / (tn + fp) if (tn + fp) else None
    denom = ((tp+fp) * (tp+fn) * (tn+fp) * (tn+fn)) ** 0.5
    mcc   = (tp*tn - fp*fn) / denom if denom > 0 else 0.0
    return dict(total=total, correct=correct,
                tp=tp, fn=fn, tn=tn, fp=fp,
                acc=acc, sens=sens, spec=spec, mcc=mcc)


def _combo_label(combo):
    """Short label for a model combination, e.g. 'fast+hac'."""
    abbrev = {
        "dorado_fast": "fast",
        "dorado_hac":  "hac",
        "dorado_sup":  "sup",
        "minknow":     "mkn",
    }
    return "+".join(abbrev.get(m, m[:4]) for m in combo)


def _combo_colour(idx):
    return COMBO_PALETTE[idx % len(COMBO_PALETTE)]


def _apply_gate(gate, vectors):
    """
    Apply AND or OR gate across a list of 0/1 vectors (one per model).
    Returns a single 0/1 vector.
    Missing models are excluded from the vote (they don't count as 0 or 1).
    """
    if not vectors:
        return []
    length = len(vectors[0])
    result = []
    for i in range(length):
        vals = [v[i] for v in vectors]
        if gate == "AND":
            result.append(1 if all(v == 1 for v in vals) else 0)
        else:   # OR
            result.append(1 if any(v == 1 for v in vals) else 0)
    return result


# ── Synthesise gate calls for every sample and combo ──────────────────────────

def _build_gate_data(data, all_samples, gate, combo):
    """
    Returns a dict: gate_data[sample] = {
        "cg_sites":      list[int],
        "actual_status": list[int],   # gated call
        "desired_sites": list[int],
    }
    Returns None if no sample has all models in the combo present.
    """
    gate_data = {}
    for sample in all_samples:
        present = [m for m in combo if m in data[sample]]
        if len(present) < 2:          # need at least 2 for a meaningful gate
            continue

        # Use the CG sites from the first available model (they're all aligned)
        ref_model = present[0]
        p_ref     = data[sample][ref_model]
        cg_sites  = p_ref["cg_sites"]
        desired   = p_ref["desired_sites"]

        vectors = [data[sample][m]["actual_status"] for m in present
                   if len(data[sample][m]["actual_status"]) == len(cg_sites)]

        if not vectors:
            continue

        gated = _apply_gate(gate, vectors)
        gate_data[sample] = dict(
            cg_sites=cg_sites,
            actual_status=gated,
            desired_sites=desired,
            models_present=present,
        )

    return gate_data if gate_data else None


# ── Plot helpers (parallel to the main script's plots) ────────────────────────

def _plot_overall_gate(ax, combo_stats, combo_list, gate):
    """Bar chart: accuracy / sensitivity / specificity / MCC per combination."""
    metrics = [
        ("acc",  "Accuracy"),
        ("sens", "Sensitivity (1→1)"),
        ("spec", "Specificity (0→0)"),
        ("mcc",  "MCC ×100"),
    ]
    alphas  = [1.0, 0.70, 0.40, 0.60]
    bar_w   = 0.18

    x = np.arange(len(combo_list))

    for i, ((key, label), alpha) in enumerate(zip(metrics, alphas)):
        vals = []
        for c in combo_list:
            v = combo_stats[c].get(key)
            if key == "mcc":
                vals.append((v or 0) * 100)
            else:
                vals.append(v if v is not None else 0)
        offset = (i - 1.5) * bar_w
        bars = ax.bar(x + offset, vals, bar_w, label=label,
                      color=[_combo_colour(j) for j in range(len(combo_list))],
                      alpha=alpha, edgecolor="white", linewidth=0.6)
        for bar, val in zip(bars, vals):
            if abs(val) > 2:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.8,
                        f"{val:.1f}", ha="center", va="bottom",
                        fontsize=5.5, rotation=40)

    ax.set_xticks(x)
    ax.set_xticklabels([_combo_label(c).replace("+", "\n+") for c in combo_list],
                       fontsize=7)
    ax.set_ylabel("Percentage / Score (%)")
    ax.set_ylim(-20, 135)
    ax.axhline(0, color="black", linewidth=0.6, linestyle="--", alpha=0.4)
    ax.set_title(f"{gate} Gate — Overall Performance\n"
                 f"(MCC ×100 shown on same axis; MCC=0 → chance, MCC=100 → perfect)")
    metric_patches = [mpatches.Patch(facecolor="grey", alpha=a, label=l)
                      for a, (_, l) in zip(alphas, metrics)]
    ax.legend(handles=metric_patches, title="Metric", fontsize=8, loc="upper right")


def _plot_per_sample_gate(ax, gate_results, all_samples, combo_list, gate):
    """Grouped bar chart: per-sample accuracy for each combination."""
    n_combos = len(combo_list)
    bar_w    = 0.75 / n_combos
    x        = np.arange(len(all_samples))

    for i, combo in enumerate(combo_list):
        accs = []
        for sample in all_samples:
            gd = gate_results[combo].get(sample)
            if gd is None:
                accs.append(0); continue
            des = _desired_vector(gd["cg_sites"], gd["desired_sites"])
            st  = _compute_stats(gd["actual_status"], des)
            accs.append(st["acc"])
        offset = (i - (n_combos - 1) / 2) * bar_w
        bars = ax.bar(x + offset, accs, bar_w,
                      label=_combo_label(combo),
                      color=_combo_colour(i),
                      edgecolor="white", linewidth=0.5)
        for bar, val in zip(bars, accs):
            if val > 2:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.5,
                        f"{val:.0f}%", ha="center", va="bottom",
                        fontsize=5.5, rotation=40)

    short_samples = []
    for s in all_samples:
        m = re.search(r"(AHEAD[_.\d]+)", s, re.IGNORECASE)
        short_samples.append(m.group(1) if m else s[-30:])

    ax.set_xticks(x)
    ax.set_xticklabels(short_samples, rotation=18, ha="right", fontsize=8)
    ax.set_ylabel("Accuracy (%)")
    ax.set_ylim(0, 125)
    ax.set_title(f"{gate} Gate — Per-Sample Accuracy")
    ax.legend(title="Combination", fontsize=7,
              bbox_to_anchor=(1.01, 1), loc="upper left", ncol=max(1, n_combos // 8))


def _plot_flip_gate(ax, gate_results, all_samples, combo_list, gate):
    """Stacked bar: correct | 1→0 | 0→1 per combination."""
    totals = {c: dict(correct=0, fn=0, fp=0) for c in combo_list}

    for combo in combo_list:
        for sample in all_samples:
            gd = gate_results[combo].get(sample)
            if gd is None:
                continue
            des = _desired_vector(gd["cg_sites"], gd["desired_sites"])
            st  = _compute_stats(gd["actual_status"], des)
            totals[combo]["correct"] += st["correct"]
            totals[combo]["fn"]      += st["fn"]
            totals[combo]["fp"]      += st["fp"]

    x     = np.arange(len(combo_list))
    bar_w = 0.55

    corrects = [totals[c]["correct"] for c in combo_list]
    fns      = [totals[c]["fn"]      for c in combo_list]
    fps      = [totals[c]["fp"]      for c in combo_list]

    ax.bar(x, corrects, bar_w, label="Correct",
           color=[_combo_colour(i) for i in range(len(combo_list))],
           alpha=0.88, edgecolor="white", linewidth=0.5)
    ax.bar(x, fns, bar_w, bottom=corrects,
           label="1→0  missed methylation",
           color=COLOUR_FN, alpha=0.90, edgecolor="white", linewidth=0.5)
    ax.bar(x, fps, bar_w,
           bottom=[c + f for c, f in zip(corrects, fns)],
           label="0→1  false positive",
           color=COLOUR_FP, alpha=0.90, edgecolor="white", linewidth=0.5)

    for i, (c, fn_v, fp_v) in enumerate(zip(corrects, fns, fps)):
        total = c + fn_v + fp_v
        if total == 0:
            continue
        ax.text(i, c / 2, f"{100*c/total:.0f}%",
                ha="center", va="center",
                fontsize=7.5, fontweight="bold", color="white")
        if fn_v > 0:
            ax.text(i, c + fn_v / 2, f"{100*fn_v/total:.0f}%",
                    ha="center", va="center", fontsize=7, color="black")
        if fp_v > 0:
            ax.text(i, c + fn_v + fp_v / 2, f"{100*fp_v/total:.0f}%",
                    ha="center", va="center", fontsize=7, color="white")

    ax.set_xticks(x)
    ax.set_xticklabels([_combo_label(c).replace("+", "\n+") for c in combo_list],
                       fontsize=7)
    ax.set_ylabel("Bricks (all samples)")
    ax.set_title(f"{gate} Gate — Error Type Breakdown (all samples combined)")
    ax.legend(fontsize=8, loc="upper right")


def _plot_mcc_scatter_gate(ax, gate_results, all_samples, combo_list, gate):
    """
    Scatter plot: x = combination (jittered), y = MCC per sample.
    Mean bar shown per combination.  Mean MCC label placed ABOVE the bar.
    Dashed baseline at MCC = 0.
    """
    np.random.seed(42)

    n_combos = len(combo_list)

    for ci, combo in enumerate(combo_list):
        mccs    = []
        samples = []
        for sample in all_samples:
            gd = gate_results[combo].get(sample)
            if gd is None:
                continue
            des = _desired_vector(gd["cg_sites"], gd["desired_sites"])
            st  = _compute_stats(gd["actual_status"], des)
            mccs.append(st["mcc"])
            samples.append(sample)

        if not mccs:
            continue

        jitter = np.random.uniform(-0.12, 0.12, size=len(mccs))
        xs     = [ci + j for j in jitter]

        ax.scatter(xs, mccs,
                   color=_combo_colour(ci), s=60, zorder=3,
                   alpha=0.85, edgecolors="white", linewidths=0.6)

        # --- sample date labels: always below each dot so they don't
        #     crowd the mean-bar label that sits above the bar
        for x, mcc_val, sample in zip(xs, mccs, samples):
            m = re.search(r"AHEAD[_.](\d+[_.]\d+[_.]\d+)", sample, re.IGNORECASE)
            lbl = m.group(1).replace("_", "/") if m else ""
            ax.annotate(
                lbl, (x, mcc_val),
                textcoords="offset points", xytext=(0, -8),
                ha="center", va="top",
                fontsize=5.5, color="grey",
            )

        mean_mcc = np.mean(mccs)
        bar_x0   = ci - 0.25
        bar_x1   = ci + 0.25

        # Draw the mean bar
        ax.plot([bar_x0, bar_x1], [mean_mcc, mean_mcc],
                color=_combo_colour(ci), linewidth=2.5, zorder=4,
                solid_capstyle="round")

        # Mean MCC label ABOVE the bar (centred), not to the side
        ax.annotate(
            f"{mean_mcc:.3f}",
            (ci, mean_mcc),
            textcoords="offset points",
            xytext=(0, 7),          # 7 pt above the bar
            ha="center", va="bottom",
            fontsize=7.5, fontweight="bold",
            color=_combo_colour(ci),
        )

    ax.axhline(0, color="black", linewidth=1.0, linestyle="--", alpha=0.5,
               label="MCC = 0  (chance)")
    ax.axhline(1, color="green", linewidth=0.6, linestyle=":", alpha=0.4,
               label="MCC = 1  (perfect)")

    ax.set_xticks(range(n_combos))
    ax.set_xticklabels([_combo_label(c).replace("+", "\n+") for c in combo_list],
                       fontsize=7)
    ax.set_ylabel("MCC  (Matthews Correlation Coefficient)")

    # Extra headroom at the top so above-bar labels are never clipped
    ax.set_ylim(-1.05, 1.20)

    # Slightly wider x-range so edge combo labels aren't cut off
    ax.set_xlim(-0.75, n_combos - 0.25)

    ax.set_title(
        f"{gate} Gate — MCC per Combination per Sample\n"
        "Dots = individual replicates  |  Bar = mean across replicates\n"
        "MCC = 0 → no better than chance  |  MCC = 1 → perfect",
        fontsize=10,
    )
    ax.legend(fontsize=8, loc="lower right")
    ax.yaxis.grid(True, linestyle="--", alpha=0.4)
    ax.set_axisbelow(True)



def _plot_heatmap_gate(fig, all_samples, combo_list, gate_results, gate):
    """
    Heatmap: rows = model combinations, columns = bricks.
    One panel per sample.  Same colour scheme as main heatmap.
    CpG count bar track beneath each panel.
    """
    from matplotlib.colors import ListedColormap, BoundaryNorm

    n_samples  = len(all_samples)
    cmap       = ListedColormap([COLOUR_FP, COLOUR_FN, COLOUR_CORRECT])
    norm       = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    cpg_colour = "#5B8DB8"

    outer_gs = fig.add_gridspec(
        n_samples, 1, hspace=0.75,
        top=0.91, bottom=0.07, left=0.14, right=0.97,
    )

    for row_idx, sample in enumerate(all_samples):
        inner_gs = GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer_gs[row_idx],
            height_ratios=[4, 1], hspace=0.10,
        )
        ax_heat = fig.add_subplot(inner_gs[0])
        ax_cpg  = fig.add_subplot(inner_gs[1])

        # Find max bricks across all combos for this sample
        n_bricks = 0
        cg_sites = None
        for combo in combo_list:
            gd = gate_results[combo].get(sample)
            if gd and len(gd["cg_sites"]) > n_bricks:
                n_bricks = len(gd["cg_sites"])
                cg_sites = gd["cg_sites"]

        if n_bricks == 0:
            ax_heat.set_visible(False)
            ax_cpg.set_visible(False)
            continue

        matrix = np.full((len(combo_list), n_bricks), np.nan)

        for ci, combo in enumerate(combo_list):
            gd = gate_results[combo].get(sample)
            if gd is None:
                continue
            des = _desired_vector(gd["cg_sites"], gd["desired_sites"])
            act = gd["actual_status"]
            for bi, (a, d) in enumerate(zip(act, des)):
                matrix[ci, bi] = 2 if a == d else (0 if d == 0 else 1)

        ax_heat.imshow(np.ma.masked_invalid(matrix),
                       cmap=cmap, norm=norm,
                       aspect="auto", interpolation="nearest")
        ax_heat.set_yticks(range(len(combo_list)))
        ax_heat.set_yticklabels([_combo_label(c) for c in combo_list], fontsize=6.5)
        ax_heat.set_xticks([])
        ax_heat.tick_params(bottom=False)

        m_t   = re.search(r"(AHEAD[_.\d]+)", sample, re.IGNORECASE)
        short = m_t.group(1) if m_t else sample[-40:]
        ax_heat.set_title(f"Sample: {short}", fontsize=9, pad=3)

        # CpG count track
        if cg_sites:
            cpg_counts = [CPG_COUNT_BY_SITE.get(s, 0) for s in cg_sites]
            bx         = np.arange(n_bricks)
            ax_cpg.bar(bx, cpg_counts, width=0.85,
                       color=cpg_colour, alpha=0.80, linewidth=0)
            max_cpg = max(cpg_counts) if cpg_counts else 1
            for b, val in zip(bx, cpg_counts):
                if val == 0:
                    continue
                ax_cpg.text(b, val / 2, str(val),
                            ha="center", va="center",
                            fontsize=5.5, fontweight="bold",
                            color="white" if val >= 2 else "black")
            ax_cpg.set_xlim(ax_heat.get_xlim())
            ax_cpg.set_ylim(0, max_cpg + 0.6)
            ax_cpg.set_yticks([1, max_cpg] if max_cpg > 1 else [1])
            ax_cpg.set_yticklabels(
                [str(v) for v in ([1, max_cpg] if max_cpg > 1 else [1])],
                fontsize=6)
            ax_cpg.set_ylabel("CpGs", fontsize=6.5, labelpad=2)
            step = max(1, n_bricks // 24)
            ax_cpg.set_xticks(range(0, n_bricks, step))
            ax_cpg.set_xticklabels(
                [str(cg_sites[i]) for i in range(0, n_bricks, step)],
                rotation=45, ha="right", fontsize=6.5)
            ax_cpg.set_xlabel("CG site / brick position", fontsize=7.5, labelpad=2)
            ax_cpg.spines["top"].set_visible(False)
            ax_cpg.spines["right"].set_visible(False)
            ax_cpg.grid(axis="x", linestyle="--", alpha=0.3)
            ax_cpg.set_axisbelow(True)

    # Shared legends
    hm_patches = [
        mpatches.Patch(color=COLOUR_FP,      label="0→1  false positive"),
        mpatches.Patch(color=COLOUR_FN,      label="1→0  missed methylation"),
        mpatches.Patch(color=COLOUR_CORRECT, label="Correct"),
        mpatches.Patch(color=cpg_colour, alpha=0.8, label="No. of CpGs"),
    ]
    fig.legend(handles=hm_patches, loc="lower center", ncol=4,
               fontsize=8, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle(
        f"Plot 4 — {gate} Gate  Per-Brick Heatmap + CpG Track\n"
        f"(rows = model combinations, columns = bricks)",
        fontsize=11, fontweight="bold",
    )


# ── Build one full PDF for a given gate × combo-size ──────────────────────────

def _write_gate_pdf(gate, combo_size, data, all_samples, prefix):
    # Filter to combos where at least one sample has all required models
    available_models = [m for m in GATE_MODELS
                        if any(m in data[s] for s in all_samples)]
    if len(available_models) < combo_size:
        print(f"  SKIP {gate} {combo_size}-way: only {len(available_models)} "
              f"gate models available")
        return None

    all_combos = list(combinations(available_models, combo_size))

    # Pre-compute gated vectors for every combo
    gate_results = {}
    valid_combos = []
    for combo in all_combos:
        gd = _build_gate_data(data, all_samples, gate, combo)
        if gd is not None:
            gate_results[combo] = gd
            valid_combos.append(combo)

    if not valid_combos:
        print(f"  SKIP {gate} {combo_size}-way: no valid combinations")
        return None

    label     = SIZE_LABEL[combo_size]
    pdf_path  = f"{prefix}_gate_{gate}_{label}.pdf"
    print(f"\n  Generating {pdf_path}  ({len(valid_combos)} combinations)...")

    # Aggregate stats per combo across all samples
    combo_stats = {}
    for combo in valid_combos:
        agg = dict(total=0, correct=0, tp=0, fn=0, tn=0, fp=0)
        for sample in all_samples:
            gd = gate_results[combo].get(sample)
            if gd is None:
                continue
            des = _desired_vector(gd["cg_sites"], gd["desired_sites"])
            st  = _compute_stats(gd["actual_status"], des)
            for k in ("total", "correct", "tp", "fn", "tn", "fp"):
                agg[k] += st[k]
        agg["acc"]  = 100 * agg["correct"] / agg["total"] if agg["total"] else 0
        agg["sens"] = (100 * agg["tp"] / (agg["tp"] + agg["fn"])
                       if (agg["tp"] + agg["fn"]) else None)
        agg["spec"] = (100 * agg["tn"] / (agg["tn"] + agg["fp"])
                       if (agg["tn"] + agg["fp"]) else None)
        tp, fp, fn, tn = agg["tp"], agg["fp"], agg["fn"], agg["tn"]
        denom = ((tp+fp) * (tp+fn) * (tn+fp) * (tn+fn)) ** 0.5
        agg["mcc"] = (tp*tn - fp*fn) / denom if denom > 0 else 0.0
        combo_stats[combo] = agg

    with PdfPages(pdf_path) as pdf:

        # Page 1 — overall
        n_combos = len(valid_combos)
        fig1, ax1 = plt.subplots(figsize=(max(10, n_combos * 1.4), 6))
        _plot_overall_gate(ax1, combo_stats, valid_combos, gate)
        fig1.tight_layout()
        pdf.savefig(fig1, bbox_inches="tight"); plt.close(fig1)

        # Page 2 — per-sample
        fig2, ax2 = plt.subplots(
            figsize=(max(10, len(all_samples) * max(2.5, n_combos * 0.5)), 6))
        _plot_per_sample_gate(ax2, gate_results, all_samples, valid_combos, gate)
        fig2.tight_layout()
        pdf.savefig(fig2, bbox_inches="tight"); plt.close(fig2)

        # Page 3 — flip breakdown
        fig3, ax3 = plt.subplots(figsize=(max(8, n_combos * 1.4), 6))
        _plot_flip_gate(ax3, gate_results, all_samples, valid_combos, gate)
        fig3.tight_layout()
        pdf.savefig(fig3, bbox_inches="tight"); plt.close(fig3)

        # Page 4 — heatmap
        hmap_h = max(6, len(all_samples) * max(4.5, n_combos * 0.7))
        fig4   = plt.figure(figsize=(16, hmap_h))
        _plot_heatmap_gate(fig4, all_samples, valid_combos, gate_results, gate)
        pdf.savefig(fig4, bbox_inches="tight"); plt.close(fig4)

        # Page 5 — MCC scatter
        fig5, ax5 = plt.subplots(figsize=(max(8, n_combos * 1.4), 6))
        _plot_mcc_scatter_gate(ax5, gate_results, all_samples, valid_combos, gate)
        fig5.tight_layout()
        pdf.savefig(fig5, bbox_inches="tight"); plt.close(fig5)

        d = pdf.infodict()
        d["Title"]   = f"AHEAD Gate Analysis — {gate} {label}"
        d["Subject"] = "Methylation calling gate logic comparison"

    print(f"    Saved: {pdf_path}")
    return combo_stats, valid_combos, gate_results


# ── Write summary CSVs for AND and OR ─────────────────────────────────────────

def _write_gate_csv(gate, all_results_by_size, all_samples, prefix):
    """
    all_results_by_size: dict {combo_size → (combo_stats, valid_combos, gate_results)}
    """
    rows = []
    for combo_size, result in all_results_by_size.items():
        if result is None:
            continue
        combo_stats, valid_combos, gate_results = result
        for combo in valid_combos:
            st = combo_stats[combo]
            # Also per-sample rows
            for sample in all_samples:
                gd = gate_results[combo].get(sample)
                if gd is None:
                    rows.append(dict(
                        gate=gate, combo_size=combo_size,
                        combination=_combo_label(combo),
                        sample=sample,
                        models_present="",
                        total_bricks="", matches="",
                        flips_1to0="", flips_0to1="",
                        accuracy_pct="", error_pct="",
                        sensitivity_pct="", specificity_pct="",
                        mcc="",
                        note="no data",
                    ))
                    continue
                des  = _desired_vector(gd["cg_sites"], gd["desired_sites"])
                s_st = _compute_stats(gd["actual_status"], des)
                def fmt(v):
                    return round(v, 4) if v is not None else "N/A"
                rows.append(dict(
                    gate=gate,
                    combo_size=combo_size,
                    combination=_combo_label(combo),
                    sample=sample,
                    models_present=",".join(gd.get("models_present", combo)),
                    total_bricks=s_st["total"],
                    matches=s_st["correct"],
                    flips_1to0=s_st["fn"],
                    flips_0to1=s_st["fp"],
                    accuracy_pct=round(s_st["acc"], 4),
                    error_pct=round(100 - s_st["acc"], 4),
                    sensitivity_pct=fmt(s_st["sens"]),
                    specificity_pct=fmt(s_st["spec"]),
                    mcc=round(s_st["mcc"], 4),
                    note="",
                ))

    csv_path = f"{prefix}_gate_{gate}_summary.csv"
    fields   = [
        "gate", "combo_size", "combination", "sample", "models_present",
        "total_bricks", "matches", "flips_1to0", "flips_0to1",
        "accuracy_pct", "error_pct", "sensitivity_pct", "specificity_pct",
        "mcc", "note",
    ]
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  Gate {gate} CSV → {csv_path}")


# ── Console ranking table ──────────────────────────────────────────────────────

def _print_gate_ranking(gate, all_results_by_size):
    """
    Four ranked tables per gate: accuracy, sensitivity, specificity, and MCC.
    MCC is the primary headline metric — it uses all four confusion-matrix
    cells and cannot be inflated by predicting the majority class.
    """
    W = 88

    all_rows = []
    for combo_size, result in all_results_by_size.items():
        if result is None:
            continue
        combo_stats, valid_combos, _ = result
        for combo in valid_combos:
            all_rows.append((combo, combo_size, combo_stats[combo]))

    if not all_rows:
        return

    print(f"\n{'='*W}")
    print(f"GATE ANALYSIS RANKING  —  {gate} gate")
    print(f"{'='*W}")
    print(
        f"  NOTE: AND gates improve specificity (fewer false positives) but\n"
        f"  reduce sensitivity (more missed methylations). OR gates do the\n"
        f"  reverse. MCC is the primary ranking metric — it accounts for all\n"
        f"  four confusion-matrix cells and is not inflated by class imbalance.\n"
    )

    col = (f"{'Combo':<22} {'Size':>5} {'Correct/Total':>15} "
           f"{'Accuracy':>10} {'Sensitivity':>13} {'Specificity':>13} {'MCC':>8}")
    sep = "-" * W

    def fmt_row(combo, combo_size, st, tag=""):
        frac   = f"{st['correct']}/{st['total']}"
        sens_s = f"{st['sens']:>12.2f}%" if st["sens"] is not None else "         N/A"
        spec_s = f"{st['spec']:>12.2f}%" if st["spec"] is not None else "         N/A"
        mcc_s  = f"{st['mcc']:>7.4f}"
        return (f"{_combo_label(combo):<22} {combo_size:>5} {frac:>15} "
                f"{st['acc']:>9.2f}%{sens_s}{spec_s} {mcc_s}{tag}")

    # ── Ranked by MCC (primary — most honest single metric) ───────────────────
    print(f"  Ranked by MCC  (primary metric — balances all four confusion cells)")
    print(col); print(sep)
    sorted_mcc = sorted(all_rows, key=lambda r: r[2]["mcc"], reverse=True)
    best_mcc   = sorted_mcc[0][2]["mcc"]
    worst_mcc  = sorted_mcc[-1][2]["mcc"]
    for combo, combo_size, st in sorted_mcc:
        tag = "  <- highest" if st["mcc"] == best_mcc  else \
              "  <- lowest"  if st["mcc"] == worst_mcc else ""
        print(fmt_row(combo, combo_size, st, tag))

    # ── Ranked by accuracy ────────────────────────────────────────────────────
    print(f"\n  Ranked by ACCURACY  (overall correct bricks — misleading if "
          f"classes are imbalanced)")
    print(col); print(sep)
    sorted_acc = sorted(all_rows, key=lambda r: r[2]["acc"], reverse=True)
    best_acc   = sorted_acc[0][2]["acc"]
    worst_acc  = sorted_acc[-1][2]["acc"]
    for combo, combo_size, st in sorted_acc:
        tag = "  <- highest" if st["acc"] == best_acc  else \
              "  <- lowest"  if st["acc"] == worst_acc else ""
        print(fmt_row(combo, combo_size, st, tag))

    # ── Ranked by sensitivity ─────────────────────────────────────────────────
    print(f"\n  Ranked by SENSITIVITY  (correctly detected methylated bricks; "
          f"fewer 1→0 misses)")
    print(col); print(sep)
    sorted_sens = sorted(
        all_rows,
        key=lambda r: r[2]["sens"] if r[2]["sens"] is not None else -1,
        reverse=True,
    )
    best_s  = sorted_sens[0][2]["sens"]
    worst_s = sorted_sens[-1][2]["sens"]
    for combo, combo_size, st in sorted_sens:
        v   = st["sens"]
        tag = "  <- highest" if v == best_s  else \
              "  <- lowest"  if v == worst_s else ""
        print(fmt_row(combo, combo_size, st, tag))

    # ── Ranked by specificity ─────────────────────────────────────────────────
    print(f"\n  Ranked by SPECIFICITY  (correctly left unmethylated; "
          f"fewer 0→1 false positives)")
    print(col); print(sep)
    sorted_spec = sorted(
        all_rows,
        key=lambda r: r[2]["spec"] if r[2]["spec"] is not None else -1,
        reverse=True,
    )
    best_sp  = sorted_spec[0][2]["spec"]
    worst_sp = sorted_spec[-1][2]["spec"]
    for combo, combo_size, st in sorted_spec:
        v   = st["spec"]
        tag = "  <- highest" if v == best_sp  else \
              "  <- lowest"  if v == worst_sp else ""
        print(fmt_row(combo, combo_size, st, tag))

    print("=" * W)

def plot_top5_mcc_comparison(ax_or_fig,
                             data, all_samples,
                             model_stats, all_models,
                             all_results_by_size_AND,
                             all_results_by_size_OR,
                             top_n=5):
    """
    Collects MCC scores for every entry in the full candidate pool:
        • 4  single models          (from model_stats / data)
        • 6  2-way AND + 6  2-way OR  combinations
        • 4  3-way AND + 4  3-way OR  combinations
        • 1  4-way AND + 1  4-way OR  combination
    = 26 candidates total.

    Ranks all candidates by mean MCC, keeps the top `top_n`, then draws the
    same scatter style used throughout the analysis:
        • one dot per sample replicate (jittered)
        • horizontal bar at the mean
        • bold mean MCC label centred above the bar
        • sample date labels 12 pt below each dot

    Parameters
    ----------
    ax_or_fig   : a pre-created matplotlib Axes object
    data        : the same defaultdict(dict) built by run()
    all_samples : sorted list of sample keys
    model_stats : dict of aggregated stats per single model (from run())
    all_models  : list of single-model names
    all_results_by_size_AND / _OR :
                  dicts {combo_size → (combo_stats, valid_combos, gate_results)}
                  as returned by _write_gate_pdf()
    top_n       : how many top candidates to display (default 5)

    How to wire in
    --------------
    Call this at the end of run_gate_analysis(), passing the accumulated
    all_results_by_size dicts for both gates, e.g.:

        and_results = {}
        or_results  = {}
        for gate in ("AND", "OR"):
            store = and_results if gate == "AND" else or_results
            for combo_size in combo_sizes:
                store[combo_size] = _write_gate_pdf(gate, combo_size, ...)

        fig, ax = plt.subplots(figsize=(12, 6))
        plot_top5_mcc_comparison(
            ax, data, all_samples,
            model_stats, all_models,
            and_results, or_results,
        )
    """
    import re
    import numpy as np

    # ── 1. Build the full candidate pool ──────────────────────────────────────
    # Each entry: {"label": str, "colour": str, "per_sample_mccs": list[float]}

    candidates = []

    # --- single models --------------------------------------------------------
    for model in all_models:
        mccs = []
        for sample in all_samples:
            if model not in data[sample]:
                continue
            p   = data[sample][model]
            des = _desired_vector(p["cg_sites"], p["desired_sites"])
            st  = _compute_stats(p["actual_status"], des)
            mccs.append(st["mcc"])
        if mccs:
            candidates.append({
                "label":           model.replace("dorado_", "").replace("_", "\n"),
                "full_label":      model,
                "colour":          _combo_colour(all_models.index(model)),
                "per_sample_mccs": mccs,
                "mean_mcc":        np.mean(mccs),
                "kind":            "single",
            })

    # --- gate combinations (AND and OR) ---------------------------------------
    gate_colour_offset = len(all_models)   # so combo colours don't clash with singles

    for gate_label, results_by_size in (("AND", all_results_by_size_AND),
                                         ("OR",  all_results_by_size_OR)):
        colour_idx = gate_colour_offset
        for combo_size, result in results_by_size.items():
            if result is None:
                continue
            combo_stats, valid_combos, gate_results = result
            for combo in valid_combos:
                mccs = []
                for sample in all_samples:
                    gd = gate_results[combo].get(sample)
                    if gd is None:
                        continue
                    des = _desired_vector(gd["cg_sites"], gd["desired_sites"])
                    st  = _compute_stats(gd["actual_status"], des)
                    mccs.append(st["mcc"])
                if mccs:
                    lbl = _combo_label(combo).replace("+", "\n+")
                    candidates.append({
                        "label":           f"{lbl}\n({gate_label})",
                        "full_label":      f"{_combo_label(combo)} [{gate_label}]",
                        "colour":          _combo_colour(colour_idx),
                        "per_sample_mccs": mccs,
                        "mean_mcc":        np.mean(mccs),
                        "kind":            gate_label,
                    })
                colour_idx += 1

    if not candidates:
        return

    # ── 2. Rank and keep top N ─────────────────────────────────────────────────
    ranked     = sorted(candidates, key=lambda c: c["mean_mcc"], reverse=True)
    top        = ranked[:top_n]

    # ── 3. Draw the scatter ────────────────────────────────────────────────────
    ax = ax_or_fig
    np.random.seed(42)

    for ci, cand in enumerate(top):
        mccs    = cand["per_sample_mccs"]
        colour  = cand["colour"]
        jitter  = np.random.uniform(-0.12, 0.12, size=len(mccs))
        xs      = [ci + j for j in jitter]

        # individual sample dots
        ax.scatter(xs, mccs,
                   color=colour, s=70, zorder=3,
                   alpha=0.85, edgecolors="white", linewidths=0.7)

        # sample date labels — always 12 pt below each dot
        for x, mcc_val, sample in zip(xs, mccs, all_samples):
            m = re.search(r"AHEAD[_.](\d+[_.]\d+[_.]\d+)", sample, re.IGNORECASE)
            lbl = m.group(1).replace("_", "/") if m else ""
            ax.annotate(
                lbl, (x, mcc_val),
                textcoords="offset points", xytext=(20, 0),
                ha="center", va="top",
                fontsize=5.5, color="grey",
            )

        # mean bar
        mean_mcc = cand["mean_mcc"]
        ax.plot([ci - 0.25, ci + 0.25], [mean_mcc, mean_mcc],
                color=colour, linewidth=2.5, zorder=4,
                solid_capstyle="round")

        # mean MCC label centred above the bar
        ax.annotate(
            f"{mean_mcc:.3f}",
            (ci, mean_mcc),
            textcoords="offset points", xytext=(0, 7),
            ha="center", va="bottom",
            fontsize=8, fontweight="bold",
            color=colour,
        )

    # ── 4. Reference lines and axes ───────────────────────────────────────────
    ax.axhline(0, color="black", linewidth=1.0, linestyle="--", alpha=0.5,
               label="MCC = 0  (chance)")
    ax.axhline(1, color="green", linewidth=0.6, linestyle=":", alpha=0.4,
               label="MCC = 1  (perfect)")

    ax.set_xticks(range(top_n))
    ax.set_xticklabels([c["label"] for c in top], fontsize=8)
    ax.set_ylabel("MCC  (Matthews Correlation Coefficient)")
    ax.set_ylim(-1.05, 1.25)
    ax.set_xlim(-0.75, top_n - 0.25)

    # rank badges on x-axis labels
    ax.set_title(
        f"Top {top_n} Performing Combinations — MCC per Sample\n"
        "Ranked by mean MCC across all replicates  |  "
        "Singles + AND/OR gate combos compared\n"
        "Dots = individual replicates  |  Bar = mean  |  "
        "MCC = 0 → chance  |  MCC = 1 → perfect",
        fontsize=10,
    )

    # kind badges in legend
    kind_patches = [
        mpatches.Patch(color="#888", label="Single model"),
        mpatches.Patch(color="#555", label="AND gate combo"),
        mpatches.Patch(color="#333", label="OR gate combo"),
    ]
    ax.legend(handles=kind_patches + [
        plt.Line2D([0], [0], color="black", linewidth=1, linestyle="--",
                   label="MCC = 0  (chance)"),
        plt.Line2D([0], [0], color="green", linewidth=0.8, linestyle=":",
                   label="MCC = 1  (perfect)"),
    ], fontsize=8, loc="lower right")

    ax.yaxis.grid(True, linestyle="--", alpha=0.4)
    ax.set_axisbelow(True)

    # ── 5. Print ranking to console ───────────────────────────────────────────
    print(f"\n{'='*70}")
    print(f"TOP {top_n} CANDIDATES BY MEAN MCC  (all singles + AND/OR combos)")
    print(f"{'='*70}")
    print(f"  {'Rank':<5} {'Candidate':<28} {'Type':<8} {'Mean MCC':>9}  Per-sample MCCs")
    print(f"  {'-'*70}")
    for rank, cand in enumerate(top, 1):
        per = "  ".join(f"{v:.3f}" for v in cand["per_sample_mccs"])
        print(f"  {rank:<5} {cand['full_label']:<28} {cand['kind']:<8} "
              f"{cand['mean_mcc']:>9.4f}  [{per}]")
    print(f"{'='*70}\n")


# ── Public entry point ─────────────────────────────────────────────────────────

# ── WIRING SNIPPET ────────────────────────────────────────────────────────────
# Replace the existing run_gate_analysis() loop with this version.
# The only change is that AND/OR results are stored separately so they can
# be passed to plot_top5_mcc_comparison() at the end.
#
# Also add this import at the top of gate_analysis.py (if not already present):
#   import matplotlib.patches as mpatches
#   import matplotlib.pyplot as plt   (already present)

def run_gate_analysis(data, all_samples, prefix,
                      model_stats=None, all_models=None):
    """
    Call this from the end of run() in per_brick_model_comparison.py:

        run_gate_analysis(data, all_samples, prefix, model_stats, all_models)

    `data`        : the same defaultdict(dict) built by run()
    `all_samples` : sorted list of sample keys
    `prefix`      : output file prefix (same as main script)
    `model_stats` : dict of aggregated stats per single model (optional,
                    needed for the top-5 comparison plot)
    `all_models`  : list of single-model names (optional, same caveat)
    """
    print("\n" + "=" * 70)
    print("GATE ANALYSIS  (AND / OR  ×  2-way / 3-way / 4-way)")
    print("=" * 70)

    combo_sizes = [s for s in [2, 3, 4] if s <= len(GATE_MODELS)]

    # Store results for both gates so we can feed the top-5 plot
    and_results = {}
    or_results  = {}

    for gate in ("OR", "AND"):
        print(f"\n── {gate} gate ──────────────────────────────────────────────")
        store = and_results if gate == "AND" else or_results
        all_results_by_size = {}
        for combo_size in combo_sizes:
            result = _write_gate_pdf(gate, combo_size, data, all_samples, prefix)
            all_results_by_size[combo_size] = result
            store[combo_size] = result

        _print_gate_ranking(gate, all_results_by_size)
        _write_gate_csv(gate, all_results_by_size, all_samples, prefix)

    # ── Top-5 comparison plot (requires model_stats + all_models) ─────────────
    if model_stats is not None and all_models is not None:
        fig_top, ax_top = plt.subplots(figsize=(14, 6))
        plot_top5_mcc_comparison(
            ax_top,
            data, all_samples,
            model_stats, all_models,
            and_results, or_results,
            top_n=5,
        )
        fig_top.tight_layout()
        top5_path = f"{prefix}_top5_mcc_comparison.pdf"
        fig_top.savefig(top5_path, bbox_inches="tight")
        plt.close(fig_top)
        print(f"\n  Top-5 comparison plot saved: {top5_path}")

    print("\nGate analysis complete.")


# ── Also update the call-site in updated_AHEAD_model_comparison.py ────────────
# In run(), change the last line from:
#
#   run_gate_analysis(data, all_samples, prefix)
#
# to:
#
#   run_gate_analysis(data, all_samples, prefix, model_stats, all_models)