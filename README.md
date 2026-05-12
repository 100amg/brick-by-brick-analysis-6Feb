# Cross-Model ASCII Accuracy Comparison and Gate-Based Ensemble Analysis

Pipeline for cross-model methylation accuracy comparison, ASCII reconstruction benchmarking, and logical gate-based ensemble analysis across nanopore methylation-calling workflows.

---

# Overview

This workflow compares methylation-calling performance across multiple nanopore methylation-calling models using ASCII reconstruction logs generated from previous methylation-decoding analyses.

The workflow evaluates:

* per-brick methylation accuracy
* binary reconstruction accuracy
* ASCII decoding consistency
* false positive and false negative rates
* sensitivity and specificity
* Matthews Correlation Coefficient (MCC)
* ensemble model performance using logical AND/OR gate combinations

The workflow integrates outputs generated from:

* Dorado FAST methylation calling
* Dorado HAC methylation calling
* Dorado SUP methylation calling
* Guppy HAC methylation calling
* DeepMod2 methylation calling
* MinKNOW methylation outputs

The workflow corresponds specifically to:

| Step | Script                              | Purpose                       |
| ---- | ----------------------------------- | ----------------------------- |
| 1    | `updated_AHEAD_model_comparison.py` | Cross-model accuracy analysis |
| 2    | `gate_analysis.py`                  | Ensemble AND/OR gate analysis |

---

# Workflow Purpose

The workflow compares reconstructed methylation patterns across multiple models and evaluates:

* how accurately each model reconstructs encoded methylation information
* whether combining models improves decoding performance
* which model or ensemble yields the most reliable methylation detection

The workflow operates on ASCII log files generated from the methylation-decoding pipeline.

---

# Repository Structure

```text id="jlwmgt"
cross-model-ascii-comparison/
│
├── README.md
│
├── updated_AHEAD_model_comparison.py
├── gate_analysis.py
├── outputs/
│
└── .gitignore
```

---

# Input Requirements

The workflow requires ASCII log files generated from the methylation decoding pipeline.

Each log corresponds to:

* one sample
* one methylation-calling model

Example filenames:

```text id="jlwm3c"
ASCII_Log_AHEAD_10_11_25_dorado_fast.log
ASCII_Log_AHEAD_10_11_25_dorado_hac.log
ASCII_Log_AHEAD_10_11_25_dorado_sup.log
ASCII_Log_AHEAD_10_11_25_guppy_hac.log
ASCII_Log_AHEAD_10_11_25_deepmod2.log
```

The workflow automatically parses model names from log suffixes.

---

# Supported Models

Recognized models include:

| Model         |
| ------------- |
| `dorado_fast` |
| `dorado_hac`  |
| `dorado_sup`  |
| `guppy_hac`   |
| `deepmod2`    |
| `minknow`     |

Additional models can be added by editing:

```text id="jlwmjg"
KNOWN_MODELS
```

---

# Recommended Directory Structure

```text id="jlwm5w"
final_model_comparison/
├── ASCII_Log_AHEAD_10_11_25_dorado_fast.log
├── ASCII_Log_AHEAD_10_11_25_dorado_hac.log
├── ASCII_Log_AHEAD_10_11_25_dorado_sup.log
├── ASCII_Log_AHEAD_10_11_25_guppy_hac.log
├── ASCII_Log_AHEAD_10_11_25_deepmod2.log
├── updated_AHEAD_model_comparison.py
├── gate_analysis.py
└── outputs/
```

---

# Main Workflow Scripts

## Main Analysis

```text id="jlwmj8"
updated_AHEAD_model_comparison.py
```

---

## Gate Analysis Module

```text id="jlwmbo"
gate_analysis.py
```

This module performs:

* logical AND gating
* logical OR gating
* ensemble methylation analysis
* combination ranking
* MCC-based ensemble evaluation

---

# Workflow Logic

The workflow:

1. Parses ASCII logs
2. Extracts methylation vectors
3. Reconstructs desired binary states
4. Compares observed vs expected methylation
5. Computes confusion-matrix statistics
6. Generates plots and PDFs
7. Evaluates ensemble gate combinations

---

# Metrics Computed

| Metric          | Description                         |
| --------------- | ----------------------------------- |
| Accuracy        | Overall correct calls               |
| Sensitivity     | Correct methylated-site detection   |
| Specificity     | Correct unmethylated-site detection |
| MCC             | Balanced correlation metric         |
| False positives | `0 → 1` flips                       |
| False negatives | `1 → 0` flips                       |

---

# Generated Outputs

The workflow produces:

```text id="jlwm9y"
ahead_comparison_per_brick.csv
ahead_comparison_summary.csv
ahead_comparison_figures.pdf
```

---

# Summary CSV

Contains:

* total bricks
* matches
* error rates
* sensitivity
* specificity
* MCC
* decoded words

---

# Per-Brick CSV

Contains per-brick comparisons:

* desired state
* actual state
* correctness
* flip type

---

# Multi-Page PDF Outputs

The workflow generates:

```text id="jlwmx2"
ahead_comparison_figures.pdf
```

containing:

| Page | Plot                      |
| ---- | ------------------------- |
| 1    | Overall model performance |
| 2    | Per-sample accuracy       |
| 3    | Flip-type breakdown       |
| 4    | Per-brick heatmap         |
| 5    | MCC scatter plot          |

---

# Heatmap Visualization

Heatmaps encode:

| Color | Meaning                      |
| ----- | ---------------------------- |
| Green | Correct                      |
| Amber | Missed methylation (`1 → 0`) |
| Red   | False positive (`0 → 1`)     |
| White | Missing data                 |

CpG count tracks are plotted beneath heatmaps.

---

# Gate-Based Ensemble Analysis

The workflow additionally evaluates ensemble combinations using:

* AND gates
* OR gates

across:

* 2-model combinations
* 3-model combinations
* 4-model combinations

---

# Gate Logic

## AND Gate

A brick is methylated only if:

```text id="jlwmpr"
ALL models predict methylation
```

This improves:

* specificity
* false-positive reduction

but may reduce sensitivity.

---

## OR Gate

A brick is methylated if:

```text id="jlwm9z"
ANY model predicts methylation
```

This improves:

* sensitivity

but may increase false positives.

---

# Gate Analysis Outputs

Generated outputs include:

```text id="jlwmzt"
ahead_comparison_gate_OR_2way.pdf
ahead_comparison_gate_OR_3way.pdf
ahead_comparison_gate_OR_4way.pdf

ahead_comparison_gate_AND_2way.pdf
ahead_comparison_gate_AND_3way.pdf
ahead_comparison_gate_AND_4way.pdf
```

and summary CSVs:

```text id="jlwmw4"
ahead_comparison_gate_OR_summary.csv
ahead_comparison_gate_AND_summary.csv
```

---

# Top-Performing Combination Analysis

The workflow ranks:

* individual models
* AND combinations
* OR combinations

using:

```text id="jlwm3g"
mean MCC across samples
```

The top-performing combinations are visualized in:

```text id="jlwmxr"
ahead_comparison_top5_mcc_comparison.pdf
```

---

# Running the Workflow

## Step 1 — Move to Analysis Directory

```bash id="jlwm1x"
cd final_model_comparison
```

---

## Step 2 — Run Model Comparison

```bash id="jlwmu8"
python updated_AHEAD_model_comparison.py \
    . \
    ahead_comparison
```

Arguments:

| Argument        | Description                       |
| --------------- | --------------------------------- |
| First argument  | Directory containing `.log` files |
| Second argument | Output prefix                     |

---

# Verifying Outputs

List generated files:

```bash id="jlwm0y"
ls -lh ahead_comparison*
```

Inspect CSV summaries:

```bash id="jlwm9j"
head ahead_comparison_summary.csv
```

---

# Documentation

# Full Documentation

Detailed workflow documentation is available here:

[Google Docs Documentation](https://docs.google.com/document/d/1Wj-gkxO755uF2FdEJx0VJSViN6hjUYA6_AcrjSxTIXk/edit?tab=t.8rlwy8hh1c9h)
