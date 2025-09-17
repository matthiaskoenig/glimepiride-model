"""Module for plotting simulated CYP2C9 PK parameters for glimepiride and its metabolites."""

from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator
import matplotlib.transforms as transforms
from scipy.stats import gaussian_kde, ks_2samp
from statsmodels.stats.multitest import multipletests
from sbmlutils.console import console
from pkdb_models.models.glimepiride import DATA_PATH_GLIMEPIRIDE
import pkdb_models.models.glimepiride as glimepiride
from pkdb_models.models.glimepiride.experiments.base_experiment import GlimepirideSimulationExperiment
import itertools
import os
from pathlib import Path

_override = os.getenv("GLIMEPIRIDE_RESULTS_PATH")
if _override:
    p = Path(_override).resolve()
    p.mkdir(parents=True, exist_ok=True)
    glimepiride.RESULTS_PATH = p
    glimepiride.RESULTS_PATH_SIMULATION = p / "simulation"
    glimepiride.RESULTS_PATH_ANALYSES = p / "cyp2c9_analyses"
    glimepiride.RESULTS_PATH_FIT = p / "fit"


# Configuration
parameter_labels = {
    "cmax": "Cmax [µM]",
    "tmax": "Tmax [hr]",
    "auc": "AUC [ng*hr/ml]"
}

axis_limits = {
    "auc": {"glimepiride": (0, 5000), "M1": (0, 700), "M2": (0, 450)},
    "cmax": {"glimepiride": (0.2, 1.0), "M1": (0, 0.3), "M2": (0, 0.15)},
    "tmax": {"glimepiride": (1.25, 3.5), "M1": (2.2, 5.2), "M2": (3.1, 7.5)}
}

genotypes_of_interest = ["*1/*1", "*1/*2", "*1/*3", "*3/*3"]

biogeographical_groups_order = [
    "African American/Afro-Caribbean",
    "American",
    "Central/South Asian",
    "East Asian",
    "European",
    "Latino",
    "Near Eastern",
    "Oceanian",
    "Sub-Saharan African",
]

# Colors for biogeographical groups
biogeographical_colors = [
    '#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf',
    '#e0f3f8', '#abd9e9', '#74add1', '#4575b4'
]

# Figure sizes
figure_size_boxplot = (14, 6)
figure_size_ridgeline = (16, 8)
figure_size_volcano = (11, 7)

# Plot styling
boxplot_alpha = 0.9
boxplot_linewidth = 1.5
scatter_alpha = 0.3
scatter_size_small = 3
scatter_size_large = 7
ridgeline_overlap = 0.3
ridgeline_alpha = 0.5

# Statistical thresholds
min_samples_per_group = 5
fdr_alpha = 0.05
volcano_p_threshold = 0.01
volcano_percentile_threshold = 75

# Font sizes
font_size_title = 22
font_size_axis_label = 18
font_size_tick_label = 16
font_size_genotype = 15
font_size_legend = 14

# Base experiment for molecular weights
base_experiment = GlimepirideSimulationExperiment()

# Unit conversion info
conversion_info = {
    "tmax": {"factor": 1/60, "old_unit": "minute", "new_unit": "hour"},
    "cmax": {"factor": 1000, "old_unit": "millimole/liter", "new_unit": "µM"},
    "auc": {"old_unit": "millimole*minute/liter", "new_unit": "ng*hr/ml"}
}

# Random state for reproducibility
rng = np.random.RandomState(42)


def get_conversion_info(param: str, compound: str) -> Tuple[Optional[float], str, str]:
    """Get conversion factor and unit info for PK parameter."""
    if param == "auc":
        # Get molecular weight
        if compound == "glimepiride":
            mw = base_experiment.Mr.gli.magnitude
        elif compound == "M1":
            mw = base_experiment.Mr.m1.magnitude
        elif compound == "M2":
            mw = base_experiment.Mr.m2.magnitude
        else:
            mw = None

        # AUC conversion needs molecular weight
        factor = mw * 1000 / 60 if mw else None
        return factor, conversion_info[param]["old_unit"], conversion_info[param]["new_unit"]
    else:
        # Fixed conversion factors
        info = conversion_info[param]
        return info["factor"], info["old_unit"], info["new_unit"]


def convert_pk_dataframe(df: pd.DataFrame, params: List[str],
                        compounds: List[str]) -> pd.DataFrame:
    """Convert PK parameters in dataframe to appropriate units."""
    df_converted = df.copy()

    for param in params:
        converted_col = f"{param}_converted"
        df_converted[converted_col] = np.nan

        for compound in compounds:
            mask = df_converted["compound"] == compound
            if mask.any():
                factor, _, _ = get_conversion_info(param, compound)
                if factor is not None:
                    df_converted.loc[mask, converted_col] = (
                        df_converted.loc[mask, param] * factor
                    )

    return df_converted


def calculate_summary_statistics(df: pd.DataFrame, value_column: str,
                               group_column: str) -> pd.DataFrame:
    """Calculate summary statistics by group."""
    return df.groupby(group_column)[value_column].agg(
        mean="mean",
        median="median",
        min="min",
        Q1=lambda x: x.quantile(0.25),
        Q3=lambda x: x.quantile(0.75),
        max="max"
    )


def plot_study_data(ax, row: pd.Series, y_pos: float, compound: str, param: str,
                   alpha: float = 0.6, marker_size: int = 3) -> None:
    """Plot study data with error bars or single points."""
    # Add reproducible jitter
    y_jitter = 0.6 * (rng.random() - 0.5)
    y_pos_jittered = y_pos + y_jitter

    # Get values from study row
    val_single = row.get("value_adjusted_converted", np.nan)
    val_mean = row.get("mean_adjusted_converted", np.nan)
    val_sd = row.get("sd_adjusted_converted", np.nan)
    val_min = row.get("min_adjusted_converted", np.nan)
    val_max = row.get("max_adjusted_converted", np.nan)

    # Get conversion factor
    conv_factor, _, _ = get_conversion_info(param, compound)

    if pd.isna(conv_factor):
        return

    # Plot according to available data
    if not pd.isna(val_single):
        ax.plot(val_single * conv_factor, y_pos_jittered,
               marker="s", markersize=marker_size,
               color="black", lw=1.0, alpha=alpha)
    elif not pd.isna(val_mean):
        if not pd.isna(val_sd):
            # Mean ± SD
            ax.errorbar(val_mean * conv_factor, y_pos_jittered,
                       xerr=val_sd * conv_factor,
                       fmt="s", markersize=marker_size,
                       color="black", ecolor="black",
                       capsize=3, lw=1.0, alpha=alpha)
        elif not pd.isna(val_min) and not pd.isna(val_max):
            # Mean with min/max
            left_err = (val_mean - val_min) * conv_factor
            right_err = (val_max - val_mean) * conv_factor
            ax.errorbar(val_mean * conv_factor, y_pos_jittered,
                       xerr=[[left_err], [right_err]],
                       fmt="s", markersize=marker_size,
                       color="black", ecolor="black",
                       capsize=3, lw=1.0, alpha=alpha)
    elif not pd.isna(val_min) and not pd.isna(val_max):
        # Only min/max - plot midpoint
        mid = (val_min + val_max) / 2.0
        left_err = (mid - val_min) * conv_factor
        right_err = (val_max - mid) * conv_factor
        ax.errorbar(mid * conv_factor, y_pos_jittered,
                   xerr=[[left_err], [right_err]],
                   fmt="s", markersize=marker_size,
                   color="black", ecolor="black",
                   capsize=3, lw=1.0, alpha=alpha)


def genotype_boxplots(df: pd.DataFrame, pk_params: List[str], compounds: List[str],
                     output_dir: Path) -> None:
    """Create boxplots by genotype with overlaid study data."""
    sns.set_style("whitegrid", {'axes.grid': False})

    # Set up genotype positions and colors
    genotype_positions = {g: idx for idx, g in enumerate(genotypes_of_interest)}
    cyp2c9_colors = GlimepirideSimulationExperiment.cyp2c9_colors_analyses

    # Load experimental study data
    df_studies = pd.read_excel(
        DATA_PATH_GLIMEPIRIDE / "cyp2c9_data.xlsx",
        sheet_name="pk_params_studies",
        comment="#",
        engine="openpyxl"
    )

    # Filter simulation data
    df_filtered = df[df["genotype"].isin(genotypes_of_interest)].copy()

    # Print summary statistics
    for compound in compounds:
        for param in pk_params:
            df_comp = df_filtered[df_filtered["compound"] == compound].copy()
            if not df_comp.empty:
                factor, _, _ = get_conversion_info(param, compound)
                if factor is not None:
                    df_comp[param] = df_comp[param] * factor
                summary = calculate_summary_statistics(df_comp, param, "genotype")
                console.print(f"\n{compound} - {param} summary by genotype:\n{summary}\n")

    # Create figure for each parameter
    for param in pk_params:
        fig, axes = plt.subplots(1, len(compounds), figsize=figure_size_boxplot, sharey=False)

        # Process each compound
        for ax_idx, (ax, compound) in enumerate(zip(axes, compounds)):
            ax.grid(True, linestyle='--', alpha=0.3, color='lightgrey')

            # Filter and convert data
            df_comp = df_filtered[df_filtered["compound"] == compound].copy()
            if df_comp.empty:
                ax.set_visible(False)
                continue

            factor, _, _ = get_conversion_info(param, compound)
            if factor is not None:
                df_comp[param] = df_comp[param] * factor

            # Plot boxplots
            sns.boxplot(
                data=df_comp,
                x=param,
                y="genotype",
                ax=ax,
                order=genotypes_of_interest,
                medianprops={"color": "black"},
                whiskerprops={"linewidth": boxplot_linewidth},
                capprops={"linewidth": boxplot_linewidth},
                boxprops={"linewidth": boxplot_linewidth},
                showfliers=False
            )

            # Color boxes by genotype
            boxes = [p for p in ax.patches if isinstance(p, mpatches.PathPatch)]
            for box_idx, genotype in enumerate(genotypes_of_interest):
                if box_idx < len(boxes):
                    boxes[box_idx].set_facecolor(cyp2c9_colors[genotype])
                    boxes[box_idx].set_alpha(boxplot_alpha)

            # Add study data points
            df_studies_comp = df_studies[
                (df_studies["measurement_type"] == param) &
                (df_studies["Substance"] == compound) &
                (df_studies["Genotype"].isin(genotypes_of_interest))
            ]

            # Track values by genotype for aggregation
            genotype_study_values = {}

            for _, row_study in df_studies_comp.iterrows():
                genotype = row_study["Genotype"]
                y_pos = genotype_positions[genotype]

                # Plot individual points
                plot_study_data(ax, row_study, y_pos, compound, param,
                               alpha=scatter_alpha, marker_size=scatter_size_small)

                # Collect values for aggregation
                val = None
                if pd.notna(row_study.get("value_adjusted_converted")):
                    val = row_study.get("value_adjusted_converted")
                elif pd.notna(row_study.get("mean_adjusted_converted")):
                    val = row_study.get("mean_adjusted_converted")
                elif pd.notna(row_study.get("min_adjusted_converted")) and pd.notna(row_study.get("max_adjusted_converted")):
                    val = (row_study.get("min_adjusted_converted") + row_study.get("max_adjusted_converted")) / 2.0

                if val is not None:
                    if genotype not in genotype_study_values:
                        genotype_study_values[genotype] = []
                    genotype_study_values[genotype].append(val)

            # Plot aggregated study points
            conv_factor, _, _ = get_conversion_info(param, compound)
            for genotype, values in genotype_study_values.items():
                if len(values) > 0:
                    y_pos = genotype_positions[genotype] + 0.1
                    mean_val = np.mean(values)

                    if len(values) > 1:
                        sd_val = np.std(values, ddof=1)
                        ax.errorbar(
                            mean_val * conv_factor, y_pos,
                            xerr=sd_val * conv_factor,
                            fmt='s', markersize=scatter_size_large,
                            markerfacecolor="black", markeredgecolor='black',
                            ecolor="black", capsize=5
                        )
                    else:
                        ax.plot(
                            mean_val * conv_factor, y_pos,
                            marker='s', markersize=scatter_size_large,
                            markerfacecolor="black", markeredgecolor='black'
                        )

            # Set axis labels
            ax.set_yticks(range(len(genotypes_of_interest)))
            if ax_idx == 0:
                ax.set_yticklabels(genotypes_of_interest, fontsize=font_size_genotype,
                                  fontweight='bold')
                ax.set_ylabel("Genotype", fontsize=font_size_axis_label, fontweight='bold')
            else:
                ax.set_yticklabels([])
                ax.set_ylabel("")

            ax.set_xlabel(parameter_labels[param], fontsize=font_size_tick_label, fontweight='bold')
            if param == "cmax":
                ax.set_title(compound.title(), fontsize=font_size_title, fontweight='bold')
            ax.tick_params(axis="x", labelrotation=45)

        # Set title and save
        plt.suptitle(f"Simulated vs. Observed CYP2C9 pharmacokinetic parameter - Glimepiride",
                    fontsize=font_size_title, fontweight='bold')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        fig.savefig(output_dir / f"simulated_{param}_by_genotype_boxplots.png", bbox_inches="tight")
        plt.show()


def ridgeline(ax, data: pd.DataFrame, value_col: str, group_col: str = "biogeographical_group",
              color_map: Optional[Dict] = None, label_left: bool = True,
              overlap: float = 0.6, category_order: Optional[List] = None) -> None:
    """Draw KDE ridgelines on an axis."""
    categories = category_order if category_order else sorted(data[group_col].dropna().unique())
    n = len(categories)
    cat_offsets = list(enumerate(reversed(categories)))

    # Process each category from bottom to top
    for k, cat in cat_offsets:
        subset = data.loc[data[group_col] == cat, value_col].dropna()
        if subset.empty:
            continue

        y_offset = k * overlap
        z_order = n - k

        # Calculate kernel density estimation
        subset_min, subset_max = subset.min(), subset.max()
        range_width = subset_max - subset_min
        grid_min = max(0, subset_min - 0.5 * range_width)
        grid_max = subset_max + 0.5 * range_width
        x_grid = np.linspace(grid_min, grid_max, 300)
        kde_vals = gaussian_kde(subset)(x_grid)
        kde_vals = (kde_vals / kde_vals.max()) * 0.8  # Normalize height

        # Plot filled area under KDE curve
        face_color = color_map[cat] if (color_map and cat in color_map) else "blue"
        ax.fill_between(
            x_grid, y_offset, y_offset + kde_vals,
            facecolor=face_color, edgecolor="none",
            alpha=ridgeline_alpha, zorder=z_order
        )
        ax.plot(x_grid, y_offset + kde_vals, color="black", linewidth=1, zorder=z_order)

        # Add group label
        if label_left:
            trans = transforms.blended_transform_factory(ax.transAxes, ax.transData)
            ax.text(-0.05, y_offset + 0.4, cat, va="center", ha="right",
                    fontsize=font_size_genotype, fontweight='bold',
                    transform=trans, zorder=z_order)

    ax.set_yticks([])
    ax.set_ylim(0, (n - 1) * overlap + 1.0)

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.5)


def pk_ridgelines_biogeographical_groups(df: pd.DataFrame, pk_params: List[str],
                                        compounds: List[str], output_dir: Path) -> None:
    """Generate ridgeline plots by biogeographical groups."""
    sns.set_style("whitegrid")

    # Create color mapping
    group_to_color = dict(zip(biogeographical_groups_order, biogeographical_colors))

    # Convert data to appropriate units
    df_converted = convert_pk_dataframe(df, pk_params, compounds)

    # Rename column for display
    df_converted["biogeographical_group"] = df_converted["ethnicity"]

    # Create plots for each parameter
    for param in pk_params:
        fig, axes = plt.subplots(1, len(compounds), figsize=figure_size_ridgeline, sharey=False)
        fig.subplots_adjust(left=0.2, top=0.85)

        for ax_idx, (ax, compound) in enumerate(zip(axes, compounds)):
            ax.grid(True, linestyle='--', alpha=0.3, color='lightgrey')

            # Filter data for current compound
            df_comp = df_converted[df_converted["compound"] == compound].copy()
            if df_comp.empty:
                ax.set_visible(False)
                continue

            value_col = f"{param}_converted"

            # Set axis limits and labels
            x_min, x_max = axis_limits[param][compound]
            ax.set_xlim(x_min, x_max)
            ax.set_xlabel(parameter_labels[param], fontsize=font_size_axis_label, fontweight='bold')
            ax.tick_params(axis='x', labelsize=font_size_tick_label)

            # Log KDE statistics
            for group in biogeographical_groups_order:
                subset = df_comp[df_comp["biogeographical_group"] == group][value_col].dropna()
                if len(subset) > 1:
                    kde = gaussian_kde(subset)
                    x_grid = np.linspace(subset.min(), subset.max(), 200)
                    y_grid = kde(x_grid)
                    mode_val = x_grid[np.argmax(y_grid)]

                    console.print(f"{compound} - {param} - {group}:", style="bold cyan")
                    console.print(
                        f"  Mode: {mode_val:.4f} | Mean: {np.mean(subset):.4f} | "
                        f"Median: {np.median(subset):.4f} | SD: {np.std(subset):.4f}"
                    )

            # Create ridgeline plot
            ridgeline(
                ax=ax,
                data=df_comp,
                value_col=value_col,
                group_col="biogeographical_group",
                color_map=group_to_color,
                label_left=(ax_idx == 0),
                overlap=ridgeline_overlap,
                category_order=biogeographical_groups_order,
            )

            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.set_title(compound.title(), fontsize=font_size_title, fontweight='bold')

        # Save plot
        plt.tight_layout(rect=[0, 0, 1, 0.93])
        fig.savefig(output_dir / f"ridgeline_{param}_by_biogeographical_groups.png",
                   dpi=600, bbox_inches="tight")
        plt.show()


def calculate_biogeographical_pk_differences(df_data: pd.DataFrame, pk_parameters_to_test: List[str],
                                            compounds_to_test: List[str],
                                            groups_to_test: List[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate differences in PK parameters between biogeographical groups."""
    console.print("Calculating biogeographical differences in PK parameters...")

    # Convert data
    df_converted = convert_pk_dataframe(df_data, pk_parameters_to_test, compounds_to_test)

    all_comparisons = []

    # Perform all pairwise comparisons
    for compound in compounds_to_test:
        for param in pk_parameters_to_test:
            value_col = f"{param}_converted"
            for group1, group2 in itertools.combinations(groups_to_test, 2):
                # Extract data for each group
                data1 = df_converted[
                    (df_converted["compound"].str.lower() == compound.lower()) &
                    (df_converted["ethnicity"] == group1)
                ][value_col].dropna()

                data2 = df_converted[
                    (df_converted["compound"].str.lower() == compound.lower()) &
                    (df_converted["ethnicity"] == group2)
                ][value_col].dropna()

                # Only compare if both groups have minimum sample size
                if len(data1) >= min_samples_per_group and len(data2) >= min_samples_per_group:
                    # Kolmogorov-Smirnov test
                    ks_stat, p_value = ks_2samp(data1, data2)
                    median1 = data1.median()
                    median2 = data2.median()

                    # Calculate effect sizes
                    perc_diff = np.nan
                    log2_fc = np.nan

                    if pd.notna(median1) and pd.notna(median2) and median1 != 0:
                        perc_diff = ((median2 - median1) / median1) * 100
                        if median2 > 0 and median1 > 0:
                            log2_fc = np.log2(median2 / median1)

                    # Simplify group names for display
                    group1_display = group1.split('/')[0] if '/' in group1 else group1
                    group2_display = group2.split('/')[0] if '/' in group2 else group2

                    all_comparisons.append({
                        "Compound": compound.title(), "Parameter": param.upper(),
                        "Group 1": group1_display, "Group 2": group2_display,
                        "KS Stat": ks_stat, "Raw P-value": p_value,
                        "% Diff": perc_diff, "Log2FC": log2_fc,
                        "_group1_full": group1, "_group2_full": group2,
                        "_median1": median1, "_median2": median2
                    })

    stats_df_all_pairs = pd.DataFrame(all_comparisons)

    # Apply multiple testing correction
    if stats_df_all_pairs.empty or stats_df_all_pairs["Raw P-value"].isna().all():
        console.print("No valid raw p-values to adjust. Skipping FDR.")
        stats_df_all_pairs["Adj. p-value"] = np.nan
        stats_df_all_pairs["-log10 Adj. p-value"] = np.nan
        stats_df_all_pairs["Sig."] = ""
    else:
        _, pvals_corrected, _, _ = multipletests(
            stats_df_all_pairs["Raw P-value"].fillna(1.0),
            method='fdr_bh',
            alpha=fdr_alpha
        )
        stats_df_all_pairs["Adj. p-value"] = pvals_corrected

        # Calculate -log10 p-values
        epsilon = 1e-300
        stats_df_all_pairs['-log10 Adj. p-value'] = -np.log10(
            stats_df_all_pairs['Adj. p-value'].replace(0, epsilon).astype(float)
        )

        # Mark significance levels
        stats_df_all_pairs["Sig."] = ""
        stats_df_all_pairs.loc[stats_df_all_pairs["Adj. p-value"] < 0.001, "Sig."] = "***"
        stats_df_all_pairs.loc[(stats_df_all_pairs["Adj. p-value"] >= 0.001) &
                              (stats_df_all_pairs["Adj. p-value"] < 0.01), "Sig."] = "**"
        stats_df_all_pairs.loc[(stats_df_all_pairs["Adj. p-value"] >= 0.01) &
                              (stats_df_all_pairs["Adj. p-value"] < 0.05), "Sig."] = "*"

    # Extract highly significant results
    if "Sig." not in stats_df_all_pairs.columns:
        significant_results_table_df = pd.DataFrame()
    else:
        significant_results_table_df = stats_df_all_pairs[
            stats_df_all_pairs["Sig."].isin(["***", "**"])
        ].copy()

    final_columns_table = ["Compound", "Parameter", "Group 1", "Group 2",
                          "KS Stat", "Adj. p-value", "Sig.", "% Diff"]

    if not significant_results_table_df.empty:
        # Format table for display
        significant_results_table_df = significant_results_table_df.sort_values(
            by=["Compound", "Parameter", "Adj. p-value"]
        )

        # Format numeric columns
        table_format_cols = {
            "KS Stat": '{:.2f}',
            "Adj. p-value": '{:.2e}',
            "% Diff": '{:.1f}%'
        }
        for col, fmt_str in table_format_cols.items():
            if col in significant_results_table_df:
                significant_results_table_df[col] = significant_results_table_df[col].map(
                    lambda x: fmt_str.format(x) if pd.notna(x) else ''
                )

        significant_results_table_df = significant_results_table_df[final_columns_table]

        console.print("Most Significant Biogeographical Differences:")
        console.print(significant_results_table_df.to_string(index=False))
    else:
        console.print("No significant differences found (Adj. p < 0.01).")
        significant_results_table_df = pd.DataFrame(columns=final_columns_table)

    return stats_df_all_pairs, significant_results_table_df


def generate_volcano_plot(all_pairs_df: pd.DataFrame, compound: str, param: str,
                         output_dir: Path) -> None:
    """Generate volcano plot for specific compound and parameter."""
    if all_pairs_df.empty:
        console.print("[Warning] No data available for volcano plot generation.")
        return

    # Filter data
    volcano_data = all_pairs_df[
        (all_pairs_df['Compound'] == compound.title()) &
        (all_pairs_df['Parameter'] == param.upper())
    ].copy()

    if volcano_data.empty:
        console.print(f"No data found for {compound} {param} to generate volcano plot.")
        return

    # Prepare data for plotting
    plot_data = volcano_data.copy()
    cols_to_check = ['Adj. p-value', '% Diff', '_group1_full', '_group2_full']

    for col in cols_to_check:
        if col in plot_data.columns:
            if col not in ['_group1_full', '_group2_full']:
                plot_data[col] = pd.to_numeric(plot_data[col], errors='coerce')
        else:
            console.print(f"[Warning] Required column '{col}' missing for volcano plot.")
            plot_data[col] = np.nan

    # Remove invalid data
    plot_data.dropna(subset=['Adj. p-value', '% Diff', '_group1_full', '_group2_full'], inplace=True)
    plot_data = plot_data[plot_data['Adj. p-value'] > 0]

    if plot_data.empty:
        console.print(f"Not enough valid data for volcano plot ({compound} {param}).")
        return

    plot_data['Abs % Diff'] = abs(plot_data['% Diff'])

    # Create figure with significance thresholds
    plt.figure(figsize=figure_size_volcano)

    # Quantile-based threshold
    abs_perc_diff_thresh = plot_data['Abs % Diff'].quantile(volcano_percentile_threshold / 100)
    console.print(f"Using {volcano_percentile_threshold}th percentile as threshold: {abs_perc_diff_thresh:.2f}%")
    p_value_thresh = volcano_p_threshold

    # Plot non-significant points
    non_sig_points = plot_data[
        (plot_data['Abs % Diff'] < abs_perc_diff_thresh) |
        (plot_data['Adj. p-value'] > p_value_thresh)
    ]
    plt.scatter(non_sig_points['Abs % Diff'], non_sig_points['Adj. p-value'],
                alpha=0.9, s=30, color='black')

    # Plot and label significant points
    sig_points = plot_data[
        (plot_data['Abs % Diff'] >= abs_perc_diff_thresh) &
        (plot_data['Adj. p-value'] <= p_value_thresh)
    ]

    if not sig_points.empty:
        # Use square markers with tab:blue fill
        plt.scatter(sig_points['Abs % Diff'], sig_points['Adj. p-value'],
                    alpha=1, s=50, marker='s', facecolor='tab:blue',
                    edgecolor='black', linewidth=1)

        # Add labels for significant points
        texts = []
        for _, row in sig_points.iterrows():
            group1_full = str(row['_group1_full'])
            group2_full = str(row['_group2_full'])
            label_text = f"{group1_full} & {group2_full}"

            # Default label positioning
            x_val = row['Abs % Diff']
            y_val = row['Adj. p-value']
            ha_val = 'left'
            va_val = 'bottom'
            fontsize = 14
            color = 'black'

            # Custom positioning for specific pairs
            if ((group1_full == "Latino" and group2_full == "Oceanian") or
                (group1_full == "Oceanian" and group2_full == "Latino") or
                (group1_full == "Central/South Asian" and group2_full == "Oceanian")):
                ha_val = 'right'
                va_val = 'bottom'
            elif group1_full == "African American/Afro-Caribbean":
                ha_val = 'right'

            texts.append(plt.text(x_val, y_val, label_text,
                                fontsize=fontsize, fontweight='normal',
                                color=color, ha=ha_val, va=va_val))

    # Add threshold lines
    plt.axvline(x=abs_perc_diff_thresh, color='black', linestyle='--', lw=0.7)
    plt.axhline(y=p_value_thresh, color='black', linestyle='--', lw=0.7)

    # Format plot with log scale and inverted y-axis
    plt.yscale('log')
    plt.gca().invert_yaxis()

    # Set axis limits and ticks
    plt.xlim(-0.5, 13)
    plt.xticks([0, 2, 4, 6, 8, 10, 12])

    # Fix y-axis labels to avoid glyph error
    plt.yticks([1e0, 1e-1, 1e-2], ['1', '0.1', '0.01'])

    # Axis labels
    plt.xlabel(f"Absolute Difference {compound} {param} [%]",
              fontsize=font_size_axis_label, fontweight='bold')
    plt.ylabel("p-value", fontsize=font_size_axis_label, fontweight='bold')

    # Grid
    plt.grid(True, which="both", linestyle=':', alpha=0.4)

    plt.tight_layout(rect=[0, 0.05, 1, 0.95])

    # Save plot
    vp_path = output_dir / f"volcano_plot_{compound}_{param}.png"
    plt.savefig(vp_path, dpi=600)

    plt.show()
    plt.close()


if __name__ == "__main__":
    # Define parameters and compounds to analyze
    pk_params = ["auc", "tmax", "cmax"]
    compounds = ["glimepiride", "M1", "M2"]

    # Set up output directory
    output_figure_dir = glimepiride.RESULTS_PATH_ANALYSES / "_figures"
    output_figure_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    tsv_path_all = Path(__file__).parent / "GlimepirideCYP2C9Scan_cyp2c9_pharmacokinetics.tsv"
    df_all = pd.read_csv(tsv_path_all, sep="\t")

    tsv_path_bio = Path(__file__).parent / "GlimepirideCYP2C9Scan_ethnicities_pharmacokinetics.tsv"
    df_bio = pd.read_csv(tsv_path_bio, sep="\t")

    # ===== GENOTYPE-BASED ANALYSIS =====
    genotype_boxplots(df_all, pk_params, compounds, output_dir=output_figure_dir)

    # ===== BIOGEOGRAPHICAL GROUP-BASED ANALYSIS =====
    pk_ridgelines_biogeographical_groups(df_bio, pk_params, compounds, output_dir=output_figure_dir)

    # ===== STATISTICAL COMPARISONS =====
    all_pairs_df, significant_table_df = calculate_biogeographical_pk_differences(
        df_data=df_bio,
        pk_parameters_to_test=pk_params,
        compounds_to_test=compounds,
        groups_to_test=biogeographical_groups_order
    )

    # ===== SAVE RESULTS =====
    # Save all pairwise comparisons
    # supplemental_data_path = output_figure_dir / "all_biogeographical_pairwise_comparisons_supplement.csv"
    cols_for_supplement = ["Compound", "Parameter", "Group 1", "Group 2",
                          "_group1_full", "_group2_full", "_median1", "_median2",
                          "KS Stat", "Raw P-value", "Adj. p-value",
                          "Sig.", "% Diff", "Log2FC", "-log10 Adj. p-value"]
    actual_cols = [col for col in cols_for_supplement if col in all_pairs_df.columns]

    # all_pairs_df[actual_cols].to_csv(supplemental_data_path, index=False, float_format='%.4g')

    # Save significant differences table
    # table_path = output_figure_dir / "significant_biogeographical_pk_differences.csv"
    # significant_table_df.to_csv(table_path, index=False)

    # ===== VOLCANO PLOT =====
    generate_volcano_plot(all_pairs_df, "Glimepiride", "AUC", output_figure_dir)