"""Module for sampling CYP2C9 allele and diplotype activities."""

from pathlib import Path
import numpy as np
import scipy.stats as stats
from scipy.stats import gaussian_kde
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from sbmlutils.console import console
from pkdb_models.models.glimepiride import RESULTS_PATH
from pkdb_models.models.glimepiride.experiments.base_experiment import GlimepirideSimulationExperiment
from pkdb_models.models.glimepiride import DATA_PATH_GLIMEPIRIDE
import pkdb_models.models.glimepiride as glimepiride
from pkdb_models.models.glimepiride.experiments.cyp2c9.analyze_cyp2c9_activity import get_lognorm_pars
import os

_override = os.getenv("GLIMEPIRIDE_RESULTS_PATH")
if _override:
    p = Path(_override).resolve()
    p.mkdir(parents=True, exist_ok=True)
    glimepiride.RESULTS_PATH = p
    glimepiride.RESULTS_PATH_SIMULATION = p / "simulation"
    glimepiride.RESULTS_PATH_ANALYSES = p / "cyp2c9_analyses"
    glimepiride.RESULTS_PATH_FIT = p / "fit"

np.random.seed(111)
cyp2c9_allele_activity = GlimepirideSimulationExperiment.cyp2c9_allele_activity
cyp2c9_colors = GlimepirideSimulationExperiment.cyp2c9_colors_analyses
populations = [
    "African American/Afro-Caribbean", "American", "Central/South Asian",
    "East Asian", "European", "Latino", "Near Eastern", "Oceanian", "Sub-Saharan African"
]
genotypes_of_interest = ["*1/*1", "*1/*2", "*1/*3", "*3/*3"]

def calculate_kde_mode(data: np.ndarray, grid_min=0, grid_max=4, num_points=200):
    """Calculate distribution mode using KDE."""
    if len(data) <= 1:
        return None, None, None, None
    x_grid = np.linspace(grid_min, grid_max, num_points)
    kde = gaussian_kde(data)
    y_grid = kde(x_grid)
    mode_idx = np.argmax(y_grid)
    return x_grid[mode_idx], kde, x_grid, y_grid


def setup_ethnicity_plots(sharey=True, sharex=True):
    """Setup figure and axes for ethnicity plots."""
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, axes = plt.subplots(3, 3, figsize=(16, 12), sharey=sharey, sharex=sharex)
    axes = axes.flatten()

    for k, ax in enumerate(axes.flat):
        ax.grid(True, linestyle='--', alpha=0.5, color='lightgrey')
        row, col = k // 3, k % 3

        if row >= 2:
            ax.set_xlabel("Activity", fontsize=16, fontweight='bold')

        if col == 0:
            ax.set_ylabel("Density", fontsize=16, fontweight='bold')

        ax.set_xlim(0, 3.0)
        ax.tick_params(labelbottom=(row == 2))

    return fig, axes


def filter_diplotypes(df_diplotype):
    """Filter diplotypes based on available allele activities."""
    def filter_diplotype(diplotype_str):
        a1, a2 = diplotype_str.split("/")
        return (a1 in cyp2c9_allele_activity) and (a2 in cyp2c9_allele_activity)

    mask_removed = ~df_diplotype["CYP2C9_allele"].apply(filter_diplotype)
    return df_diplotype[~mask_removed].copy(), df_diplotype[mask_removed].copy()


def format_scientific(value, precision=4):
    """Format numbers based on magnitude."""
    if pd.isna(value) or value is None or value == '':
        return "0"

    value = float(value)

    if value == 0:
        return "0"

    if abs(value) < 10 ** (-precision):
        return f"{value:.2e}"

    return f"{value:.{precision}f}".rstrip('0').rstrip('.')


def save_excluded_genotypes(df_diplotype_excluded, populations):
    """Save summary of excluded genotypes."""
    excluded_df = df_diplotype_excluded.sort_values(by="mean_freq", ascending=False)

    # Calculate totals
    population_totals = {pop: excluded_df[pop].sum() for pop in populations}
    total_excluded_freq = excluded_df["mean_freq"].sum()

    # Create summary row
    summary_data = {"CYP2C9_allele": "TOTAL", "mean_freq": total_excluded_freq}
    summary_data.update({pop: population_totals[pop] for pop in populations})

    # Add summary row
    summary_df = pd.DataFrame([summary_data])
    excluded_df_with_summary = pd.concat([excluded_df, summary_df], ignore_index=True)

    # Format values
    for col in populations + ["mean_freq"]:
        excluded_df_with_summary[col] = excluded_df_with_summary[col].apply(format_scientific)

    # Save to CSV
    excluded_output_path = Path(__file__).parent / "excluded_genotypes_summary.csv"
    excluded_df_with_summary.to_csv(excluded_output_path, index=False)

    # Log summary
    console.print(f"Total excluded frequency: {total_excluded_freq:.4f}")
    console.print(f"Top excluded diplotypes:")
    for idx, row in excluded_df.head(5).iterrows():
        console.print(f"  • {row['CYP2C9_allele']}: {row['mean_freq']:.4f}")

    console.rule(style="red")


def read_and_filter_excel(excel_path: Path, export_csv: bool = True):
    df_activity = pd.read_excel(excel_path, sheet_name="activity_studies", comment="#")
    df_allele = pd.read_excel(excel_path, sheet_name="allele_frequencies", comment="#")
    df_diplotype_original = pd.read_excel(excel_path, sheet_name="diplotype_frequencies", comment="#")

    df_filled = df_diplotype_original[populations].fillna(0.0)
    row_sums = df_filled.sum(axis=1)
    df_diplotype_nonempty = df_diplotype_original[row_sums != 0].copy()
    df_diplotype_nonempty["mean_freq"] = df_diplotype_nonempty[populations].fillna(0).mean(axis=1)
    df_diplotype_nonempty = df_diplotype_nonempty.sort_values(by="mean_freq", ascending=False)

    # optionally export a reproducible snapshot to results
    if export_csv:
        outdir = glimepiride.RESULTS_PATH_ANALYSES
        outdir.mkdir(parents=True, exist_ok=True)
        df_diplotype_nonempty.to_csv(outdir / "diplotype_nonempty.csv", index=False)

    df_allele_filtered = df_allele[df_allele["CYP2C9_allele"].isin(cyp2c9_allele_activity.keys())].copy()
    df_diplotype_filtered, df_diplotype_excluded = filter_diplotypes(df_diplotype_nonempty)

    console.rule("[bold red]EXCLUDED DIPLOTYPES[/bold red]", style="red")
    console.print(
        f"Original: {df_diplotype_nonempty.shape[0]} | Included: {df_diplotype_filtered.shape[0]} | Excluded: {df_diplotype_excluded.shape[0]}"
    )

    if not df_diplotype_excluded.empty:
        save_excluded_genotypes(df_diplotype_excluded, populations)

    return df_activity, df_allele_filtered, df_diplotype_filtered



def sample_allele(allele: str, shape: float, n=100000) -> np.ndarray:
    """Sample from allele distribution with mean-based scaling."""
    mean_val = cyp2c9_allele_activity[allele]
    scale = mean_val / np.exp(shape**2 / 2)
    return stats.lognorm.rvs(s=shape, scale=scale, size=n)


def sample_diplotype(allele1: str, allele2: str, shape: float, n=100000) -> np.ndarray:
    """Sample genotype activity by averaging allele samples."""
    s1 = sample_allele(allele1, shape, n=n)
    s2 = sample_allele(allele2, shape, n=n)
    return (s1 + s2) / 2


def collect_population_samples(df_dipl_f: pd.DataFrame, shape: float, sample_diplotype_func, total_samples: int, ethnicity: str):
    """Collect activity samples based on genotype frequencies."""
    population_data = []
    for _, row in df_dipl_f.iterrows():
        genotype = row["CYP2C9_allele"]
        freq = row.get(ethnicity, 0.0)
        if pd.isna(freq) or freq == 0 or total_samples * freq < 1:
            continue

        count = int(total_samples * freq)
        a1, a2 = genotype.split("/")
        genotype_samples = sample_diplotype_func(a1, a2, shape, n=count)

        for val in genotype_samples:
            population_data.append({"activity": val, "genotype": genotype})

    return pd.DataFrame(population_data)


def analyze_sampled_genotypes(df_dipl_f: pd.DataFrame, shape: float, sample_diplotype_func,
                              total_samples):
    """Analyze actual genotypes sampled during simulation by ethnicity
    (may differ from theoretical frequencies due to random sampling and
    rounding effects, especially for rare genotypes).
    """
    console.rule("[bold yellow]SAMPLED GENOTYPES BY ETHNICITY[/bold yellow]", style="yellow")
    sampling_counts = []

    for ethnicity in populations:
        # Get population samples
        pop_df = collect_population_samples(df_dipl_f, shape, sample_diplotype_func, total_samples, ethnicity)
        if pop_df.empty:
            continue

        # Count genotypes
        genotype_counts = pop_df["genotype"].value_counts().reset_index()
        genotype_counts.columns = ["genotype", "count"]
        genotype_counts["frequency"] = genotype_counts["count"] / total_samples
        genotype_counts["ethnicity"] = ethnicity
        sampling_counts.append(genotype_counts)

        # Print summary
        console.print(f"{ethnicity}: {len(genotype_counts)} genotypes sampled")

    # Create summary table
    if sampling_counts:
        sampled_df = pd.concat(sampling_counts, ignore_index=True)
        pivot_df = sampled_df.pivot_table(
            index="genotype",
            columns="ethnicity",
            values="frequency",
            fill_value=0
        ).reset_index()

        # Add mean frequency
        pivot_df["mean_freq"] = pivot_df.drop("genotype", axis=1).mean(axis=1)
        pivot_df = pivot_df.sort_values("mean_freq", ascending=False)

        # Calculate totals for each ethnicity
        population_totals = {pop: pivot_df[pop].sum() for pop in populations if pop in pivot_df.columns}
        total_mean_freq = pivot_df["mean_freq"].mean()

        # Create summary row
        summary_data = {"genotype": "TOTAL", "mean_freq": total_mean_freq}
        summary_data.update({pop: population_totals.get(pop, 0) for pop in populations})

        # Add summary row to pivot_df
        summary_df = pd.DataFrame([summary_data])
        pivot_df_with_summary = pd.concat([pivot_df, summary_df], ignore_index=True)

        # Apply format_scientific to all numeric columns
        formatted_df = pivot_df_with_summary.copy()
        for col in populations + ["mean_freq"]:
            if col in formatted_df.columns:
                formatted_df[col] = formatted_df[col].apply(format_scientific)

        # Save to CSV
        pivot_path = Path(__file__).parent / "sampled_genotypes_by_ethnicity.csv"
        formatted_df.to_csv(pivot_path, index=False)


def stacked_hist_ethnicites(
        df_dipl_f: pd.DataFrame,
        shape: float,
        sample_diplotype_func,
        total_samples,
        genotypes_of_interest=None,
        fig_path: Path = None
):
    """Plot stacked histograms highlighting key genotypes."""
    color_dict = cyp2c9_colors.copy()
    color_dict["other"] = "darkgrey"

    fig, axes = setup_ethnicity_plots()

    console.rule("[bold blue]ETHNICITY ACTIVITY STATISTICS[/bold blue]", style="blue")
    console.print(f"{'Ethnicity':<30} {'Mode':<10} {'Mean':<10} {'Median':<10} {'SD':<10}")
    console.rule(style="blue")

    for k, ethnicity in enumerate(populations):
        pop_df = collect_population_samples(df_dipl_f, shape, sample_diplotype_func, total_samples, ethnicity)
        if pop_df.empty:
            continue

        pop_df["genotype_interest"] = pop_df["genotype"].apply(
            lambda g: g if g in genotypes_of_interest else "other"
        )

        # Calculate statistics
        x_data = pop_df["activity"].dropna().values
        mode_val, _, _, _ = calculate_kde_mode(x_data)
        mean_val = np.mean(x_data)
        median_val = np.median(x_data)
        sd_val = np.std(x_data)

        if mode_val is not None:
            console.print(f"{ethnicity:<30} {mode_val:.3f}     {mean_val:.3f}     {median_val:.3f}     {sd_val:.3f}")

        # Plot histogram and KDE
        sns.histplot(
            data=pop_df,
            x="activity",
            hue="genotype_interest",
            multiple="stack",
            hue_order=["other", "*3/*3", "*1/*3", "*1/*2", "*1/*1"],
            palette=color_dict,
            alpha=0.7,
            bins=30,
            binrange=(0, 4),
            stat="density",
            legend=False,
            ax=axes[k]
        )

        sns.kdeplot(
            data=pop_df,
            x="activity",
            color="black",
            linewidth=2,
            ax=axes[k]
        )

        if mode_val is not None:
            axes[k].axvline(mode_val, color="black", linestyle="--", label="KDE mode")

        axes[k].set_title(ethnicity, fontsize=14, fontweight='bold')

    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor=color_dict["*1/*1"], edgecolor='black', label="*1/*1"),
        mpatches.Patch(facecolor=color_dict["*1/*2"], edgecolor='black', label="*1/*2"),
        mpatches.Patch(facecolor=color_dict["*1/*3"], edgecolor='black', label="*1/*3"),
        mpatches.Patch(facecolor=color_dict["*3/*3"], edgecolor='black', label="*3/*3"),
        mpatches.Patch(facecolor=color_dict["other"], edgecolor='black', label="Other genotypes"),
        plt.Line2D([0], [0], color="black", linestyle="--", label="KDE mode")
    ]

    fig.legend(
        handles=legend_elements,
        loc='lower center',
        ncol=len(legend_elements),
        bbox_to_anchor=(0.5, 0.02),
        frameon=True,
        fontsize=12
    )

    fig.suptitle("CYP2C9 Key Genotypes - Activity Distribution by Ethnicity",
                 fontsize=22, fontweight='bold')

    plt.tight_layout(rect=[0, 0.08, 1, 0.95])
    if fig_path is not None:
        fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close()


def boxplots_ethnicities(
        df_dipl_f: pd.DataFrame,
        shape: float,
        sample_diplotype_func,
        total_samples,
        fig_path: Path = None
):
    """Create boxplots of CYP2C9 activity by ethnicity."""
    plt.style.use('seaborn-v0_8-whitegrid')

    # Collect samples
    boxplot_data = []
    for ethnicity in populations:
        pop_df = collect_population_samples(df_dipl_f, shape, sample_diplotype_func, total_samples, ethnicity)
        if not pop_df.empty:
            for _, row in pop_df.iterrows():
                boxplot_data.append({"activity": row["activity"], "ethnicity": ethnicity})

    df_box = pd.DataFrame(boxplot_data)

    # Log statistics
    console.rule("[bold green]ETHNICITY BOXPLOT STATISTICS[/bold green]", style="green")
    console.print(f"{'Ethnicity':<30} {'Mean':<10} {'Median':<10} {'Q1':<10} {'Q3':<10}")
    console.rule(style="green")

    for ethnicity, group in df_box.groupby("ethnicity"):
        mean = group["activity"].mean()
        median = group["activity"].median()
        q1 = group["activity"].quantile(0.25)
        q3 = group["activity"].quantile(0.75)

        console.print(f"{ethnicity:<30} {mean:.3f}     {median:.3f}     {q1:.3f}     {q3:.3f}")

    # Create boxplot
    plt.figure(figsize=(10, 6))
    palette = ['#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4']
    ax = sns.boxplot(
        x="activity",
        y="ethnicity",
        hue="ethnicity",
        data=df_box,
        showfliers=False,
        palette=palette,
        legend=False
    )

    for patch in ax.patches:
        if isinstance(patch, mpatches.PathPatch):
            patch.set_alpha(0.8)

    ax.grid(True, linestyle='--', alpha=0.5, color='lightgrey')
    ax.set_title("CYP2C9 Activity by Ethnicity", fontsize=22, fontweight='bold')
    ax.set_xlabel("Activity", fontsize=16, fontweight='bold')
    ax.set_ylabel("", fontsize=16, fontweight='bold')
    ax.tick_params(axis='y', labelsize=16)

    for ytick in ax.get_yticklabels():
        ytick.set_fontweight('bold')

    plt.tight_layout()
    if fig_path is not None:
        plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close()


def stacked_hist_ethnicities_individual(
        df_dipl_f: pd.DataFrame,
        shape: float,
        sample_diplotype_func,
        total_samples,
        genotypes_of_interest=None,
        output_dir: Path = None
):
    """Plot stacked histograms separately for world map."""
    color_dict = cyp2c9_colors.copy()
    color_dict["other"] = "darkgrey"

    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)

    for ethnicity in populations:
        plt.style.use('seaborn-v0_8-whitegrid')
        fig, ax = plt.subplots(figsize=(8, 6))

        # Set transparent background
        fig.patch.set_alpha(0.0)
        ax.patch.set_alpha(0.0)

        pop_df = collect_population_samples(df_dipl_f, shape, sample_diplotype_func, total_samples, ethnicity)
        if pop_df.empty:
            plt.close(fig)
            continue

        pop_df["genotype_interest"] = pop_df["genotype"].apply(
            lambda g: g if g in genotypes_of_interest else "other"
        )

        x_data = pop_df["activity"].dropna().values
        mean_val = np.mean(x_data)  # Calculate mean instead of mode

        sns.histplot(
            data=pop_df,
            x="activity",
            hue="genotype_interest",
            multiple="stack",
            hue_order=["other", "*3/*3", "*1/*3", "*1/*2", "*1/*1"],
            palette=color_dict,
            alpha=1.0,  # Changed from 0.7 to 1.0
            bins=30,
            binrange=(0, 4),
            stat="density",
            legend=False,
            ax=ax,
        )

        sns.kdeplot(
            data=pop_df,
            x="activity",
            color="black",
            linewidth=2,
            ax=ax
        )

        # Draw mean line
        ax.axvline(mean_val, color="black", linestyle="--", linewidth=2, label=None)

        # Add mean annotation at the top of the plot
        y_max = ax.get_ylim()[1]
        ax.text(mean_val, y_max * 0.9, f'Mean: {mean_val:.2f}',
                horizontalalignment='left', verticalalignment='top', fontsize=30, fontweight='bold')

        ax.set_ylabel("", fontsize=20, fontweight='bold')
        ax.set_xlabel("", fontsize=20, fontweight='bold')

        # Remove grid
        ax.grid(False)
        ax.tick_params(axis='both', which='major', labelsize=25)
        ax.set_xlim(0, 2.5)

        for spine in ax.spines.values():
            spine.set_visible(False)

        plt.tight_layout()

        if output_dir is not None:
            fig_path = output_dir / f"ethnicity_activity_{ethnicity.replace('/', '_').replace(' ', '_')}.svg"
            fig.savefig(fig_path, format='svg', bbox_inches="tight", transparent=True)

        plt.close(fig)


def create_genotype_legend_svg(output_path: Path = None):
    """Create a standalone SVG legend for CYP2C9 genotypes."""
    # Set up colors
    color_dict = cyp2c9_colors.copy()
    color_dict["other"] = "darkgrey"

    # Create figure with transparent background
    fig = plt.figure(figsize=(2.5, 3))
    fig.patch.set_alpha(0.0)

    # Create legend elements
    legend_elements = [
        mpatches.Patch(facecolor=color_dict["*1/*1"], edgecolor='black', label="*1/*1"),
        mpatches.Patch(facecolor=color_dict["*1/*2"], edgecolor='black', label="*1/*2"),
        mpatches.Patch(facecolor=color_dict["*1/*3"], edgecolor='black', label="*1/*3"),
        mpatches.Patch(facecolor=color_dict["*3/*3"], edgecolor='black', label="*3/*3"),
        mpatches.Patch(facecolor=color_dict["other"], edgecolor='black', label="Other")
    ]

    # Add title
    plt.figtext(0.5, 0.95, "Genotypes", ha="center", fontsize=14, fontweight="bold")

    # Create the legend
    fig.legend(
        handles=legend_elements,
        loc='center',
        frameon=False,
        fontsize=12
    )

    # Save the figure
    if output_path is not None:
        plt.savefig(output_path, format='svg', bbox_inches="tight", transparent=True)

    plt.close(fig)


if __name__ == "__main__":
    # Load parameters and data
    lognorm_pars = get_lognorm_pars()
    shape_all = lognorm_pars["All"]["shape"]

    excel_path = DATA_PATH_GLIMEPIRIDE / "cyp2c9_data.xlsx"
    df_act, df_allele_f, df_dipl_f = read_and_filter_excel(excel_path)
    df_dipl_f = df_dipl_f.sort_values(by="mean_freq", ascending=False)

    output_figure_dir = glimepiride.RESULTS_PATH_ANALYSES / "_figures"
    output_figure_dir.mkdir(parents=True, exist_ok=True)

    analyze_sampled_genotypes(
        df_dipl_f,
        shape_all,
        sample_diplotype,
        total_samples=100000
    )

    stacked_hist_ethnicites(
        df_dipl_f=df_dipl_f,
        shape=shape_all,
        sample_diplotype_func=sample_diplotype,
        total_samples=100000,
        genotypes_of_interest=genotypes_of_interest,
        fig_path = output_figure_dir / "ethnicity_activity_distribution.png"
    )

    boxplots_ethnicities(
        df_dipl_f=df_dipl_f,
        shape=shape_all,
        sample_diplotype_func=sample_diplotype,
        total_samples=100000,
        fig_path = output_figure_dir / "ethnicity_activity_boxplots.png"
    )

    # Individual stacked histograms for world map
    # stacked_hist_ethnicities_individual(
    #     df_dipl_f=df_dipl_f,
    #     shape=shape_all,
    #     sample_diplotype_func=sample_diplotype,
    #     total_samples=1000000,
    #     genotypes_of_interest=genotypes_of_interest,
    #     output_dir=output_figure_dir
    # )
    # create_genotype_legend_svg(output_figure_dir / "genotype_legend.svg")