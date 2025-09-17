"""Module for determining distributions of CYP2C9 activities."""

from pathlib import Path
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from sbmlutils.console import console
from pkdb_models.models.glimepiride import DATA_PATH_GLIMEPIRIDE
import pkdb_models.models.glimepiride as glimepiride
from pkdb_models.models.glimepiride.experiments.base_experiment import GlimepirideSimulationExperiment
import os
from pathlib import Path
cyp2c9_allele_activity = GlimepirideSimulationExperiment.cyp2c9_allele_activity

_override = os.getenv("GLIMEPIRIDE_RESULTS_PATH")
if _override:
    p = Path(_override).resolve()
    p.mkdir(parents=True, exist_ok=True)
    glimepiride.RESULTS_PATH = p
    glimepiride.RESULTS_PATH_SIMULATION = p / "simulation"
    glimepiride.RESULTS_PATH_ANALYSES = p / "cyp2c9_analyses"
    glimepiride.RESULTS_PATH_FIT = p / "fit"

# Set random seed for reproducibility
np.random.seed(1234)

def load_data(verbose=True) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    """Load CYP2C9 data and return dataframe and clearance values."""
    df = pd.read_excel(DATA_PATH_GLIMEPIRIDE / "cyp2c9_data.xlsx", sheet_name="Yang2012_Tab5", comment="#")

    # Combined data
    all_data = df["CLint"].dropna().values
    data = {"All": all_data}

    if verbose:
        console.rule("RAW DATA STATISTICS (Yang2012)", align="left", style="bold white")
        d = data["All"]
        console.print(f"Sample size: {len(d)}")
        stats_table = {
            "Statistic": ["Mean", "Median", "SD", "Min", "Max", "Q1", "Q3"],
            "Value": [
                f"{np.mean(d):.2f}",
                f"{np.median(d):.2f}",
                f"{np.std(d):.2f}",
                f"{np.min(d):.2f}",
                f"{np.max(d):.2f}",
                f"{np.percentile(d, 25):.2f}",
                f"{np.percentile(d, 75):.2f}"
            ]
        }
        console.print(pd.DataFrame(stats_table).to_string(index=False))
        console.print("\n")

    return df, data


def fit_lognorm(data: dict[str, np.ndarray]) -> dict[str, dict]:
    """Fit lognormal distribution for combined data."""
    lognormal_pars = {}
    ethnicity = "All"
    shape, loc, scale = stats.lognorm.fit(data[ethnicity], floc=0)
    lognormal_pars[ethnicity] = {"shape": shape, "loc": loc, "scale": scale}
    return lognormal_pars


def lognorm_stats(shape: float, scale: float):
    """Compute arithmetic mean, SD, and CV for a lognormal distribution."""
    mean_val = scale * np.exp((shape ** 2) / 2)
    sd_val = mean_val * np.sqrt(np.exp(shape ** 2) - 1)
    cv_val = sd_val / mean_val
    return mean_val, sd_val, cv_val


def get_lognorm_pars() -> dict[str, dict]:
    """Return lognormal parameters for use by other scripts."""
    df, data = load_data(verbose=False)
    return fit_lognorm(data)


def plot_fit(data: dict[str, np.ndarray], bins: int, fig_path: Path) -> dict[str, dict]:
    """Fit and plot lognormal distribution for 'All' data."""
    lognormal_pars = fit_lognorm(data)

    fig, ax = plt.subplots(1, 1, figsize=(10, 7), sharey=True, sharex=True, layout="constrained")

    # Set grid style
    plt.style.use('seaborn-v0_8-whitegrid')
    ax.grid(True, linestyle='-', alpha=0.6, color='lightgrey')

    # "All" data only
    ethnicity = "All"
    d = data[ethnicity]
    shape = lognormal_pars[ethnicity]["shape"]
    loc = lognormal_pars[ethnicity]["loc"]
    scale = lognormal_pars[ethnicity]["scale"]

    # Calculate statistics
    mean_v, sd_v, cv_v = lognorm_stats(shape, scale)
    mode_v = scale * np.exp(-shape**2)

    # Kolmogorov-Smirnov test
    D, p_value = stats.kstest(d, lambda x: stats.lognorm.cdf(x, shape, loc=0, scale=scale))

    # Print statistics
    console.rule("LOGNORMAL DISTRIBUTION PARAMETERS", align="left", style="bold white")
    params_table = {
        "Parameter": ["Shape (σ)", "Scale (s)", "Location", "Arithmetic Mean", "SD", "CV", "Mode", "K-S test D",
                      "K-S test p-value"],
        "Value": [
            f"{shape:.2f}",
            f"{scale:.2f}",
            f"{loc:.2f}",
            f"{mean_v:.2f}",
            f"{sd_v:.2f}",
            f"{cv_v:.2f}",
            f"{mode_v:.2f}",
            f"{D:.2f}",
            f"{p_value:.2f}"
        ]
    }
    console.print(pd.DataFrame(params_table).to_string(index=False))
    console.print("\n")

    # Setup plot parameters
    x_max = 1000
    bin_edges = np.linspace(0, x_max, bins + 1)
    x_vals = np.linspace(0, x_max, 500)

    # Plot histogram
    ax.hist(
        d,
        bins=bin_edges,
        color="black",
        density=True,
        label="Combined Data",
        histtype='step',
        linewidth=1.5
    )

    # Plot fitted PDF
    pdf_vals = stats.lognorm.pdf(x_vals, shape, loc, scale)
    ax.plot(
        x_vals,
        pdf_vals,
        color="black",
        linestyle="solid",
        linewidth=3,
        label=f"Lognormal Fit"
    )

    # Set labels and styling
    ax.set_title("Intrinsic Clearance Distribution", fontsize=22, fontweight='bold')
    ax.set_xlabel("Intrinsic Clearance [mL/min]", fontsize=16, fontweight='bold')
    ax.set_ylabel("Probability Density", fontsize=16, fontweight='bold')
    ax.set_xlim(0, x_max)

    # Add legend
    legend = ax.legend(fontsize=14, frameon=True, framealpha=0.9, edgecolor='black')
    for text in legend.get_texts():
        text.set_fontweight('bold')

    # Set tick parameters
    ax.tick_params(axis='both', labelsize=12)

    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close()

    return lognormal_pars


def plot_activity(lognormal_pars: dict[str, dict], allele_activity: dict[str, float], fig_path: Path) -> None:
    """Plot activity distributions for alleles using lognormal model."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 7), sharey=True)

    # Set style
    plt.style.use('seaborn-v0_8-whitegrid')
    ax.grid(True, linestyle='--', alpha=0.5, color='lightgrey')

    x_vals = np.linspace(0, 2.5, 1000)
    color_map = {"*1": "black", "*2": "#bc66ff", "*3": "#440c6f"}
    label_map = {"*1": "*1", "*2": "*2", "*3": "*3"}


    # Table to store statistics for all alleles
    allele_stats = {"Allele": [], "Activity": [], "Mean": [], "Median": [], "SD": [], "Mode": []}

    # Use shape from combined data
    shape = lognormal_pars["All"]["shape"]

    # Plot PDFs and histograms for each allele
    for allele, activity in allele_activity.items():
        color = color_map[allele]
        scale_new = activity / np.exp(shape ** 2 / 2)

        # Plot fitted PDF
        pdf_vals = stats.lognorm.pdf(x_vals, s=shape, loc=0, scale=scale_new)
        ax.plot(x_vals, pdf_vals, color=color, linewidth=3, label=f"{label_map[allele]} - Lognorm Fit")

        # Generate and plot samples
        samples = stats.lognorm.rvs(s=shape, loc=0, scale=scale_new, size=10000)
        ax.hist(
            samples,
            bins=np.linspace(0, 4, num=20),
            color=color,
            density=True,
            label=f"{label_map[allele]} - Samples",
            histtype='step',
            linewidth=1.5
        )

        # Calculate statistics
        mean = np.mean(samples)
        median = np.median(samples)
        sd = np.std(samples)
        mode = scale_new * np.exp(-shape ** 2)

        # Add to statistics table
        allele_stats["Allele"].append(allele)
        allele_stats["Activity"].append(f"{activity:.2f}")
        allele_stats["Mean"].append(f"{mean:.2f}")
        allele_stats["Median"].append(f"{median:.2f}")
        allele_stats["SD"].append(f"{sd:.2f}")
        allele_stats["Mode"].append(f"{mode:.2f}")

        # Add mean line and annotation
        ax.axvline(mean, color=color, linestyle="dotted", linewidth=2)
        y_pos = stats.lognorm.pdf(mean, s=shape, loc=0, scale=scale_new) * 0.7
        ax.annotate(f"Mean: {mean:.2f}", xy=(mean, y_pos), xytext=(mean + 0.05, y_pos + 0.7),
                    fontsize=13, color=color, fontweight='bold')

    # Print allele statistics table
    console.rule("ALLELE ACTIVITY STATISTICS", align="left", style="bold white")
    console.print(f"Sample Size: {len(samples)}")
    console.print(pd.DataFrame(allele_stats).to_string(index=False))
    console.print("\n")

    # Set labels and styling
    ax.set_title("CYP2C9 Allele Activity Distribution", fontsize=22, fontweight='bold')
    ax.set_xlabel("Relative Enzymatic Activity", fontsize=16, fontweight='bold')
    ax.set_ylabel("Probability Density", fontsize=16, fontweight='bold')
    ax.set_xlim(0, 2.5)
    ax.tick_params(axis='both', labelsize=12)

    # Add legend
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(),
              fontsize=12, frameon=True, framealpha=0.9,
              edgecolor='black', loc='upper right')


    plt.tight_layout()
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close()


def visualize_genotype_distributions(shape: float, n_samples, fig_path: Path = None):
    """Visualize activity distributions for genotypes of interest."""
    plt.style.use('seaborn-v0_8-whitegrid')

    # Define genotypes and color mapping
    genotypes_of_interest = ["*1/*1", "*1/*2", "*1/*3", "*3/*3"]
    cyp2c9_colors = GlimepirideSimulationExperiment.cyp2c9_colors

    # Create table for genotype statistics
    genotype_stats = {"Genotype": [], "Mean": [], "Median": [], "SD": [], "Q1": [], "Q3": [], "Mode": []}

    def sample_allele(allele: str, shape: float, n=10000) -> np.ndarray:
        """Sample from an allele's distribution using mean-based scaling."""
        mean_val = cyp2c9_allele_activity[allele]
        scale = mean_val / np.exp(shape ** 2 / 2)
        return stats.lognorm.rvs(s=shape, scale=scale, size=n)

    def sample_diplotype(allele1: str, allele2: str, shape: float, n=10000) -> np.ndarray:
        """Sample genotype-level activity by averaging allele activities."""
        s1 = sample_allele(allele1, shape, n=n)
        s2 = sample_allele(allele2, shape, n=n)
        return (s1 + s2) / 2

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    # Plot distributions for each genotype
    for genotype in genotypes_of_interest:
        a1, a2 = genotype.split("/")
        samples = sample_diplotype(a1, a2, shape, n=n_samples)

        # Calculate statistics
        mean_val = np.mean(samples)
        median_val = np.median(samples)
        sd_val = np.std(samples)
        q1_val = np.percentile(samples, 25)
        q3_val = np.percentile(samples, 75)

        # Fit lognormal to samples
        shape_fit, loc_fit, scale_fit = stats.lognorm.fit(samples, floc=0)
        mode_val = scale_fit * np.exp(-shape_fit ** 2)

        # Add to statistics table
        genotype_stats["Genotype"].append(genotype)
        genotype_stats["Mean"].append(f"{mean_val:.2f}")
        genotype_stats["Median"].append(f"{median_val:.2f}")
        genotype_stats["SD"].append(f"{sd_val:.2f}")
        genotype_stats["Q1"].append(f"{q1_val:.2f}")
        genotype_stats["Q3"].append(f"{q3_val:.2f}")
        genotype_stats["Mode"].append(f"{mode_val:.2f}")

        # Plot histogram
        ax.hist(samples, bins=30, density=True, alpha=1.0,
                color=cyp2c9_colors[genotype],
                label=f"{genotype} - Samples", histtype='step', linewidth=2.0)

        # Plot lognormal fit
        x_vals = np.linspace(0, 3, 500)
        y_vals = stats.lognorm.pdf(x_vals, shape_fit, loc=loc_fit, scale=scale_fit)
        ax.plot(x_vals, y_vals, color=cyp2c9_colors[genotype],
                linewidth=3, label=f"{genotype} - Lognorm Fit")

        # Add mean line and annotation
        y_at_mean = stats.lognorm.pdf(mean_val, shape_fit, loc=loc_fit, scale=scale_fit)
        ax.axvline(x=mean_val, color=cyp2c9_colors[genotype],
                   linestyle='--', linewidth=1.5)
        ax.text(mean_val + 0.02, y_at_mean + 1.0,
                f'Mean: {mean_val:.2f}',
                color=cyp2c9_colors[genotype],
                fontweight='bold', fontsize=11,
                verticalalignment='center')

    # Print genotype statistics table
    console.rule("GENOTYPE ACTIVITY STATISTICS", align="left", style="bold white")
    console.print(f"Sample Size: {len(samples)}")
    console.print(pd.DataFrame(genotype_stats).to_string(index=False))
    console.print("\n")

    # Set labels and styling
    ax.set_title('Genotype Activity Distributions', fontsize=20, fontweight='bold')
    ax.set_xlabel('Activity', fontsize=16, fontweight='bold')
    ax.set_ylabel('Probability Density', fontsize=16, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5, color='lightgrey')
    ax.legend(fontsize=12)
    ax.set_xlim(0, 2.3)

    plt.tight_layout()
    if fig_path is not None:
        fig.savefig(fig_path, dpi=300, bbox_inches="tight")

    plt.show()


if __name__ == "__main__":
    # Load data
    df, data = load_data()
    output_figure_dir = glimepiride.RESULTS_PATH_ANALYSES / "_figures"
    output_figure_dir.mkdir(parents=True, exist_ok=True)

    # Fit and plot intrinsic clearance
    lognorm_pars = plot_fit(
        data=data,
        bins=8,
        fig_path = output_figure_dir / "cyp2c9_clearance_distribution.png"
    )

    # Plot allele activity distributions
    allele_activity = {
        "*1": 1.00,
        "*2": 0.68,
        "*3": 0.23
    }
    plot_activity(
        lognorm_pars,
        allele_activity,
        fig_path = output_figure_dir / "cyp2c9_allele_activity_distribution.png"
    )

    # Plot genotype distributions
    visualize_genotype_distributions(
        shape=lognorm_pars["All"]["shape"],
        n_samples=100000,
        fig_path = output_figure_dir / "cyp2c9_genotype_activity_distributions.png"
    )