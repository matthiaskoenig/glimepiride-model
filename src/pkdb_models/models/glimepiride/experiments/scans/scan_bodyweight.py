"""Bodyweight parameter scan for glimepiride pharmacokinetics."""

from typing import Dict
import pandas as pd
import numpy as np
from sbmlsim.plot.serialization_matplotlib import plt, FigureMPL
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from sbmlsim.simulation import Timecourse, TimecourseSim, ScanSim, Dimension
from pkdb_models.models.glimepiride.experiments.base_experiment import GlimepirideSimulationExperiment
from pkdb_models.models.glimepiride.helpers import run_experiments
from pkdb_models.models.glimepiride.glimepiride_pk import calculate_glimepiride_pk

class GlimepirideBodyweightScan(GlimepirideSimulationExperiment):
    """Scans bodyweight effects on glimepiride pharmacokinetics."""

    bodyweight_range = np.linspace(45, 170, 100)

    # Plotting
    figure_size = (16, 12)
    dpi = 300
    marker_size = 8
    line_width = 2
    error_cap_size = 5
    scatter_size = 25
    scatter_alpha = 0.8
    y_axis_padding = 1.2
    x_axis_padding = 0.05

    # Colors
    color_palette = ['#F5C6D0', '#BB4681', '#601038']

    # Compounds
    compounds = {
        "glimepiride": "gli",
        "M1": "m1",
        "M2": "m2"
    }

    # PK parameters
    pk_parameters = ["cmax", "tmax", "auc", "thalf"]
    pk_labels = ["Cmax", "Tmax", "AUC", "Half-life"]

    # Units
    unit_conversions = {
        "tmax": {"factor": 1/60, "unit": "hr"},
        "thalf": {"factor": 1/60, "unit": "hr"},
        "cmax": {"factor_multiplier": 1000, "unit": "ng/ml"},  # Will be multiplied by MW
        "auc": {"factor_multiplier": 1000/60, "unit": "ng*hr/ml"},  # Will be multiplied by MW
        "cl": {"factor": 60, "unit": "L/hr"},
    }

    # Shukla2004 study data (8mg dose)
    shukla_data = {
        "normal": {
            "weight": (72, 10.0),  # (mean, std)
            "glimepiride": {
                "cmax": (547, 218), "tmax": (2.89, 0.90),
                "auc": (3205, 1033), "thalf": (12.6, 12.8)
            },
            "M1": {
                "cmax": (180, 63), "tmax": (4.50, 0.78),
                "auc": (1514, 573), "thalf": (9.76, 3.24)
            },
            "M2": {
                "cmax": (50.3, 15.7), "tmax": (5.00, 1.14),
                "auc": (454, 130), "thalf": (7.09, 3.89)
            }
        },
        "obese": {
            "weight": (130, 35.5),
            "glimepiride": {
                "cmax": (410, 124), "tmax": (2.90, 0.89),
                "auc": (2818, 1112), "thalf": (8.89, 3.91)
            },
            "M1": {
                "cmax": (135, 52), "tmax": (4.33, 0.82),
                "auc": (1253, 344), "thalf": (11.6, 5.2)
            },
            "M2": {
                "cmax": (36.5, 14.7), "tmax": (5.36, 1.45),
                "auc": (316, 119), "thalf": (6.37, 4.87)
            }
        }
    }

    # Gu2010 AUC data (bodyweight vs AUC, originally for 2mg dose)
    gu2010_raw = np.array([
        [46.89, 1170.33], [47.88, 1008.24], [46.94, 909.34], [46.76, 788.46],
        [47.16, 793.96], [47.88, 755.49], [46.67, 684.07], [50.50, 1178.57],
        [50.99, 1096.15], [50.99, 1046.70], [50.77, 980.77], [50.32, 947.80],
        [50.68, 901.10], [50.41, 846.15], [53.83, 1192.31], [53.83, 1167.58],
        [56.17, 1126.37], [57.84, 1079.67], [53.69, 1016.48], [53.87, 975.27],
        [52.57, 917.58], [52.57, 879.12], [53.06, 760.99], [52.03, 747.25],
        [52.03, 708.79], [52.93, 670.33], [55.59, 799.45], [55.63, 667.58],
        [55.72, 640.11], [56.89, 615.38], [56.98, 538.46], [58.06, 717.03],
        [55.68, 1016.48], [56.17, 1010.99], [58.96, 1142.86], [59.01, 1115.38],
        [66.85, 1255.49], [66.89, 1197.80], [63.87, 1120.88], [63.92, 810.44],
        [63.92, 703.30], [63.92, 532.97], [66.94, 1065.93], [68.83, 997.25],
        [66.80, 780.22], [68.78, 714.29], [69.91, 892.86], [70.05, 870.88],
        [69.91, 832.42], [72.39, 934.07], [72.30, 865.38], [72.79, 796.70],
        [72.97, 530.22], [77.97, 898.35], [78.02, 851.65], [84.91, 725.27],
        [84.82, 695.06], [80.95, 1989.01], [81.00, 1881.87]
    ])

    def simulations(self) -> Dict[str, ScanSim]:
        """Create bodyweight scan simulation."""
        Q_ = self.Q_
        return {
            "bodyweight_scan": ScanSim(
                simulation=TimecourseSim(
                    Timecourse(
                        start=0,
                        end=24 * 60,
                        steps=5000,
                        changes={
                            **self.default_changes(),
                            "KI__f_renal_function": Q_(1.0, "dimensionless"),
                            "f_cirrhosis": Q_(0.0, "dimensionless"),
                            "PODOSE_gli": Q_(8, "mg"),
                        },
                    )
                ),
                dimensions=[
                    Dimension(
                        "dim_scan",
                        changes={"BW": Q_(self.bodyweight_range, "kg")},
                    )
                ],
            )
        }

    def data(self) -> Dict:
        """Specify model variables to record."""
        self.add_selections_data(
            selections=[
                "time", "PODOSE_gli", "BW",
                "[Cve_gli]", "[Cve_m1]", "[Cve_m2]",
                "Aurine_m1", "Aurine_m2",
                "Afeces_gli", "Afeces_m1", "Afeces_m2"
            ]
        )
        return {}

    def figures_mpl(self) -> Dict[str, FigureMPL]:
        """Create figures for bodyweight scan results."""
        # Calculate PK parameters
        pk_df = self._calculate_pk_parameters()

        # Create figure
        figure = self._create_pk_vs_bodyweight_figure(pk_df)

        return {"pk_parameters_bw": figure}

    def _calculate_pk_parameters(self) -> pd.DataFrame:
        """Calculate and prepare PK parameters."""
        # Get results
        xres = self.results["task_bodyweight_scan"]

        # Calculate PK parameters
        pk_df = calculate_glimepiride_pk(experiment=self, xres=xres)

        # Add bodyweight values
        num_compounds = len(pk_df['compound'].unique())
        pk_df["bodyweight"] = np.tile(self.bodyweight_range, num_compounds)

        # Convert units
        pk_df = self._convert_units(pk_df)

        return pk_df

    def _convert_units(self, pk_df: pd.DataFrame) -> pd.DataFrame:
        """Convert PK parameters to standard units."""
        df = pk_df.copy()

        for compound, mr_attr in self.compounds.items():
            mask = df["compound"] == compound

            # Get molecular weight
            mr = getattr(self.Mr, mr_attr).magnitude

            # Apply conversions
            for param, conversion in self.unit_conversions.items():
                # Calculate conversion factor
                if "factor_multiplier" in conversion:
                    factor = mr * conversion["factor_multiplier"]
                else:
                    factor = conversion["factor"]
                # Apply conversion
                df.loc[mask, param] = df.loc[mask, param] * factor
                df.loc[mask, f"{param}_unit"] = conversion["unit"]

        return df

    def _get_dose_adjusted_gu2010_data(self) -> np.ndarray:
        """Get Gu2010 AUC data adjusted for simulation dose."""
        dose_ratio = 8 / 2.0  # Simulation uses 8mg, original data was for 2mg
        adjusted_data = self.gu2010_raw.copy()
        adjusted_data[:, 1] *= dose_ratio
        return adjusted_data

    def _create_pk_vs_bodyweight_figure(self, pk_df: pd.DataFrame) -> FigureMPL:
        """Create main figure showing PK parameters vs bodyweight."""
        fig, axes = plt.subplots(
            nrows=3,
            ncols=4,
            figsize=self.figure_size,
            dpi=self.dpi
        )

        fig.suptitle('Glimepiride Pharmacokinetics - Bodyweight', fontsize=20, fontweight='bold')

        # Setup colormap
        cmap = LinearSegmentedColormap.from_list("bodyweight", self.color_palette)
        norm = plt.Normalize(vmin=self.bodyweight_range.min(), vmax=self.bodyweight_range.max())

        # Plot each compound
        compounds = ["glimepiride", "M1", "M2"]
        for row_idx, compound in enumerate(compounds):
            self._plot_compound_row(
                axes[row_idx], pk_df, compound, cmap, norm
            )

        # Add legend
        self._add_figure_legend(fig, cmap, norm)

        # Adjust layout
        plt.tight_layout(rect=[0, 0.05, 1, 0.95])

        return fig

    def _plot_compound_row(self, axes_row, pk_df: pd.DataFrame, compound: str,
                           cmap, norm) -> None:
        """Plot PK parameters for a compound."""
        # Filter data for this compound
        df_compound = pk_df[pk_df["compound"] == compound].sort_values("bodyweight")

        for col_idx, (param, label) in enumerate(zip(self.pk_parameters, self.pk_labels)):
            ax = axes_row[col_idx]

            # Plot simulation line
            if param in df_compound.columns:
                ax.plot(
                    df_compound["bodyweight"],
                    df_compound[param],
                    '-',
                    color="black",
                    linewidth=self.line_width,
                    zorder=100
                )

            # Plot study data
            self._add_study_data_points(ax, compound, param, cmap, norm)

            # Add scatter data for glimepiride AUC
            if compound == "glimepiride" and param == "auc":
                self._add_gu2010_scatter(ax, cmap, norm)

            # Configure axis
            self._configure_axis(ax, compound, param, label, col_idx)

            # Set y-limits
            self._set_y_limits(ax, df_compound, compound, param)

    def _add_study_data_points(self, ax, compound: str, param: str,
                              cmap, norm) -> None:
        """Add Shukla2004 study data points."""
        for group_name, group_data in self.shukla_data.items():
            # Get weight
            weight_mean, weight_std = group_data["weight"]
            # Get parameter value
            if compound in group_data and param in group_data[compound]:
                value_mean, value_std = group_data[compound][param]
                # Get color
                color = cmap(norm(weight_mean))
                # Plot point with error bars
                ax.errorbar(
                    weight_mean, value_mean,
                    xerr=weight_std, yerr=value_std,
                    fmt='s',
                    markersize=self.marker_size,
                    color=color,
                    ecolor=color,
                    capsize=self.error_cap_size,
                    markeredgecolor='black',
                    markeredgewidth=1.5,
                    zorder=150
                )

    def _add_gu2010_scatter(self, ax, cmap, norm) -> None:
        """Add Gu2010 AUC scatter data."""
        gu2010_data = self._get_dose_adjusted_gu2010_data()
        ax.scatter(
            gu2010_data[:, 0],
            gu2010_data[:, 1],
            c=gu2010_data[:, 0],
            cmap=cmap,
            norm=norm,
            marker='s',
            s=self.scatter_size,
            alpha=self.scatter_alpha,
            edgecolor='black',
            linewidth=0.5,
            zorder=125
        )

    def _configure_axis(self, ax, compound: str, param: str, label: str,
                       col_idx: int) -> None:
        """Configure axis labels and appearance."""
        # Title (first row only)
        if ax.get_subplotspec().rowspan.start == 0:
            ax.set_title(label, fontsize=14, fontweight='bold')

        # X-label (last row only)
        if ax.get_subplotspec().rowspan.start == 2:
            ax.set_xlabel("Bodyweight [kg]", fontsize=12, fontweight='bold')

        # Y-label
        unit = self.unit_conversions.get(param, {}).get("unit", "")
        ylabel = f"{compound.capitalize()} {label}"
        if unit:
            ylabel += f" [{unit}]"
        ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')

        # Grid
        ax.grid(True, linestyle='--', alpha=0.4)

        # X-axis limits
        x_range = self.bodyweight_range.max() - self.bodyweight_range.min()
        x_padding = x_range * self.x_axis_padding
        ax.set_xlim(
            self.bodyweight_range.min() - x_padding,
            self.bodyweight_range.max() + x_padding
        )

    def _set_y_limits(self, ax, df_compound: pd.DataFrame, compound: str,
                     param: str) -> None:
        """Set appropriate y-axis limits."""
        y_values = []

        # Add simulation values
        if param in df_compound.columns:
            y_values.extend(df_compound[param].dropna().tolist())

        # Add study values
        for group_data in self.shukla_data.values():
            if compound in group_data and param in group_data[compound]:
                mean, std = group_data[compound][param]
                y_values.extend([mean - std, mean + std])

        # Add Gu2010 data if applicable
        if compound == "glimepiride" and param == "auc":
            gu2010_data = self._get_dose_adjusted_gu2010_data()
            y_values.extend(gu2010_data[:, 1])

        # Set limits
        if y_values:
            ax.set_ylim(0, max(y_values) * self.y_axis_padding)
        else:
            ax.set_ylim(0, 1)

    def _add_figure_legend(self, fig, cmap, norm) -> None:
        """Add legend to the figure."""
        # Create legend elements
        handles = [
            Line2D([0], [0], color='black', lw=self.line_width, label='Simulation'),
            Line2D([0], [0], marker='s', color='w',
                   markerfacecolor=cmap(norm(self.shukla_data["normal"]["weight"][0])),
                   markeredgecolor='black', markersize=self.marker_size,
                   linestyle='None', label='Normal weight (Shukla2004)'),
            Line2D([0], [0], marker='s', color='w',
                   markerfacecolor=cmap(norm(self.shukla_data["obese"]["weight"][0])),
                   markeredgecolor='black', markersize=self.marker_size,
                   linestyle='None', label='Morbidly obese (Shukla2004)'),
            Line2D([0], [0], marker='s', color='w', markersize=6,
                   linestyle='None', label='AUC data (Gu2010)')
        ]

        # Add legend
        fig.legend(
            handles=handles,
            loc='lower center',
            bbox_to_anchor=(0.5, 0.01),
            ncol=4,
            fontsize=12,
            frameon=True
        )


if __name__ == "__main__":
    output_dir = "GlimepirideBodyweightScan"
    run_experiments(
        experiment_classes=GlimepirideBodyweightScan,
        output_dir=output_dir,
        save_results=True
    )