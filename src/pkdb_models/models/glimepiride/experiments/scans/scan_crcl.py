"""Glimepiride parameter scan over creatinine clearance (CrCl)."""

from typing import Dict, Tuple, List
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from sbmlsim.simulation import Timecourse, TimecourseSim, ScanSim, Dimension
from sbmlsim.data import load_pkdb_dataframe
from pkdb_models.models.glimepiride import DATA_PATH_GLIMEPIRIDE
from pkdb_models.models.glimepiride.experiments.base_experiment import GlimepirideSimulationExperiment
from pkdb_models.models.glimepiride.helpers import run_experiments
from pkdb_models.models.glimepiride.glimepiride_pk import calculate_glimepiride_pk


class GlimepirideCrClScan(GlimepirideSimulationExperiment):
    """Scan from CrCl=0 to 150 ml/min by varying f_renal_function."""

    n_points = 100
    crcl_range = np.linspace(0, 150, n_points)
    normal_crcl = 110.0  # ml/min

    # Plotting
    figure_size_single = (28, 7)
    figure_size_extended = (32, 7)
    marker_size = 70
    line_width = 2
    capsize = 3
    font_size_label = 20
    font_size_legend = 14

    # Renal function categories
    renal_categories = [
        {"name": "Normal", "range": (90, 200), "color": "white"},
        {"name": "Mild impairment", "range": (60, 90), "color": "#66c2a4"},
        {"name": "Moderate impairment", "range": (30, 60), "color": "#2ca25f"},
        {"name": "Severe impairment", "range": (15, 30), "color": "#006d2c"},
        {"name": "Kidney failure", "range": (0, 15), "color": "#003917"}
    ]

    # Compounds
    compounds = {
        "glimepiride": {"display": "Glimepiride", "mr_attr": "gli"},
        "M1": {"display": "M1", "mr_attr": "m1"},
        "M2": {"display": "M2", "mr_attr": "m2"},
    }

    # Parameter labels
    parameter_units = {
        "cl": "Clearance [ml/min]",
        "cl_renal": "Renal Clearance [ml/min]",
        "vd": "Volume of Distribution [L]",
        "thalf": "Half-life [hr]",
        "auc": "AUC [ng*hr/ml]",
        "tmax": "Tmax [hr]",
        "cmax": "Cmax [ng/ml]",
    }

    # Regression plot indices for each figure
    regression_plots = {
        "Fig2": [0],
        "Fig3": [0, 1],
        "Fig4": [0, 1]
    }

    # Study data
    error_bar_data = {
        "Fig2": {
            3: ("creatinine_clearance", "glimepiride_cmax", [
                (77.7, 21.9, 359.2, 98.3),
                (27.4, 8.0, 205.3, 29.0),
                (9.4, 5.9, 194.0, 42.4)
            ]),
            4: ("creatinine_clearance", "glimepiride_tmax", [
                (77.7, 21.9, 1.9, 0.2),
                (27.4, 8.0, 2.7, 1.3),
                (9.4, 5.9, 2.2, 1.0)
            ]),
            5: ("creatinine_clearance", "glimepiride_auc", [
                (77.7, 21.9, 1357, 452),
                (27.4, 8.0, 622, 106),
                (9.4, 5.9, 622, 226)
            ])
        },
        "Fig3": {
            3: ("creatinine_clearance", "glimepiride-M1_cmax", [
                (77.7, 21.9, 70.8, 14),
                (27.4, 8.0, 93.0, 12.5),
                (9.4, 5.9, 103.6, 24.1)
            ]),
            4: ("creatinine_clearance", "glimepiride-M1_tmax", [
                (77.7, 21.9, 3.8, 0.8),
                (27.4, 8.0, 3.2, 0.8),
                (9.4, 5.9, 4.1, 2.3)
            ])
        },
        "Fig4": {
            3: ("creatinine_clearance", "glimepiride-M2_cmax", [
                (77.7, 21.9, 21.8, 8.5),
                (27.4, 8.0, 42.0, 2.6),
                (9.4, 5.9, 61.7, 25.9)
            ]),
            4: ("creatinine_clearance", "glimepiride-M2_tmax", [
                (77.7, 21.9, 4.8, 1.5),
                (27.4, 8.0, 5.7, 1.2),
                (9.4, 5.9, 7.0, 1.2)
            ])
        }
    }

    # Figure configuration
    figure_config = {
        "Fig2": {"title": "Impact of Renal Function on Glimepiride Pharmacokinetics", "n_plots": 6},
        "Fig3": {"title": "Impact of Renal Function on M1 Pharmacokinetics", "n_plots": 5},
        "Fig4": {"title": "Impact of Renal Function on M2 Pharmacokinetics", "n_plots": 5}
    }

    def simulations(self) -> Dict[str, ScanSim]:
        """Create scan over f_renal_function."""
        Q_ = self.Q_

        # Calculate f_renal_function values from CrCl
        f_renal_values = self.crcl_range / self.normal_crcl

        return {
            "scan_crcl": ScanSim(
                simulation=TimecourseSim(
                    Timecourse(
                        start=0,
                        end=24 * 60,
                        steps=5000,
                        changes={
                            **self.default_changes(),
                            "PODOSE_gli": Q_(3, "mg"),
                            "f_cirrhosis": Q_(0.0, "dimensionless"),
                        },
                    )
                ),
                dimensions=[
                    Dimension(
                        "dim_crcl_scan",
                        changes={
                            "KI__f_renal_function": Q_(f_renal_values, "dimensionless")
                        },
                    )
                ],
            )
        }

    def data(self) -> Dict:
        """Select model variables to store."""
        self.add_selections_data(
            selections=[
                "time", "PODOSE_gli",
                "[Cve_gli]", "[Cve_m1]", "[Cve_m2]",
                "Aurine_m1", "Aurine_m2",
                "Afeces_gli", "Afeces_m1", "Afeces_m2",
                "KI__crcl",
            ]
        )
        return {}

    def figures_mpl(self) -> Dict[str, plt.Figure]:
        """Generate and return figures."""
        # Calculate PK parameters
        pk_df = self._calculate_pk_parameters()

        # Convert units
        pk_df_converted = self._convert_units(pk_df)

        # Generate plots for each compound
        figures = {}
        figure_map = {"Fig2": "glimepiride", "Fig3": "M1", "Fig4": "M2"}

        for fig_id, compound in figure_map.items():
            fig = self._create_crcl_figure(fig_id, pk_df_converted)
            figures[f"pk_parameters_crcl_{compound}"] = fig

        return figures

    def _calculate_pk_parameters(self) -> pd.DataFrame:
        """Calculate PK parameters from simulation results."""
        xres = self.results["task_scan_crcl"]

        # Calculate PK parameters
        pk_df = calculate_glimepiride_pk(experiment=self, xres=xres)

        # Add CrCl values
        num_compounds = len(pk_df['compound'].unique())
        pk_df["crcl"] = np.tile(self.crcl_range, num_compounds)

        return pk_df

    def _convert_units(self, pk_df: pd.DataFrame) -> pd.DataFrame:
        """Convert PK parameters to standard units for plotting."""
        df = pk_df.copy()

        # Add converted columns for each compound
        for compound, info in self.compounds.items():
            mask = df['compound'] == compound
            if not mask.any():
                continue

            # Get molecular weight
            mr = getattr(self.Mr, info["mr_attr"]).magnitude

            # Clearance conversions
            if "cl" in df.columns:
                if compound == 'glimepiride':
                    # liter*milligram/millimole/minute to ml/min
                    df.loc[mask, "cl_mlmin"] = (df.loc[mask, "cl"] / mr) * 1000
                else:
                    # l/min to ml/min
                    df.loc[mask, "cl_mlmin"] = df.loc[mask, "cl"] * 1000

            # Renal clearance (M1/M2 only)
            if "cl_renal" in df.columns and compound in ['M1', 'M2']:
                df.loc[mask, "cl_renal_mlmin"] = df.loc[mask, "cl_renal"] * 1000

            # Volume of distribution
            if "vd" in df.columns:
                df.loc[mask, "vd_L"] = df.loc[mask, "vd"] / mr

            # Time conversions (minutes to hours)
            for param in ["thalf", "tmax"]:
                if param in df.columns:
                    df.loc[mask, f"{param}_h"] = df.loc[mask, param] / 60

            # Concentration conversions
            if "auc" in df.columns:
                df.loc[mask, "auc_ng_h_ml"] = df.loc[mask, "auc"] * mr * 1000 / 60

            if "cmax" in df.columns:
                df.loc[mask, "cmax_ng_ml"] = df.loc[mask, "cmax"] * mr * 1000

        return df

    def _create_crcl_figure(self, fig_id: str, pk_df: pd.DataFrame) -> plt.Figure:
        """Create figure for CrCl dependency."""
        # Get configuration
        config = self.figure_config[fig_id]
        n_plots = config["n_plots"]

        # Create figure
        fig_width = self.figure_size_extended[0] if fig_id == "Fig2" else self.figure_size_single[0]
        fig, axes = plt.subplots(1, n_plots, figsize=(fig_width, 7))
        fig.subplots_adjust(bottom=0.15)

        # Load experimental data
        df_exp = load_pkdb_dataframe(f"Rosenkranz1996a_{fig_id}", data_path=DATA_PATH_GLIMEPIRIDE)
        df_single = df_exp[df_exp["intervention"] == "GLI3"]
        df_multi = df_exp[df_exp["intervention"] == "GLI_multi"]

        # Plot data from dataframe first
        grouped = df_single.groupby(["x_label", "y_label"])
        for idx, ((x_label, y_label), group_single) in enumerate(grouped):
            if idx >= len(axes):
                break

            ax = axes[idx]

            # Get multi-dose data
            group_multi = df_multi[
                (df_multi["x_label"].str.replace("_multi", "") == x_label) &
                (df_multi["y_label"].str.replace("_multi", "") == y_label)
            ]

            # Plot the subplot
            self._plot_dataframe_data(ax, x_label, y_label, group_single, group_multi,
                                     pk_df, fig_id, idx)

        # Add error bar data at specific indices
        error_data = self.error_bar_data.get(fig_id, {})
        for idx, (x_label, y_label, data) in error_data.items():
            if idx < len(axes):
                ax = axes[idx]
                # Clear and replot if this index was already used
                ax.clear()
                self._plot_error_bar_subplot(ax, x_label, y_label, data, pk_df, fig_id, idx)

        # Add renal function legend
        self._add_renal_legend(fig)

        plt.tight_layout(rect=[0, 0.13, 1, 0.95])
        return fig

    def _plot_dataframe_data(self, ax, x_label: str, y_label: str, group_single: pd.DataFrame,
                            group_multi: pd.DataFrame, pk_df: pd.DataFrame, fig_id: str, idx: int) -> None:
        """Plot subplot using data from dataframe."""
        # Add renal bands
        if "creatinine_clearance" in x_label:
            self._add_renal_bands(ax)

        # Plot scatter data
        ax.scatter(group_single["x"], group_single["y"],
                  color="black", marker="D", edgecolor="k", alpha=0.8,
                  label="Single dose", s=self.marker_size, zorder=5)
        ax.scatter(group_multi["x"], group_multi["y"],
                  color="black", marker="s", edgecolor="k", alpha=0.8,
                  label="Multiple dose", s=self.marker_size, zorder=5)

        # Add simulation line
        self._add_simulation_line(ax, pk_df, x_label, y_label)

        # Add regression if needed
        if idx in self.regression_plots.get(fig_id, []):
            self._add_regression_line(ax, group_single, group_multi, y_label)

        # Format axis
        self._format_axis(ax, x_label, y_label)

        # Add legend
        ax.legend(fontsize=self.font_size_legend)

    def _plot_error_bar_subplot(self, ax, x_label: str, y_label: str, error_data: List[Tuple],
                               pk_df: pd.DataFrame, fig_id: str, idx: int) -> None:
        """Plot subplot with error bar data."""
        # Add renal bands
        if "creatinine_clearance" in x_label:
            self._add_renal_bands(ax)

        # Add simulation line
        self._add_simulation_line(ax, pk_df, x_label, y_label)

        # Add error bar data
        self._add_study_data(ax, error_data)

        # Format axis
        self._format_axis(ax, x_label, y_label)

        # Add legend
        ax.legend(fontsize=self.font_size_legend)

    def _add_renal_bands(self, ax) -> None:
        """Add colored bands for renal function categories."""
        for category in self.renal_categories:
            lower, upper = category["range"]
            ax.axvspan(lower, upper, alpha=0.2, color=category["color"], zorder=0)

    def _add_simulation_line(self, ax, pk_df: pd.DataFrame, x_label: str, y_label: str) -> None:
        """Add simulation line to plot."""
        if "creatinine_clearance" not in x_label:
            return

        # Parse y_label to get compound and parameter
        compound, param = self._parse_label(y_label)
        if not compound or not param:
            return

        # Get column name
        col_name = self._get_column_name(compound, param)
        if col_name not in pk_df.columns:
            return

        # Plot data
        df_filtered = pk_df[pk_df["compound"] == compound].sort_values("crcl")
        if not df_filtered.empty:
            ax.plot(df_filtered["crcl"], df_filtered[col_name],
                   color="black", linewidth=self.line_width,
                   label="Simulation", zorder=4)

    def _add_regression_line(self, ax, group_single: pd.DataFrame,
                           group_multi: pd.DataFrame, y_label: str) -> None:
        """Add regression line to plot."""
        # Combine data
        all_x = pd.concat([group_single["x"], group_multi["x"]])
        all_y = pd.concat([group_single["y"], group_multi["y"]])

        if len(all_x) < 2:
            return

        # Determine if regression should go through origin
        force_origin = not ("clearance" in y_label and "renal" not in y_label)

        # Calculate regression
        if force_origin:
            slope = np.sum(all_x * all_y) / np.sum(all_x ** 2)
            intercept = 0
            formula = f"y={slope:.2f}x"
        else:
            slope, intercept, _, _, _ = stats.linregress(all_x, all_y)
            formula = f"y={slope:.2f}x+{intercept:.2f}"

        # Plot line
        x_range = np.linspace(0, all_x.max(), 100)
        y_range = slope * x_range + intercept
        ax.plot(x_range, y_range, color="black",
               label=f"Regression ({formula})",
               linestyle="--", linewidth=2)

    def _add_study_data(self, ax, study_data: List[Tuple]) -> None:
        """Add study data error bars."""
        for idx, (x, xerr, y, yerr) in enumerate(study_data):
            label = "Study data" if idx == 0 else None
            ax.errorbar(x, y, xerr=xerr, yerr=yerr,
                       fmt='s', markersize=8,
                       markerfacecolor='black', markeredgecolor='black',
                       ecolor='black', capsize=self.capsize,
                       label=label)

    def _format_axis(self, ax, x_label: str, y_label: str) -> None:
        """Format axis with labels and grid."""
        ax.set_xlabel(self._get_axis_label(x_label),
                     fontsize=self.font_size_label, fontweight='bold')
        ax.set_ylabel(self._get_axis_label(y_label),
                     fontsize=self.font_size_label, fontweight='bold')
        ax.grid(True, linestyle='--', alpha=0.3, color='gray')

        if "creatinine_clearance" in x_label:
            ax.set_xlim(0, 130)

    def _add_renal_legend(self, fig) -> None:
        """Add legend for renal function categories."""
        patches = []
        for category in self.renal_categories:
            lower, upper = category["range"]
            patch = mpatches.Patch(
                color=category["color"],
                alpha=0.2,
                label=f"{category['name']} ({lower}-{upper} ml/min)"
            )
            patches.append(patch)

        fig.legend(handles=patches, loc='lower center',
                  ncol=len(patches), bbox_to_anchor=(0.5, 0.03),
                  frameon=True, fontsize=20)

    def _parse_label(self, label: str) -> Tuple[str, str]:
        """Parse compound and parameter from label."""
        # Handle special compound names
        if "glimepiride-M1" in label:
            compound = "M1"
            param = label.replace("glimepiride-M1_", "")
        elif "glimepiride-M2" in label:
            compound = "M2"
            param = label.replace("glimepiride-M2_", "")
        elif "glimepiride" in label:
            compound = "glimepiride"
            param = label.replace("glimepiride_", "")
        else:
            return None, None

        return compound, param

    def _get_column_name(self, compound: str, param: str) -> str:
        """Get the appropriate column name for a compound/parameter combination."""
        # Map parameter names to column suffixes
        param_map = {
            "clearance": "cl_mlmin",
            "clearance_renal": "cl_renal_mlmin",
            "vd": "vd_L",
            "thalf": "thalf_h",
            "auc": "auc_ng_h_ml",
            "tmax": "tmax_h",
            "cmax": "cmax_ng_ml",
        }

        return param_map.get(param, param)

    def _get_axis_label(self, label: str) -> str:
        """Get formatted axis label with units."""
        if "creatinine_clearance" in label:
            return "Creatinine Clearance [ml/min]"

        # Parse compound and parameter
        compound, param = self._parse_label(label)
        if not compound or not param:
            return label

        # Get parameter label
        param_label = self.parameter_units.get(param, param)

        # Format with compound name
        compound_name = self.compounds.get(compound, {}).get("display", compound)
        return f"{compound_name} {param_label}"


if __name__ == "__main__":
    output_dir = "GlimepirideCrClScan"
    run_experiments(
        experiment_classes=GlimepirideCrClScan,
        output_dir=output_dir,
        save_results=True
    )