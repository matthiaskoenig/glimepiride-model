"""Cirrhosis parameter scan for glimepiride pharmacokinetics."""

from typing import Dict, Tuple
import numpy as np
import pandas as pd
from sbmlsim.plot.serialization_matplotlib import plt, FigureMPL
from sbmlsim.simulation import Timecourse, TimecourseSim, ScanSim, Dimension
from pkdb_models.models.glimepiride.experiments.base_experiment import GlimepirideSimulationExperiment
from pkdb_models.models.glimepiride.helpers import run_experiments
from pkdb_models.models.glimepiride.glimepiride_pk import calculate_glimepiride_pk


class GlimepirideCirrhosisScan(GlimepirideSimulationExperiment):
    """Scan of cirrhosis effects on glimepiride pharmacokinetics."""

    cirrhosis_range = np.linspace(0, 0.99, 70)

    # Cirrhosis severity regions
    cirrhosis_regions = [
        (0.3, 0.6, '#74a9cf', "Mild cirrhosis (CTP A)"),
        (0.6, 0.8, '#2b8cbe', "Moderate cirrhosis (CTP B)"),
        (0.8, 1.01, '#045a8d', "Severe cirrhosis (CTP C)")
    ]

    # PK parameters
    pk_parameters = ["cmax", "tmax", "auc", "thalf"]
    pk_labels = {
        "cmax": "Cmax [ng/ml]",
        "tmax": "Tmax [hr]",
        "auc": "AUC [ng*hr/ml]",
        "thalf": "Half-life [hr]"
    }

    # Unit conversions
    unit_conversions = {
        "tmax": {"factor": 1/60, "unit": "hr"},
        "thalf": {"factor": 1/60, "unit": "hr"},
        "cmax": {"factor_multiplier": 1000, "unit": "ng/ml"},
        "auc": {"factor_multiplier": 1000/60, "unit": "ng*hr/ml"},
    }

    # Compounds
    compounds = {
        "glimepiride": "gli",
        "M1": "m1",
        "M2": "m2"
    }

    # Clinical reference data
    clinical_data = {
        "cmax": {"min": 72, "max": 187},
        "tmax": {"min": 0.5, "max": 4.0},
        "auc": {"min": 213, "max": 529},
        "aurine_m1": {"min": 59, "max": 442},
        "aurine_m2": {"min": 80, "max": 227}
    }

    healthy_data = {
        "cmax": {"mean": 103, "sd": 34},
        "tmax": {"mean": 2.3, "sd": 0.5},
        "auc": {"mean": 326, "sd": 97},
        "aurine_m1": {"mean": 320, "sd": 40},
        "aurine_m2": {"mean": 160, "sd": 50}
    }

    # Plot configuration
    figure_size = (20, 15)
    line_width = 2.0
    marker_size = 8
    capsize = 7
    capthick = 1.5
    reference_color = 'black'

    # Reference cirrhosis data positioning
    cirrhosis_vertical_position = 0.65
    cirrhosis_horizontal_start = 0.3
    cirrhosis_horizontal_end = 1.0

    def simulations(self) -> Dict[str, ScanSim]:
        """Create cirrhosis scan simulation."""
        Q_ = self.Q_
        return {
            "cirrhosis_scan": ScanSim(
                simulation=TimecourseSim(
                    Timecourse(
                        start=0,
                        end=24 * 60,
                        steps=5000,
                        changes={
                            **self.default_changes(),
                            "KI__f_renal_function": Q_(1.0, "dimensionless"),
                            "f_cirrhosis": Q_(0.0, "dimensionless"),
                            "PODOSE_gli": Q_(1, "mg"),
                        },
                    )
                ),
                dimensions=[
                    Dimension(
                        "dim_scan",
                        changes={"f_cirrhosis": Q_(self.cirrhosis_range, "dimensionless")},
                    )
                ],
            )
        }

    def data(self) -> Dict:
        """Define model variables to record."""
        self.add_selections_data(
            selections=[
                "time", "PODOSE_gli", "f_cirrhosis",
                "[Cve_gli]", "[Cve_m1]", "[Cve_m2]",
                "Aurine_m1", "Aurine_m2",
                "Afeces_gli", "Afeces_m1", "Afeces_m2",
            ]
        )
        return {}

    def figures_mpl(self) -> Dict[str, FigureMPL]:
        """Calculate PK & urinary data, create and return figure object."""
        # Calculate PK parameters and urinary excretion
        pk_df, urinary_df = self._calculate_pk_and_urinary_data()

        # Create figure
        fig = self._create_cirrhosis_figure(pk_df, urinary_df)

        return {"pk_parameters_hepatic": fig}

    def _calculate_pk_and_urinary_data(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Calculate PK parameters and urinary excretion from simulation results."""
        # Get results
        xres = self.results["task_cirrhosis_scan"]

        # Calculate PK parameters
        pk_df = calculate_glimepiride_pk(experiment=self, xres=xres)

        # Add cirrhosis values
        num_compounds = len(pk_df['compound'].unique())
        if num_compounds > 0:
            pk_df["cirrhosis"] = np.tile(self.cirrhosis_range, num_compounds)

        # Convert units
        pk_df = self._convert_units(pk_df)

        # Calculate urinary excretion
        urinary_df = self._calculate_urinary_excretion(xres.xds)

        return pk_df, urinary_df

    def _convert_units(self, pk_df: pd.DataFrame) -> pd.DataFrame:
        """Convert PK parameters to standard units."""
        df = pk_df.copy()

        for compound, mr_attr in self.compounds.items():
            mask = df["compound"] == compound
            if not mask.any():
                continue

            # Get molecular weight
            mr = getattr(self.Mr, mr_attr).magnitude

            # Apply conversions
            for param, conversion in self.unit_conversions.items():
                if param not in df.columns:
                    continue

                # Calculate conversion factor
                if "factor_multiplier" in conversion:
                    factor = mr * conversion["factor_multiplier"]
                else:
                    factor = conversion["factor"]

                # Apply conversion
                df.loc[mask, param] = df.loc[mask, param] * factor
                df.loc[mask, f"{param}_unit"] = conversion["unit"]

        return df

    def _calculate_urinary_excretion(self, xds) -> pd.DataFrame:
        """Calculate urinary excretion data from xarray dataset."""
        if "Aurine_m1" not in xds or "Aurine_m2" not in xds or "_time" not in xds.dims:
            return pd.DataFrame(columns=["cirrhosis", "aurine_m1", "aurine_m2"])

        # Get final time point values
        aurine_m1_end = xds["Aurine_m1"].values[-1, :]
        aurine_m2_end = xds["Aurine_m2"].values[-1, :]

        # Convert to µg
        mr_m1 = self.Mr.m1.magnitude
        mr_m2 = self.Mr.m2.magnitude

        urinary_data = [
            {
                "cirrhosis": cirr_value,
                "aurine_m1": aurine_m1_end[idx] * mr_m1 * 1000,
                "aurine_m2": aurine_m2_end[idx] * mr_m2 * 1000,
            }
            for idx, cirr_value in enumerate(self.cirrhosis_range)
        ]

        return pd.DataFrame(urinary_data)

    def _create_cirrhosis_figure(self, pk_df: pd.DataFrame, urinary_df: pd.DataFrame) -> FigureMPL:
        """Create main cirrhosis dependency figure."""
        fig, axes = plt.subplots(nrows=3, ncols=4, figsize=self.figure_size)

        # Plot each compound
        compounds = ["glimepiride", "M1", "M2"]
        legend_handles = {"simulation": None, "reference": None, "regions": []}

        for row_idx, compound in enumerate(compounds):
            self._plot_compound_row(
                axes[row_idx], pk_df, urinary_df, compound, row_idx, legend_handles
            )

        # Add figure elements
        self._add_figure_legend(fig, legend_handles)
        fig.suptitle("Glimepiride Pharmacokinetics - Cirrhosis", fontsize=20, fontweight='bold')
        fig.tight_layout(rect=[0, 0.05, 1, 0.95])

        return fig

    def _plot_compound_row(self, axes_row, pk_df: pd.DataFrame, urinary_df: pd.DataFrame,
                           compound: str, row_idx: int, legend_handles: dict) -> None:
        """Plot all parameters for a single compound."""
        # Filter PK data for this compound
        df_compound = pk_df[pk_df["compound"] == compound].sort_values("cirrhosis")

        # Define which parameters to plot for each column
        col_params = self._get_column_parameters(compound, row_idx)

        for col_idx, param_info in enumerate(col_params):
            ax = axes_row[col_idx]

            # Add cirrhosis regions
            if row_idx == 0 and col_idx == 0:
                legend_handles["regions"] = self._add_cirrhosis_regions(ax)
            else:
                self._add_cirrhosis_regions(ax)

            # Plot data
            if param_info["type"] == "pk":
                self._plot_pk_parameter(ax, df_compound, compound, param_info, legend_handles)
            else:  # urinary
                self._plot_urinary_parameter(ax, urinary_df, compound, param_info, legend_handles)

            # Configure axis
            self._configure_axis(ax, compound, param_info, row_idx)

    def _get_column_parameters(self, compound: str, row_idx: int) -> list:
        """Get parameter configuration for each column."""
        # Default parameters for columns 0-3
        params = [
            {"type": "pk", "name": "cmax", "ylim": (0, 200) if row_idx == 0 else None},
            {"type": "pk", "name": "tmax", "ylim": (0, 5) if row_idx == 0 else (0, 6)},
            {"type": "pk", "name": "auc", "ylim": None},
            {"type": "pk", "name": "thalf", "ylim": None}
        ]

        # M1 and M2 have urinary excretion in column 3
        if compound in ["M1", "M2"]:
            params[3] = {
                "type": "urinary",
                "name": f"aurine_{compound.lower()}",
                "ylim": (0, 500) if compound == "M1" else (0, 250)
            }

        return params

    def _plot_pk_parameter(self, ax, df_compound: pd.DataFrame, compound: str,
                          param_info: dict, legend_handles: dict) -> None:
        """Plot PK parameter data."""
        param = param_info["name"]
        if df_compound.empty or param not in df_compound.columns:
            return

        # Plot simulation line
        line, = ax.plot(
            df_compound["cirrhosis"],
            df_compound[param],
            color=self.reference_color,
            linewidth=self.line_width,
            linestyle="-"
        )

        if legend_handles["simulation"] is None:
            legend_handles["simulation"] = line

        # Add reference data for glimepiride
        if compound == "glimepiride" and self._reference_data(param):
            ref_handle = self._add_clinical_data_reference(ax, param)
            if legend_handles["reference"] is None and ref_handle is not None:
                legend_handles["reference"] = ref_handle

    def _plot_urinary_parameter(self, ax, urinary_df: pd.DataFrame, compound: str,
                               param_info: dict, legend_handles: dict) -> None:
        """Plot urinary excretion data."""
        param = param_info["name"]
        if urinary_df.empty or param not in urinary_df.columns:
            return

        # Plot simulation line
        line, = ax.plot(
            urinary_df["cirrhosis"],
            urinary_df[param],
            color=self.reference_color,
            linewidth=self.line_width,
            linestyle="-"
        )

        if legend_handles["simulation"] is None:
            legend_handles["simulation"] = line

        # Add reference data
        if self._reference_data(param):
            ref_handle = self._add_clinical_data_reference(ax, param)
            if legend_handles["reference"] is None and ref_handle is not None:
                legend_handles["reference"] = ref_handle

    def _reference_data(self, param: str) -> bool:
        """Check if reference data should be shown for this parameter."""
        return param in self.clinical_data or param in self.healthy_data

    def _add_cirrhosis_regions(self, ax):
        """Add cirrhosis region markers."""
        patches = []
        for start, end, color, label in self.cirrhosis_regions:
            patch = ax.axvspan(start, end, alpha=0.4, color=color, zorder=-10)
            patches.append((patch, label))
        ax.set_xlim(-0.02, 1.01)
        return patches

    def _add_clinical_data_reference(self, ax, param_name):
        """Add clinical data reference lines/points."""
        clinical_ref = self.clinical_data.get(param_name)
        healthy_ref = self.healthy_data.get(param_name)
        legend_handle = None

        # Add clinical range (cirrhosis patients)
        if clinical_ref:
            min_val, max_val = clinical_ref["min"], clinical_ref["max"]
            mid_val = (min_val + max_val) / 2

            # Vertical error bar
            yerr_vertical = [[mid_val - min_val], [max_val - mid_val]]
            ax.errorbar(
                self.cirrhosis_vertical_position, mid_val,
                yerr=yerr_vertical, fmt='none',
                color=self.reference_color, capsize=self.capsize,
                capthick=self.capthick, linewidth=self.capthick,
                alpha=1.0, zorder=5
            )

            # Horizontal error bar
            x_mid = (self.cirrhosis_horizontal_start + self.cirrhosis_horizontal_end) / 2
            x_err = [[x_mid - self.cirrhosis_horizontal_start],
                     [self.cirrhosis_horizontal_end - x_mid]]
            ax.errorbar(
                x_mid, mid_val,
                xerr=x_err, fmt='none',
                color=self.reference_color, capsize=self.capsize,
                capthick=self.capthick, linewidth=self.capthick,
                alpha=1.0, zorder=5
            )

        # Add healthy reference
        if healthy_ref:
            mean_val, sd_val = healthy_ref["mean"], healthy_ref["sd"]
            healthy_point = ax.errorbar(
                0, mean_val, yerr=sd_val, fmt='s',
                color=self.reference_color, capsize=self.capsize,
                capthick=self.capthick, linewidth=self.capthick,
                alpha=1.0, markersize=self.marker_size,
                markerfacecolor=self.reference_color,
                markeredgecolor=self.reference_color, zorder=5
            )
            if legend_handle is None:
                legend_handle = healthy_point

        return legend_handle

    def _configure_axis(self, ax, compound: str, param_info: dict, row_idx: int) -> None:
        """Configure axis labels and appearance."""
        # X-label (bottom row only)
        if row_idx == 2:
            ax.set_xlabel("Cirrhosis Severity", fontsize=16, fontweight='bold')

        # Y-label
        if param_info["type"] == "pk":
            ylabel = f"{compound.capitalize()} {self.pk_labels[param_info['name']]}"
        else:
            ylabel = f"{compound} Urinary Excretion [µg]"
        ax.set_ylabel(ylabel, fontsize=16, fontweight='bold')

        # Grid
        ax.grid(True, linestyle='--', alpha=0.7)

        # Y-limits
        ylim = param_info.get("ylim")
        if ylim:
            ax.set_ylim(*ylim)
        else:
            ax.set_ylim(0, None)

    def _add_figure_legend(self, fig, legend_handles: dict) -> None:
        """Add legend to the figure."""
        handles, labels = [], []

        # Add simulation line
        if legend_handles["simulation"]:
            handles.append(legend_handles["simulation"])
            labels.append("Simulation")

        # Add reference data
        if legend_handles["reference"]:
            handles.append(legend_handles["reference"])
            labels.append("Study Data (Range/Mean ± SD)")

        # Add region patches
        for patch, label in legend_handles["regions"]:
            handles.append(patch)
            labels.append(label)

        # Create legend
        fig.legend(
            handles, labels,
            loc='lower center',
            bbox_to_anchor=(0.5, 0.02),
            ncol=len(handles),
            fontsize=16,
            frameon=True
        )


if __name__ == "__main__":
    output_dir = "GlimepirideCirrhosisScan"
    run_experiments(
        experiment_classes=GlimepirideCirrhosisScan,
        output_dir=output_dir,
        save_results=True
    )