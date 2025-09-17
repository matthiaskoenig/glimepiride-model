"""Dose parameter scan for glimepiride pharmacokinetics."""

from typing import Dict
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from sbmlsim.simulation import Timecourse, TimecourseSim, ScanSim, Dimension
from sbmlsim.plot.serialization_matplotlib import plt, FigureMPL
from pkdb_models.models.glimepiride.experiments.base_experiment import GlimepirideSimulationExperiment
from pkdb_models.models.glimepiride.helpers import run_experiments
from pkdb_models.models.glimepiride.glimepiride_pk import calculate_glimepiride_pk
from pkdb_models.models.glimepiride import DATA_PATH_GLIMEPIRIDE


class GlimepirideDoseScan(GlimepirideSimulationExperiment):
    """Scans dose dependency in glimepiride pharmacokinetics."""

    # Configuration
    dose_range = np.linspace(0.1, 9.0, num=50)

    # Plotting configuration
    plot_config = {
        "figure_size": (16, 12),
        "dpi": 600,
        "compounds": ["glimepiride", "M1", "M2"],
        "pk_params": ["cmax", "tmax", "auc", "thalf"],
        "param_labels": ["Cmax", "Tmax", "AUC", "Half-life"],
    }

    def simulations(self) -> Dict[str, ScanSim]:
        """Create dose scan simulation."""
        Q_ = self.Q_

        return {
            "dose_scan": ScanSim(
                simulation=TimecourseSim(
                    Timecourse(
                        start=0,
                        end=24 * 60,
                        steps=2500,
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
                        changes={"PODOSE_gli": Q_(self.dose_range, "mg")},
                    )
                ],
            )
        }

    def data(self) -> Dict:
        """Record model variables."""
        self.add_selections_data(
            selections=[
                "time", "PODOSE_gli",
                "[Cve_gli]", "[Cve_m1]", "[Cve_m2]",
                "Aurine_m1", "Aurine_m2",
                "Afeces_gli", "Afeces_m1", "Afeces_m2",
            ]
        )
        return {}

    def figures_mpl(self) -> Dict[str, FigureMPL]:
        """Creates and returns the dose dependency figure object."""
        # Calculate PK parameters
        pk_df = self._calculate_pk_data()

        # Load study data
        df_study = self._load_study_data()

        # Create figure
        fig = self._create_dose_dependency_figure(pk_df, df_study)

        return {"pk_parameters_dose": fig}

    def _calculate_pk_data(self) -> pd.DataFrame:
        """Calculate PK parameters from simulation results."""
        sim_key = "dose_scan"
        task_key = f"task_{sim_key}"

        xres = self.results[task_key]
        pk_df = calculate_glimepiride_pk(experiment=self, xres=xres)

        # Set doses for metabolites
        self._assign_metabolite_doses(pk_df)

        # Convert units for all compounds
        mr_values = self.Mr
        compound_mr = {
            "glimepiride": mr_values.gli.m,
            "M1": mr_values.m1.m,
            "M2": mr_values.m2.m,
        }

        for compound, mr in compound_mr.items():
            self._convert_units(pk_df, compound, mr)

        return pk_df

    def _assign_metabolite_doses(self, pk_df: pd.DataFrame) -> None:
        """Assign doses to metabolites based on glimepiride doses."""
        glimepiride_mask = pk_df["compound"] == "glimepiride"
        if not glimepiride_mask.any():
            return

        # Get sorted glimepiride doses
        glimepiride_doses = pk_df[glimepiride_mask].sort_values(by="dose")["dose"].values

        # Assign same doses to metabolites
        for metabolite in ["M1", "M2"]:
            metabolite_mask = pk_df["compound"] == metabolite
            if metabolite_mask.any() and metabolite_mask.sum() == len(glimepiride_doses):
                pk_df.loc[metabolite_mask, "dose"] = glimepiride_doses

    def _convert_units(self, pk_df: pd.DataFrame, compound: str, mr: float) -> None:
        """Convert units for a compound."""
        mask = pk_df["compound"] == compound
        if not mask.any():
            return

        conversions = {
            "tmax": {"factor": 1/60, "new_unit": "hr"},
            "thalf": {"factor": 1/60, "new_unit": "hr"},
            "cmax": {"factor": mr * 1000, "new_unit": "ng/ml"},
            "auc": {"factor": (mr * 1000) / 60, "new_unit": "ng*hr/ml"}
        }

        for column, info in conversions.items():
            if column in pk_df.columns:
                pk_df.loc[mask, column] = pk_df.loc[mask, column] * info["factor"]
                pk_df.loc[mask, f"{column}_unit"] = info["new_unit"]

    def _load_study_data(self) -> pd.DataFrame:
        """Load and prepare study data."""
        study_data_path = DATA_PATH_GLIMEPIRIDE / "dose_dependency_data.xlsx"
        df_study = pd.read_excel(study_data_path, sheet_name="Tab1", comment="#")

        # Convert numeric columns
        numeric_columns = [c for c in df_study.columns if c.endswith(('_mean', '_sd')) or c == 'Dose']
        for col in numeric_columns:
            if col in df_study.columns:
                df_study[col] = pd.to_numeric(df_study[col], errors='coerce')

        return df_study

    def _create_dose_dependency_figure(self, pk_df: pd.DataFrame, df_study: pd.DataFrame) -> FigureMPL:
        """Create the dose dependency figure."""
        config = self.plot_config

        # Create figure
        fig, axes = plt.subplots(
            nrows=len(config["compounds"]),
            ncols=len(config["pk_params"]),
            figsize=config["figure_size"],
            dpi=config["dpi"]
        )

        fig.suptitle('Glimepiride Pharmacokinetics - Dose Dependency', fontsize=20, fontweight='bold')

        # Data by compound
        df_by_compound = {c: pk_df[pk_df["compound"] == c] for c in config["compounds"]}

        # Set column titles
        for col, param_label in enumerate(config["param_labels"]):
            axes[0, col].set_title(param_label, fontsize=16, fontweight="bold")

        # Plot each compound and parameter
        for row, compound in enumerate(config["compounds"]):
            df_compound = df_by_compound.get(compound, pd.DataFrame())

            for col, param in enumerate(config["pk_params"]):
                ax = axes[row, col]
                # Add grid
                ax.grid(True, alpha=0.4, color='lightgrey', linestyle='--')
                # Plot data
                if not df_compound.empty:
                    self._plot_simulation_line(ax, df_compound, param)
                self._plot_study_data(ax, df_study, compound, param)
                self._configure_axis(ax, row, col, compound, param, df_compound)

        # Add legend
        self._add_legend(fig)
        plt.tight_layout(rect=[0, 0.05, 1, 0.96])
        return fig

    def _configure_axis(self, ax, row: int, col: int, compound: str, param: str, df_compound: pd.DataFrame) -> None:
        """Configure axis settings."""
        # Get unit
        unit = ""
        if not df_compound.empty and f"{param}_unit" in df_compound.columns:
            valid_units = df_compound[f"{param}_unit"].dropna()
            if not valid_units.empty:
                unit = valid_units.iloc[0]

        # Labels
        if row == len(self.plot_config["compounds"]) - 1:
            ax.set_xlabel("Glimepiride Dose [mg]", fontsize=14, fontweight="bold")

        param_label = self.plot_config["param_labels"][col]
        y_label = f"{compound} {param_label}"
        if unit:
            y_label += f" [{unit}]"
        ax.set_ylabel(y_label, fontsize=14, fontweight="bold")

        # Axis limits
        ax.set_xlim(0, self.dose_range.max() * 1.05)
        ax.set_ylim(bottom=0)

        # Tick parameters
        ax.tick_params(axis='both', which='major', labelsize=10)

    def _plot_simulation_line(self, ax, df_compound: pd.DataFrame, param: str) -> None:
        """Plot simulation line."""
        if 'dose' in df_compound.columns and param in df_compound.columns:
            df_sorted = df_compound.sort_values(by="dose")
            ax.plot(df_sorted["dose"], df_sorted[param], color='black', linewidth=2)

    def _plot_study_data(self, ax, df_study: pd.DataFrame, compound: str, param: str) -> None:
        """Plot study data points."""
        # Check if dose_colors is defined
        if not hasattr(self, 'dose_colors'):
            return

        compound_prefixes = {"glimepiride": "gli", "M1": "m1", "M2": "m2"}
        if compound not in compound_prefixes:
            return

        prefix = compound_prefixes[compound]
        mean_col = f"{prefix}_{param}_mean"
        sd_col = f"{prefix}_{param}_sd"

        if mean_col not in df_study.columns or 'Dose' not in df_study.columns:
            return

        valid_study_data = df_study.dropna(subset=[mean_col, 'Dose'])

        for _, row in valid_study_data.iterrows():
            color = self.dose_colors.get(row["Dose"], 'gray')
            yerr = row[sd_col] if sd_col in row and pd.notna(row[sd_col]) else None

            ax.errorbar(
                x=row["Dose"],
                y=row[mean_col],
                yerr=yerr,
                fmt='s',
                markersize=8,
                markerfacecolor=color,
                markeredgecolor='black',
                ecolor='black',
                capsize=3,
                linestyle='None'
            )

    def _add_legend(self, fig) -> None:
        """Add legend to the figure."""
        legend_elements = [
            Line2D([0], [0], color='black', lw=2, label='Simulation'),
            Line2D([0], [0], marker='s', color='w', markerfacecolor='w',
                   markeredgecolor='black', markersize=8, label='Study Data', linestyle='None')
        ]

        fig.legend(
            handles=legend_elements,
            loc='lower center',
            ncol=2,
            bbox_to_anchor=(0.5, 0.02),
            frameon=True,
            fontsize=15
        )


if __name__ == "__main__":
    output_dir_name = "GlimepirideDoseScan"
    run_experiments(
        experiment_classes=GlimepirideDoseScan,
        output_dir=output_dir_name,
        save_results=True
    )