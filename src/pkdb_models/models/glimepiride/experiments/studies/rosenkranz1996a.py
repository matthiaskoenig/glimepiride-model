from typing import Dict
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from pkdb_models.models.glimepiride.experiments.base_experiment import GlimepirideSimulationExperiment
from pkdb_models.models.glimepiride.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health, Fasting, GlimepirideMappingMetaData, Coadministration
)
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.glimepiride.helpers import run_experiments
import pkdb_models.models.glimepiride as glimepiride
from pathlib import Path

class Rosenkranz1996a(GlimepirideSimulationExperiment):
    """Simulation experiment for Rosenkranz1996a."""

    group_colors = {
        "single_group1": "#66c2a4",
        "single_group2": "#2ca25f",
        "single_group3": "#006d2c",
        "single_group_normal": "#000000",
    }
    group_renal_functions = {
        "single_group1": "mild",
        "single_group2": "moderate",
        "single_group3": "severe",
        "single_group_normal": "normal",
    }
    renal_functions = {
        "normal": 1.0,
        "mild": 0.5,
        "moderate": 0.35,
        "severe": 0.2,
    }
    info = {
        "gli": "glimepiride",
        "m1": "glimepiride-M1",
        "m2": "glimepiride-M2",
        "m1_m2": "M1+M2",
        "m1_m2_urine": "M1+M2_urine"
    }
    interventions = [
        "GLI3"
    ]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Tab2B"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                if label.startswith("glimepiride_"):
                    dset.unit_conversion("value", 1 / self.Mr.gli)
                elif label.startswith("glimepiride-M1"):
                    dset.unit_conversion("value", 1 / self.Mr.m1)
                elif label.startswith("glimepiride-M2"):
                    dset.unit_conversion("value", 1 / self.Mr.m2)
                elif label.startswith("M1+M2_urine"):
                    dset.unit_conversion("mean", 1 / self.Mr.m1_m2)
                elif label.startswith("M1+M2"):
                    dset.unit_conversion("value", 1 / self.Mr.m1_m2)
                dsets[label] = dset
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        # simulation for each renal function group
        for group, renal_function_key in self.group_renal_functions.items():
            tcsims[f"po_gli3_{group}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=50 * 60,
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "PODOSE_gli": Q_(3, "mg"),
                        "KI__f_renal_function": Q_(self.renal_functions[renal_function_key], "dimensionless"),
                    },
                )]
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        # Urinary recovery M1+M2
        for group in self.group_renal_functions.keys():
            label = f"M1+M2_urine_GLI3_{group}"

            if label not in self.datasets():
                continue

            mappings[f"fm_po_gli3_m1_m2_urine_{group}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=label,
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=f"task_po_gli3_{group}",
                    xid="time",
                    yid="Aurine_m1_m2"
                ),
                metadata=GlimepirideMappingMetaData(
                    tissue=Tissue.URINE,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.T2DM_RENAL_IMPAIRMENT,
                    fasting=Fasting.FED,
                    coadministration=Coadministration.NONE,
                ),
            )

        # Plasma concentrations S1
        for sid, metabolite_name in self.info.items():
            if sid == "m1_m2_urine":
                continue
            if label not in self.datasets():
                continue
            mappings[f"fm_po_gli3_{sid}_plasma"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"{metabolite_name}_GLI3_S1",
                    xid="time",
                    yid="value",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task="task_po_gli3_single_group3",
                    xid="time",
                    yid=f"[Cve_{sid}]"
                ),
                metadata=GlimepirideMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.T2DM_RENAL_IMPAIRMENT,
                    fasting=Fasting.FED,
                    coadministration=Coadministration.NONE,
                ),
            )
        return mappings

    def figure_plasma_urine(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1_Tab2",
            num_rows=1,
            num_cols=4,
            name=f"{self.__class__.__name__} (Renal Impairment)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )

        # Set y-axes for the plots
        plots[0].set_yaxis(self.label_gli_plasma, self.unit_gli)
        plots[1].set_yaxis(self.label_m1_plasma, self.unit_m1)
        plots[2].set_yaxis(self.label_m2_plasma, self.unit_m2)
        plots[3].set_yaxis(self.label_m1_m2_urine, unit=self.unit_m1_m2_urine)

        # Set x-axis limits
        for k in range(len(plots)):
            plots[k].xaxis.min = 0

        plots[0].xaxis.max = 10
        plots[1].xaxis.max = 15
        plots[2].xaxis.max = 15
        plots[3].xaxis.max = 50

        # For plasma plots (gli, m1, m2)
        plasma_sids = ["gli", "m1", "m2"]

        for k, sid in enumerate(plasma_sids):
            # Add simulation normal renal function
            plots[k].add_data(
                task="task_po_gli3_single_group_normal",
                xid="time",
                yid=f"[Cve_{sid}]",
                label=f"Sim 3mg PO (normal)",
                color=self.group_colors["single_group_normal"]
            )

            # Add simulation and data for severe renal impairment
            name = self.info[sid]
            plots[k].add_data(
                task="task_po_gli3_single_group3",
                xid="time",
                yid=f"[Cve_{sid}]",
                label=f"Sim 3mg PO (severe)",
                color="#006d2c"
            )

            if f"{name}_GLI3_S1" in self.datasets():
                plots[k].add_data(
                    dataset=f"{name}_GLI3_S1",
                    xid="time",
                    yid="value",
                    count="count",
                    label=f"Data 3mg PO (severe)",
                    color="#006d2c"
                )

        # Add urinary data for 4th plot
        urinary_plot = plots[3]
        for group, color in self.group_colors.items():
            renal_label = self.group_renal_functions[group]
            urinary_plot.add_data(
                task=f"task_po_gli3_{group}",
                xid="time",
                yid="Aurine_m1_m2",
                label=f"Sim 3mg ({renal_label})",
                color=color,
            )

            # Add data
            dataset_name = f"M1+M2_urine_GLI3_{group}"
            if dataset_name in self.datasets():
                urinary_plot.add_data(
                    dataset=dataset_name,
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"Data 3mg ({renal_label})",
                    color=color,
                    linestyle="",
                )

        return {fig.sid: fig}

    def figure_plasma(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=1,
            num_cols=4,
            name=f"{self.__class__.__name__} (Renal Impairment, Plasma)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )

        plots[0].set_yaxis(self.label_gli_plasma, self.unit_gli)
        plots[1].set_yaxis(self.label_m1_plasma, self.unit_m1)
        plots[2].set_yaxis(self.label_m2_plasma, self.unit_m2)
        plots[3].set_yaxis(self.label_m1_m2_plasma, self.unit_m1_m2)

        for k in range(len(plots)):
            plots[k].xaxis.min = 0

        plots[0].xaxis.max = 10
        plots[1].xaxis.max = 15
        plots[3].xaxis.max = 15

        sids = [key for key in self.info.keys() if key != "m1_m2_urine"]

        # simulation line for normal renal function
        for k, sid in enumerate(sids):
            plots[k].add_data(
                task="task_po_gli3_single_group_normal",
                xid="time",
                yid=f"[Cve_{sid}]",
                label=f"Sim 3mg PO (normal)",
                color=self.group_colors["single_group_normal"]
            )

        for k, sid in enumerate(sids):
            name = self.info[sid]
            plots[k].add_data(
                task="task_po_gli3_single_group3",
                xid="time",
                yid=f"[Cve_{sid}]",
                label=f"Sim 3mg PO (severe)",
                color="#006d2c"
                )
            plots[k].add_data(
                dataset=f"{name}_GLI3_S1",
                xid="time",
                yid="value",
                count="count",
                label=f"Data 3mg PO (severe)",
                color="#006d2c"
            )
        return {fig.sid: fig}

    def figure_urine(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Tab2",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__} (Renal Impairment, Urine)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        urinary_plot = plots[0]
        urinary_plot.set_yaxis(self.label_m1_m2_urine, unit=self.unit_m1_m2_urine)

        for group, color in self.group_colors.items():
            renal_label = self.group_renal_functions[group]
            urinary_plot.add_data(
                task=f"task_po_gli3_{group}",
                xid="time",
                yid="Aurine_m1_m2",
                label=f"Sim 3mg ({renal_label})",
                color=color,
            )
            # Add data points if dataset exists
            dataset_name = f"M1+M2_urine_GLI3_{group}"
            if dataset_name in self.datasets():
                urinary_plot.add_data(
                    dataset=dataset_name,
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"Data 3mg ({renal_label})",
                    color=color,
                    linestyle="",
            )
        return {fig.sid: fig}

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_plasma(),
            **self.figure_urine(),
            **self.figure_plasma_urine()
        }


if __name__ == "__main__":
    out = glimepiride.RESULTS_PATH_SIMULATION / Rosenkranz1996a.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Rosenkranz1996a, output_dir=out)