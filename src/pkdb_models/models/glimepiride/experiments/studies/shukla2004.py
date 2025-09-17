from typing import Dict
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from pkdb_models.models.glimepiride.experiments.base_experiment import GlimepirideSimulationExperiment
from pkdb_models.models.glimepiride.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health,
    Fasting, GlimepirideMappingMetaData, Coadministration
)
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.glimepiride.helpers import run_experiments
from sbmlutils.console import console
from pathlib import  Path
import pkdb_models.models.glimepiride as glimepiride

class Shukla2004(GlimepirideSimulationExperiment):
    """Simulation experiment of Shukla2004."""
    bodyweights = {
        "normal": 72.0,  # [kg]
        "obese": 130.0,  # [kg]
    }

    groups = list(bodyweights.keys())

    info = {
        "[Cve_gli]": "glimepiride",
        "[Cve_m1]": "glimepiride-M1",
        "[Cve_m2]": "glimepiride-M2",
        "[Cve_m1_m2]": "M1+M2",
        "Aurine_m1": "glimepiride-M1_amount",
        "Aurine_m2": "glimepiride-M2_amount",
        "Aurine_m1_m2": "M1+M2_amount",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                # concentrations
                if label.startswith("glimepiride-") and not "amount" in label:
                    if "M1" in label:
                        dset.unit_conversion("mean", 1 / self.Mr.m1)
                    elif "M2" in label:
                        dset.unit_conversion("mean", 1 / self.Mr.m2)
                elif label.startswith("glimepiride") and not "amount" in label:
                    dset.unit_conversion("mean", 1 / self.Mr.gli)
                # amounts
                elif "amount" in label:
                    dset.unit_conversion("mean", 1 / self.Mr.m1_m2)
                dsets[label] = dset
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        # simulation for each bodyweight group
        for group in ["normal", "obese"]:
            tcsims[f"GLI8_{group}"] = TimecourseSim(
                [
                    Timecourse(
                        start=0,
                        end=25 * 60,
                        steps=500,
                        changes={
                            **self.default_changes(),
                            "BW": Q_(self.bodyweights[group], "kg"),
                            "PODOSE_gli": Q_(8.0, "mg"),
                        },
                    )
                ])
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for group in self.groups:
            for sid, name in self.info.items():
                mappings[f"fm_{name}_{group}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{group}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd" if sid.startswith("Aurine") else None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_GLI8_{group}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=GlimepirideMappingMetaData(
                        tissue=Tissue.URINE if sid.startswith("Aurine") else Tissue.SERUM,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY if group == "normal" else Health.OBESE,
                        fasting=Fasting.FED, # overnight fast; from lunch on no food restriction
                        coadministration=Coadministration.NONE,
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_plasma(),
            **self.figure_urine(),
        }

    def figure_plasma(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=1,
            num_cols=4,
            name=f"{self.__class__.__name__} (Plasma, Bodyweight)",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_gli_plasma, self.unit_gli)
        plots[1].set_yaxis(self.label_m1_plasma, self.unit_m1)
        plots[2].set_yaxis(self.label_m2_plasma, self.unit_m2)
        plots[3].set_yaxis(self.label_m1_m2_plasma, self.unit_m1_m2)

        for k, sid in enumerate(self.info):
            if sid.startswith("Aurine"):  # Skip urine data
                continue
            substance = self.info[sid]
            for group in self.groups:
                plots[k].add_data(
                    task=f"task_GLI8_{group}",
                    xid="time",
                    yid=sid,
                    label=f"Sim 8mg gli po ({group})",
                    color=self.colors()[group],
                )
                plots[k].add_data(
                    dataset=f"{substance}_{group}",
                    xid="time",
                    yid="mean",
                    count="count",
                    label=f"Data 8mg gli po ({group})",
                    color=self.colors()[group],
                )
        return {fig.sid: fig}

    def figure_urine(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            num_rows=1,
            num_cols=3,
            name=f"{self.__class__.__name__} (Urine, Bodyweight)",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_m1_urine, self.unit_m1_urine)
        plots[1].set_yaxis(self.label_m2_urine, self.unit_m2_urine)
        plots[2].set_yaxis(self.label_m1_m2_urine, self.unit_m1_m2_urine)

        k = 0
        for sid in self.info:
            if not sid.startswith("Aurine"):  # Only process urine data
                continue
            substance = self.info[sid]

            for group in self.groups:
                plots[k].add_data(
                    task=f"task_GLI8_{group}",
                    xid="time",
                    yid=sid,
                    label=f"Sim 8mg PO ({group})",
                    color=self.colors()[group],
                )
                plots[k].add_data(
                    dataset=f"{substance}_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"Data 8mg PO ({group})",
                    color=self.colors()[group],
                )
            k += 1
        return {fig.sid: fig}

    def colors(self):
        """Get same colors consistent with BodyweightScanExperiment."""
        # Colors for custom color palette
        colors = ['#F5C6D0', '#BB4681', '#601038']
        cmap = LinearSegmentedColormap.from_list("custom_bw_cmap", colors)

        # Define bodyweight range (same as in BodyweightScanExperiment)
        bodyweight_range = np.linspace(45, 170, 100)

        # Normalized positions
        bw_min, bw_max = min(bodyweight_range), max(bodyweight_range)
        normal_pos = (self.bodyweights["normal"] - bw_min) / (bw_max - bw_min)
        obese_pos = (self.bodyweights["obese"] - bw_min) / (bw_max - bw_min)

        # Get colors at positions
        return {
            "normal": cmap(normal_pos),
            "obese": cmap(obese_pos),
        }

if __name__ == "__main__":
    out = glimepiride.RESULTS_PATH_SIMULATION / Shukla2004.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Shukla2004, output_dir=out)