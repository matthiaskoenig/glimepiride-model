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

class Lehr1990(GlimepirideSimulationExperiment):
    """Simulation experiment for Lehr1990."""
    info = {
        "gli": "glimepiride",
        "m1": "glimepiride-M1",
        "m2": "glimepiride-M2",
        "m1_m2": "M1+M2"
    }
    interventions = ["GLI3"]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        # Map label prefixes to Mr
        conversion_map = {
            "glimepiride_": self.Mr.gli,
            "glimepiride-M1": self.Mr.m1,
            "glimepiride-M2": self.Mr.m2,
        }
        for fig_id in ["Fig4"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                for prefix, mr in conversion_map.items():
                    if label.startswith(prefix):
                        dset.unit_conversion("mean", 1 / mr)
                        break
                dsets[label] = dset
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        tcsims[f"po_gli3"] = TimecourseSim(
            [Timecourse(
                start=0,
                end=15 * 60,
                steps=500,
                changes={
                    **self.default_changes(),
                    "PODOSE_gli": Q_(3, "mg"),
                },
            )
        ])
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for sid in self.info:
            name = self.info[sid]
            for intervention in self.interventions:
                mappings[f"fm_po_gli3_{sid}_{intervention}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_po_gli3",
                        xid="time",
                        yid=f"[Cve_{sid}]",
                    ),
                    metadata=GlimepirideMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.NR,
                        coadministration=Coadministration.NONE,
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid=f"Fig4",
            num_rows=1,
            num_cols=4,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_gli_plasma, unit=self.unit_gli)
        plots[1].set_yaxis(self.label_m1_plasma, unit=self.unit_m1)
        plots[2].set_yaxis(self.label_m2_plasma, unit=self.unit_m2)
        plots[3].set_yaxis(self.label_m1_m2_plasma, unit=self.unit_m1_m2)

        for intervention in self.interventions:
            for k, sid in enumerate(self.info):
                name = self.info[sid]
                # Simulation data
                plots[k].add_data(
                    task=f"task_po_gli3",
                    xid="time",
                    yid=f"[Cve_{sid}]",
                    label=f"Sim 3mg PO",
                    color="black"
                )
                # Experimental data
                plots[k].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    count="count",
                    label=f"Data 3mg PO",
                    color="black"
                )
        return {fig.sid: fig}


if __name__ == "__main__":
    if __name__ == "__main__":
        out = glimepiride.RESULTS_PATH_SIMULATION / Lehr1990.__name__
        out.mkdir(parents=True, exist_ok=True)
        run_experiments(Lehr1990, output_dir=out)