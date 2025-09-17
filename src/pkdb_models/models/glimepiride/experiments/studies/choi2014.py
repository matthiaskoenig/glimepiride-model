from typing import Dict
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from pkdb_models.models.glimepiride.experiments.base_experiment import GlimepirideSimulationExperiment
from pkdb_models.models.glimepiride.experiments.metadata import (
    Tissue,
    Route,
    Dosing,
    ApplicationForm,
    Health,
    Fasting,
    GlimepirideMappingMetaData,
    Coadministration,
    Genotype,
)
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.glimepiride.helpers import run_experiments
import pkdb_models.models.glimepiride as glimepiride
from pathlib import Path


class Choi2014(GlimepirideSimulationExperiment):
    """Simulation experiment of Choi2014."""

    bodyweight = 70.8  # "mean (SD) weight was 70.8 (7.8) kg"

    info = {
        "gli": "glimepiride",
        "m1": "glimepiride-M1",
    }
    interventions = [
        "GLI4",
        "GLI4, GEM50"
    ]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion
                if label.startswith("glimepiride_"):
                    dset.unit_conversion("mean", 1 / self.Mr.gli)
                elif label.startswith("glimepiride-M1"):
                    dset.unit_conversion("mean", 1 / self.Mr.m1)
                dsets[f"{label}"] = dset
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        tcsims[f"po_gli4"] = TimecourseSim(
            [Timecourse(
                start=0,
                end=25 * 60,
                steps=500,
                changes={
                    **self.default_changes(),
                    "BW": Q_(self.bodyweight, "kg"),
                    "PODOSE_gli": Q_(4, "mg"),
                },
            )]
        )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for k, sid in enumerate(self.info):
            name = self.info[sid]
            for intervention in self.interventions:
                mappings[f"fm_po_gli4_{sid}_{intervention}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_po_gli4", xid="time", yid=f"[Cve_{sid}]",
                    ),
                    metadata=GlimepirideMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FED, # overnight fast, food restriction until 1h after administration
                        coadministration=Coadministration.GEMIGLIPTIN if "GEM50" in intervention else Coadministration.NONE,
                        genotype=Genotype.NR,
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy, Coadmin)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_gli_plasma, unit=self.unit_gli)
        plots[1].set_yaxis(self.label_m1_plasma, unit=self.unit_m1)

        for k, sid in enumerate(self.info):
            # simulation
            plots[k].add_data(
            task=f"task_po_gli4",
            xid="time",
            yid=f"[Cve_{sid}]",
            label=f"Sim 4mg GLI PO",
            color="black"
        )

        for k, sid in enumerate(self.info):
            name = self.info[sid]
            for intervention in self.interventions:
                # study data
                plots[k].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"Data 4mg GLI PO" if intervention == "GLI4" else f"Data 4mg GLI & 50mg GEM PO",
                    color="tab:blue" if "GEM50" in intervention else "black"
                )
        return {
            fig.sid: fig,
        }

if __name__ == "__main__":
    out = glimepiride.RESULTS_PATH_SIMULATION / Choi2014.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Choi2014, output_dir=out)