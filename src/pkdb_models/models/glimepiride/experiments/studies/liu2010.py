from typing import Dict
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
import pkdb_models.models.glimepiride as glimepiride
from pathlib import Path

class Liu2010(GlimepirideSimulationExperiment):
    """
    Simulation experiment of Liu2010.
    food after 4h -> possible effect on pk
    """

    bodyweight = 64.0  # "64.0 [8.4] kg [range, 52.0–82.0 kg]"

    info = {
        "gli": "glimepiride",
        "m1": "glimepiride-M1",
    }
    interventions = [
        "GLI2test",
        "GLI2ref",
    ]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1A", "Fig1B"]:
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
        tcsims["po_GLI2"] = TimecourseSim([
            Timecourse(
                start=0,
                end=50 * 60,
                steps=500,
                changes={
                    **self.default_changes(),
                    "BW": Q_(self.bodyweight, "kg"),
                    "PODOSE_gli": Q_(2, "mg"),
                },
            )
        ])
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for intervention in self.interventions:
                mappings[f"fm_po_{intervention}_{sid}"] = FitMapping(
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
                        self,
                        task=f"task_po_GLI2",
                        xid="time",
                        yid=f"[Cve_{sid}]",
                    ),
                    metadata=GlimepirideMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FED, # 10-hour overnight fast, food intake was allowed 4 hours after treatment
                        coadministration=Coadministration.NONE,
                        outlier=True  # very strange peaks, probably due to food
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_gli_plasma, unit=self.unit_gli)
        plots[1].set_yaxis(self.label_m1_plasma, unit=self.unit_m1)

        for sid, name in self.info.items():
            plots[0 if sid == "gli" else 1].add_data(
                task=f"task_po_GLI2",
                xid="time",
                yid=f"[Cve_{sid}]",
                label="Sim 2mg PO",
                color="black"
            )

        for sid, name in self.info.items():
            for intervention in self.interventions:
                # data
                plots[0 if sid == "gli" else 1].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label="Data 2mg (test) PO" if "test" in intervention else "Data 2mg (ref) PO",
                    color="tab:blue" if "test" in intervention else "black"
                )
        return {fig.sid: fig}


if __name__ == "__main__":
    if __name__ == "__main__":
        out = glimepiride.RESULTS_PATH_SIMULATION / Liu2010.__name__
        out.mkdir(parents=True, exist_ok=True)
        run_experiments(Liu2010, output_dir=out)