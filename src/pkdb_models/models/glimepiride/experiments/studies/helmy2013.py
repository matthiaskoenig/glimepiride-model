from typing import Dict
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console
from pkdb_models.models.glimepiride.experiments.base_experiment import GlimepirideSimulationExperiment
from pkdb_models.models.glimepiride.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health, Fasting,
    GlimepirideMappingMetaData, Coadministration
)
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.glimepiride.helpers import run_experiments
from pathlib import Path
import pkdb_models.models.glimepiride.helpers as glimepiride


class Helmy2013(GlimepirideSimulationExperiment):
    """Simulation experiment of Helmy2013."""

    bodyweight = 72.75  # "weight range (60 - 80 kg) with a mean value of 72.75 +/- 7.75 kg"
    doses = [1, 2, 3, 4, 6]
    info = {"gli": "glimepiride"}
    interventions = ["GLI1", "GLI2", "GLI3", "GLI4", "GLI6"]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                if label.startswith("glimepiride_"):
                    dset.unit_conversion("mean", 1 / self.Mr.gli)
                dsets[label] = dset
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for dose in self.doses:
            tcsims[f"po_gli{dose}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=25 * 60,
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweight, "kg"),
                        f"PODOSE_gli": Q_(dose, "mg"),
                    },
                )]
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for dose in self.doses:
                observable_yid = f"[Cve_{sid}]"
                tissue = Tissue.PLASMA
                dataset_label = f"{name}_GLI{dose}"
                mappings[f"fm_po_gli{dose}_{sid}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=dataset_label,
                        xid="time",
                        yid="mean",
                        yid_sd=None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_po_gli{dose}",
                        xid="time",
                        yid=observable_yid,
                    ),
                    metadata=GlimepirideMappingMetaData(
                        tissue=tissue,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FED, # 12-hour overnight fast; meals were served at 2, 5, and 10 hours after drug dosing
                        coadministration=Coadministration.NONE,
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=1,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_gli_plasma, unit=self.unit_gli)
        for dose in self.doses:
            # Simulation data
            plots[0].add_data(
                task=f"task_po_gli{dose}",
                xid="time",
                yid=f"[Cve_gli]",
                label=f"Sim {dose}mg PO",
                color=self.dose_colors[dose]
            )
            # Experimental data
            plots[0].add_data(
                dataset=f"glimepiride_GLI{dose}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"Data {dose}mg PO",
                color=self.dose_colors[dose]
             )
        return {fig.sid: fig}


if __name__ == "__main__":
    out = glimepiride.RESULTS_PATH_SIMULATION / Helmy2013.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Helmy2013, output_dir=out)