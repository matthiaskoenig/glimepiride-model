from typing import Dict
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console
from pkdb_models.models.glimepiride.experiments.base_experiment import (
    GlimepirideSimulationExperiment,
)
from pkdb_models.models.glimepiride.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, GlimepirideMappingMetaData, Coadministration
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.glimepiride.helpers import run_experiments
from pathlib import Path
import pkdb_models.models.glimepiride as glimepiride

class Malerczyk1994(GlimepirideSimulationExperiment):
    """Simulation experiment of Malerczyk1994."""
    bodyweight = 78.0  # "mean weight 78 kg (64-98)"
    doses = [1, 2, 4, 8]
    info = {
        "gli": "glimepiride",
        "m1_m2": "M1+M2",
        "mean_m1_m2": "amount_mean_M1+M2",
    }
    interventions = [
        "GLI1",
        "GLI2",
        "GLI4",
        "GLI8"
    ]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2", "Fig6"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion
                if label.startswith("glimepiride_"):
                    dset.unit_conversion("mean", 1 / self.Mr.gli)
                if label.startswith("amount_mean_M1+M2"):
                    dset.unit_conversion("mean", 1 / self.Mr.m1_m2)
                dsets[label] = dset
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for dose in self.doses:
            tcsims[f"po_gli{dose}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=50 * 60,
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
        for k, sid in enumerate(self.info):
            name = self.info[sid]
            for dose in self.doses:
                if sid == "gli":
                    observable_yid = f"[Cve_{sid}]"
                    tissue = Tissue.PLASMA
                    dataset_label = f"{name}_GLI{dose}"
                elif sid == "m1_m2":
                    observable_yid = f"Aurine_{sid}"
                    tissue = Tissue.URINE
                    dataset_label = f"{name}_GLI{dose}"
                elif sid == "mean_m1_m2":
                    observable_yid = f"Aurine_m1_m2"
                    tissue = Tissue.URINE
                    dataset_label = f"{name}_GLI{dose}"
                else:
                    continue

                mappings[f"fm_po_gli{dose}_{sid}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=dataset_label,
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd" if sid == "mean_m1_m2" else None,
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
                        fasting=Fasting.FASTED, # overnight fasting period of 12 hours, and two hours after administration
                        coadministration=Coadministration.NONE,
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1_Fig2",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_gli_plasma, unit=self.unit_gli)
        plots[1].set_yaxis(self.label_m1_m2_urine, unit=self.unit_m1_m2_urine)

        plots[0].xaxis.max = 37

        for dose in self.doses:
            plots[0].add_data(
                task=f"task_po_gli{dose}",
                xid="time",
                yid=f"[Cve_gli]",
                label=f"Sim gli {dose} mg po",
                color=self.dose_colors[dose]
            )
            plots[0].add_data(
                dataset=f"glimepiride_GLI{dose}",
                xid="time",
                yid="mean",
                yid_sd=None,
                count="count",
                label=f"Data {dose}mg PO",
                color=self.dose_colors[dose]
             )

            plots[1].add_data(
                task=f"task_po_gli{dose}",
                xid="time",
                yid=f"Aurine_m1_m2",
                label=f"Sim {dose}mg PO",
                color=self.dose_colors[dose]
            )
            plots[1].add_data(
                dataset=f"M1+M2_GLI{dose}",
                xid="time",
                yid="mean",
                yid_sd=None,
                count="count",
                label=f"Data {dose}mg PO",
                color=self.dose_colors[dose]
            )
            plots[1].add_data(
                dataset=f"amount_mean_M1+M2_GLI{dose}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=None,
                color=self.dose_colors[dose]
             )
        return {fig.sid: fig,}


if __name__ == "__main__":
    if __name__ == "__main__":
        out = glimepiride.RESULTS_PATH_SIMULATION / Malerczyk1994.__name__
        out.mkdir(parents=True, exist_ok=True)
        run_experiments(Malerczyk1994, output_dir=out)