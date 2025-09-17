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

class Ratheiser1993(GlimepirideSimulationExperiment):
    """Simulation experiment of Ratheiser1993."""
    # BMI: 22.5, bodyweight not reported
    doses = [0.25, 0.5, 0.75, 1, 1.25, 1.5]
    info = {
        "gli": "glimepiride",
        "Aurine_m1_m2": "glimepiride-M1+glimepiride-M2_urine",
    }
    interventions = [
        "GLI0.25",
        "GLI0.5",
        "GLI0.75",
        "GLI1",
        "GLI1.25",
        "GLI1.5",
    ]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Tab1B"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion
                if label.startswith("glimepiride_"):
                    dset.unit_conversion("mean", 1 / self.Mr.gli)
                dsets[label] = dset
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for dose in self.doses:
            tcsims[f"iv_gli_{dose}"] = TimecourseSim(
                [
                    Timecourse(
                        start=0,
                        end=1,
                        steps=20,
                        changes={
                            **self.default_changes(),
                            #"BW": Q_(self.bodyweight, "kg"), # bodyweight not reported
                            f"Ri_gli": Q_(dose, "mg/min"),
                        },
                    ),
                    Timecourse(
                        start=0,
                        end=52 * 60,
                        steps=2000,
                        changes={
                            f"Ri_gli": Q_(0, "mg/min"),
                        },
                    ),
                ]
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
                elif sid == "Aurine_m1_m2":
                    observable_yid = "Aurine_m1_m2"
                    tissue = Tissue.URINE
                    dataset_label = f"{name}_GLI{dose}"
                else:
                    continue

                mappings[f"fm_iv_gli_{dose}_{sid}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=dataset_label,
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd" if sid == "Aurine_m1_m2" else None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_iv_gli_{dose}",
                        xid="time",
                        yid=observable_yid,
                    ),
                    metadata=GlimepirideMappingMetaData(
                        tissue=tissue,
                        route=Route.IV,
                        application_form=ApplicationForm.SOLUTION,
                        dosing=Dosing.CONSTANT_INFUSION,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,  # overnight fasting period of 10 - 12 hours
                        coadministration=Coadministration.NONE,
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        Figure.legend_fontsize=7.5
        fig = Figure(
            experiment=self,
            sid="Fig1_Tab1B",
            num_rows=1,
            num_cols=2,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        plots[0].set_yaxis(self.label_gli_serum, unit=self.unit_gli)
        plots[1].set_yaxis(self.label_m1_m2_urine, unit=self.unit_m1_m2_urine)

        plots[0].xaxis.min = -0.5
        plots[0].xaxis.max = 10

        for dose in self.doses:
            # Plot 0: Serum concentrations
            plots[0].add_data(
                task=f"task_iv_gli_{dose}",
                xid="time",
                yid=f"[Cve_gli]",
                label=f"Sim gli {dose} mg IV",
                color=self.dose_colors[dose],
            )
            plots[0].add_data(
                dataset=f"glimepiride_GLI{dose}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"Data {dose} mg IV",
                color=self.dose_colors[dose],
                linestyle="",
            )

            # Plot 1: Urine amounts
            plots[1].add_data(
                task=f"task_iv_gli_{dose}",
                xid="time",
                yid="Aurine_m1_m2",
                label=f"Sim M1+M2 {dose} mg IV",
                color=self.dose_colors[dose],
            )
            plots[1].add_data(
                dataset=f"glimepiride-M1+glimepiride-M2_urine_GLI{dose}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"Data M1+M2 {dose} mg IV",
                color=self.dose_colors[dose],
                linestyle="",
            )

        return {fig.sid: fig}


if __name__ == "__main__":
    out = glimepiride.RESULTS_PATH_SIMULATION / Ratheiser1993.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Ratheiser1993, output_dir=out)