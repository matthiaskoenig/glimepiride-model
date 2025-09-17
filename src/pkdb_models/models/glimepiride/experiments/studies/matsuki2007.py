from typing import Dict
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console
from pkdb_models.models.glimepiride.experiments.base_experiment import GlimepirideSimulationExperiment
from pkdb_models.models.glimepiride.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, GlimepirideMappingMetaData, Coadministration, Genotype
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.glimepiride.helpers import run_experiments
from pathlib import Path
import pkdb_models.models.glimepiride as glimepiride

class Matsuki2007(GlimepirideSimulationExperiment):
    """Simulation experiment of Matsuki2007."""

    info = {
        "gli": "glimepiride"
    }

    interventions = [
        "GLI2",
        "GLI1"
    ]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                if label in ["glimepiride_GLI2", "glimepiride_GLI1"]:
                    dset = DataSet.from_df(df_label, self.ureg)
                    dset.unit_conversion("mean", 1 / self.Mr.gli)
                    dsets[label] = dset
        return dsets


    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        # Single dose simulation
        tcsims[f"po_gli2"] = TimecourseSim(
            [Timecourse(
                start=0,
                end=25 * 60,
                steps=500,
                changes={
                    **self.default_changes(),
                    "PODOSE_gli": Q_(2, "mg"),
                },
            )]
        )
        # Multiple-dose simulation
        tc0 = Timecourse(
            start=0,
            end=12 * 60,
            steps=500,
            changes={
                **self.default_changes(),
                "PODOSE_gli": Q_(1, "mg"),
            },
        )
        tc1 = Timecourse(
            start=0,
            end=12 * 60,
            steps=500,
            changes={
                "PODOSE_gli": Q_(1, "mg"),
            },
        )
        tcsims["po_gli1"] = TimecourseSim([tc0] + [tc1])
        return tcsims


    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for sid in self.info:
            name = self.info[sid]
            for intervention in self.interventions:
                mappings[f"fm_po_{intervention.lower()}_{sid}"] = FitMapping(
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
                        task=f"task_po_{intervention.lower()}",
                        xid="time",
                        yid=f"[Cve_{sid}]",
                    ),
                    metadata=GlimepirideMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE if intervention == "GLI1" else Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FED, # 2 mg of glimepiride before breakfast or twice-daily with 1 mg before breakfast and dinner
                        coadministration=Coadministration.NONE,
                        genotype=Genotype.NR,
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
        plots[0].xaxis.max = 24

        for sid in self.info:
            name = self.info[sid]
            for intervention in self.interventions:
                # Simulation
                plots[0].add_data(
                    task=f"task_po_{intervention.lower()}",
                    xid="time",
                    yid=f"[Cve_{sid}]",
                    label=f"Sim 2mg (single) PO" if intervention == "GLI2" else f"Sim 1mg (multiple) PO",
                    color="black" if intervention == "GLI2" else "tab:blue",
                )
                # Data
                plots[0].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"Data 2mg (single) PO" if intervention == "GLI2" else f"Data 1mg (multiple) PO",
                    color="black" if intervention == "GLI2" else "tab:blue",
                )
        return {fig.sid: fig}


if __name__ == "__main__":
    if __name__ == "__main__":
        out = glimepiride.RESULTS_PATH_SIMULATION / Matsuki2007.__name__
        out.mkdir(parents=True, exist_ok=True)
        run_experiments(Matsuki2007, output_dir=out)