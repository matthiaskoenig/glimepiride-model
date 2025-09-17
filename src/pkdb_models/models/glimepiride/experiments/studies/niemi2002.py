from typing import Dict
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from pkdb_models.models.glimepiride.experiments.base_experiment import GlimepirideSimulationExperiment
from pkdb_models.models.glimepiride.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health,
    Fasting, GlimepirideMappingMetaData, Genotype
)
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.glimepiride.helpers import run_experiments
from pathlib import Path
import pkdb_models.models.glimepiride as glimepiride


class Niemi2002(GlimepirideSimulationExperiment):
    """Simulation experiment for Niemi2002."""

    bodyweight = 68.0  # "mean weight of 68 kg (range, 55-83 kg)"

    info = {
        "gli": "glimepiride"
    }

    interventions = [
        "GLI"
    ]

    groups = {
        "group1": "*1/*1",
        "group2": "*1/*2",
        "group3": "*1/*3",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                if label.startswith("glimepiride_"):
                    dset.unit_conversion("mean", 1 / self.Mr.gli)
                dsets[f"{label}"] = dset
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for group, genotype in self.groups.items():
            tcsims[f"po_gli_{group}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=13 * 60,
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweight, "kg"),
                        "PODOSE_gli": Q_(0.5, "mg"),
                        "LI__f_cyp2c9": Q_(self.cyp2c9_activity[genotype], "dimensionless")
                    },
                )]
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for group, genotype in self.groups.items():
            dataset_label = f"glimepiride_GLI_{group}"
            mappings[f"fm_po_gli_{group}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=dataset_label,
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=f"task_po_gli_{group}",
                    xid="time",
                    yid=f"[Cve_gli]",
                ),
                metadata=GlimepirideMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FED, # overnight fast; breakfast 15 min after administration; meal after 3 hours
                    genotype=Genotype(genotype),
                ),
            )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            num_rows=1,
            name=f"{self.__class__.__name__} (CYP2C9)",
        )
        plot = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)[0]
        plot.set_yaxis(self.label_gli_plasma, unit=self.unit_gli)

        for group, genotype in self.groups.items():
            plot.add_data(
                task=f"task_po_gli_{group}",
                xid="time",
                yid=f"[Cve_gli]",
                label=f"Sim 0.5mg PO ({genotype})",
                color=self.cyp2c9_colors[genotype],
            )
            plot.add_data(
                dataset=f"glimepiride_GLI_{group}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"Data 0.5mg PO ({genotype})",
                color=self.cyp2c9_colors[genotype],
            )
        return {fig.sid: fig}


if __name__ == "__main__":
    if __name__ == "__main__":
        out = glimepiride.RESULTS_PATH_SIMULATION / Niemi2002.__name__
        out.mkdir(parents=True, exist_ok=True)
        run_experiments(Niemi2002, output_dir=out)