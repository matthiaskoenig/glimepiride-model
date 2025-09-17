from typing import Dict
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console
from pkdb_models.models.glimepiride.experiments.base_experiment import GlimepirideSimulationExperiment
from pkdb_models.models.glimepiride.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, GlimepirideMappingMetaData, Coadministration, Genotype
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pathlib import Path
import pkdb_models.models.glimepiride as glimepiride
from pkdb_models.models.glimepiride.helpers import run_experiments


class Ahmed2016(GlimepirideSimulationExperiment):
    """Simulation experiment of Ahmed2016."""

    bodyweight = 75.0  # [kg] "75 ± 7.5 kg median weight"

    info = {"gli": "glimepiride"}

    interventions = ["GLI3"] # 1mg PO

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig5"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                if label == "glimepiride_GLI3":
                    dset = DataSet.from_df(df_label, self.ureg)
                    # Convert the 'mean' values
                    dset.unit_conversion("mean", 1 / self.Mr.gli)
                    dsets[label] = dset
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        tcsims[f"po_gli3"] = TimecourseSim(
            [Timecourse(
                start=0,
                end=25 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    "BW": Q_(self.bodyweight, "kg"),
                    "PODOSE_gli": Q_(1, "mg"),
                },
            )]
        )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for sid in self.info:
            name = self.info[sid]
            for intervention in self.interventions:
                mappings[f"fm_po_gli3_{sid}_{intervention}"] = FitMapping(
                    self,
                    # Reference experimental data (dataset)
                    reference=FitData(
                        self,
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        count="count",
                    ),
                    # Simulation observable (model output)
                    observable=FitData(
                        self,
                        task="task_po_gli3",
                        xid="time",
                        yid=f"[Cve_{sid}]",
                    ),
                    # Metadata describing the mapping
                    metadata=GlimepirideMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.NR,
                        coadministration=Coadministration.NONE,
                        genotype=Genotype.NR,
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig5",
            num_rows=1,
            name=f"{self.__class__.__name__} (Healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_gli_plasma, unit=self.unit_gli)

        for sid in self.info:
            name = self.info[sid]
            for intervention in self.interventions:
                # simulation
                plots[0].add_data(
                    task=f"task_po_gli3",
                    xid="time",
                    yid=f"[Cve_{sid}]",
                    label=f"Sim 1mg PO",
                    color="black"
                )
                # data
                plots[0].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"Data 1mg PO",
                    color="black"
                )
        return {
            fig.sid: fig,
        }

if __name__ == "__main__":
    out = glimepiride.RESULTS_PATH_SIMULATION / Ahmed2016.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Ahmed2016, output_dir=out)