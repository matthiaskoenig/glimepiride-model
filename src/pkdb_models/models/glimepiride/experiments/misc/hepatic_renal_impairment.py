from typing import Dict
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.glimepiride.experiments.base_experiment import (
    GlimepirideSimulationExperiment,
)
from pkdb_models.models.glimepiride.helpers import run_experiments


class HepaticRenalImpairmentExperiment(GlimepirideSimulationExperiment):
    """Tests hepatic and renal impairment effects on glimepiride."""

    maps = {
        "hepatic": GlimepirideSimulationExperiment.cirrhosis_map,
        "renal": GlimepirideSimulationExperiment.renal_map,
    }
    parameters = {
        "hepatic": "f_cirrhosis",
        "renal": "KI__f_renal_function",
    }
    colors = {
        "hepatic": GlimepirideSimulationExperiment.cirrhosis_colors,
        "renal": GlimepirideSimulationExperiment.renal_colors,
    }

    impairments = list(maps.keys())

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for impairment in self.impairments:
            map = self.maps[impairment]
            parameter = self.parameters[impairment]

            for group, value in map.items():
                tcsims[f"gli_{impairment}_{group}"] = TimecourseSim(
                    Timecourse(
                        start=0,
                        end=36 * 60,
                        steps=500,
                        changes={
                            **self.default_changes(),
                            f"PODOSE_gli": Q_(3, "mg"),
                            parameter: Q_(value, "dimensionless")
                        },
                    )
                )
        return tcsims

    def figures(self) -> Dict[str, Figure]:
        figures = {}

        for impairment in self.impairments:
            map = self.maps[impairment]
            colors = self.colors[impairment]

            fig = Figure(
                experiment=self,
                sid=f"timecourse_po_{impairment}",
                num_rows=3,
                num_cols=3,
                name=f"Glimepiride - {impairment.title()} Impairment"
            )
            plots = fig.create_plots(xaxis=Axis(label=self.label_time, unit=self.unit_time), legend=True)

            # Species to plot in each subplot
            species = [
                # plasma
                "[Cve_gli]",
                "[Cve_m1]",
                "[Cve_m2]",

                # urine
                "Aurine_m1_m2",
                "Aurine_m1",
                "Aurine_m2",


                # feces
                "Afeces_m1_m2",
                "Afeces_m1",
                "Afeces_m2",

            ]

            for ksid, sid in enumerate(species):
                plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

                for group, value in map.items():
                    plots[ksid].add_data(
                        task=f"task_gli_{impairment}_{group}",
                        xid="time",
                        yid=sid,
                        label=f"{group}",
                        color=colors[group],
                    )
            figures[fig.sid] = fig
        return figures


if __name__ == "__main__":
    run_experiments(HepaticRenalImpairmentExperiment, output_dir=HepaticRenalImpairmentExperiment.__name__)