from typing import Dict

from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.glimepiride.experiments.base_experiment import (
    GlimepirideSimulationExperiment,
)
from pkdb_models.models.glimepiride.helpers import run_experiments


class GeneticVariantExperiment(GlimepirideSimulationExperiment):
    """Simulation of CYP2C9 variant effects on glimepiride pharmacokinetics."""

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for variant, activity in self.cyp2c9_activity.items():
            tcsims[f"gli_{variant}"] = TimecourseSim(
                Timecourse(
                    start=0,
                    end=24 * 60,  # [min]
                    steps=2000,
                    changes={
                        **self.default_changes(),
                        "PODOSE_gli": Q_(4, "mg"),
                        "LI__f_cyp2c9": Q_(activity, "dimensionless"),
                    },
                )
            )
        return tcsims

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_pk(),
            #**self.figure_liver(),
        }

    def figure_pk(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="timecourse_po_genetic_variants",
            num_rows=3,
            num_cols=3,
            name="Glimepiride - CYP2C9 Variants",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, self.unit_time), legend=True)

        sids = [
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

        for ksid, sid in enumerate(sids):
            plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

            for variant in self.cyp2c9_activity:
                plots[ksid].add_data(
                    task=f"task_gli_{variant}",
                    xid="time",
                    yid=sid,
                    label=variant,
                    color=self.cyp2c9_colors[variant],
                )

        return {fig.sid: fig}

    def figure_liver(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig_genetic_variants_liver",
            num_rows=2,
            num_cols=3,
            name="Genetic Variants - Liver",
        )
        plots = fig.create_plots(xaxis=Axis("Time", unit="hr"), legend=True)

        sids = [
            # plasma & liver compartments
            "[Cve_gli]",
            "[Cve_m1]",
            "[Cve_m2]",
            "[LI__gli]",
            "[LI__m1]",
            "[LI__m2]",
        ]

        for ksid, sid in enumerate(sids):
            plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

            for variant in self.cyp2c9_activity:
                plots[ksid].add_data(
                    task=f"task_gli_{variant}",
                    xid="time",
                    yid=sid,
                    label=variant,
                    color=self.cyp2c9_colors[variant],
                )

        return {fig.sid: fig}


if __name__ == "__main__":
    run_experiments(GeneticVariantExperiment, output_dir=GeneticVariantExperiment.__name__)