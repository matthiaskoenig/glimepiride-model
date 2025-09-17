from typing import Dict
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.glimepiride.experiments.base_experiment import GlimepirideSimulationExperiment
from pkdb_models.models.glimepiride.helpers import run_experiments


class DoseDependencyExperiment(GlimepirideSimulationExperiment):
    """Dose-response simulation for glimepiride administration."""

    # routes = {"gli": ["PO", "IV"], "m1": ["IV"]}
    routes = {"gli": ["PO"]}

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for substance, routes in self.routes.items():
            for route in routes:
                for dose in self.dose_values:
                    tcsims[f"{substance}_{route}_{dose}"] = TimecourseSim(
                        Timecourse(
                            start=0,
                            end=24 * 60,  # [min]
                            steps=500,
                            changes={
                                **self.default_changes(),
                                f"{route}DOSE_{substance}": Q_(dose, "mg"),
                            },
                        )
                    )
        return tcsims

    def figures(self) -> Dict[str, Figure]:
        Figure.legend_fontsize = 10
        figures = {}
        for substance, routes in self.routes.items():
            for route in routes:
                fig = Figure(
                    experiment=self,
                    sid=f"timecourse_po_dose",
                    num_rows=3,
                    num_cols=3,
                    name=f"Glimepiride - Dose Dependency",
                )

                plots = fig.create_plots(xaxis=Axis(label=self.label_time, unit=self.unit_time), legend=True)

                # Define plot contents (plasma, urine, feces)
                sids = [
                    "[Cve_gli]", "[Cve_m1]", "[Cve_m2]",  # plasma
                    "Aurine_m1_m2", "Aurine_m1", "Aurine_m2",  # urine
                    "Afeces_m1_m2", "Afeces_m1", "Afeces_m2",  # feces
                ]

                # Set axis labels and units
                for ksid, sid in enumerate(sids):
                    plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

                # Add data for each dose
                for ksid, sid in enumerate(sids):
                    for kval, dose in enumerate(self.dose_values):
                        plots[ksid].add_data(
                            task=f"task_{substance}_{route}_{dose}",
                            xid="time",
                            yid=sid,
                            label=f"{dose} mg",
                            color=self.dose_colors[dose],
                        )
                figures[fig.sid] = fig
        return figures


if __name__ == "__main__":
    run_experiments(DoseDependencyExperiment, output_dir=DoseDependencyExperiment.__name__)