from typing import Dict
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from sbmlutils.console import console
from pkdb_models.models.glimepiride.experiments.base_experiment import GlimepirideSimulationExperiment
from pkdb_models.models.glimepiride.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health,
    Fasting, GlimepirideMappingMetaData, Coadministration
)
from pkdb_models.models.glimepiride.helpers import run_experiments
import pkdb_models.models.glimepiride as glimepiride
from pathlib import Path


class FDA1995(GlimepirideSimulationExperiment):
    """Simulation experiment of FDA1995."""

    info = {
        # Fig1 (Concentrations):
        "[Cve_gli]": "glimepiride",
        "[Cve_m1]": "glimepiride-M1",
        "[Cve_m2]": "glimepiride-M2",
        "[Cve_m1_m2]": "M1+M2",

        # Fig2 (Urine amounts):
        "Aurine_m1_m2": "urine_M1+M2",

        # Fig4 (Feces amounts, only PO):
        "Afeces_gli": "feces_glimepiride",
        "Afeces_m1": "feces_glimepiride-M1",
        "Afeces_m2": "feces_glimepiride-M2",
        "Afeces_m1_m2": "feces_M1+M2",
    }

    colors = {
        "iv": "tab:blue",
        "po": "black",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        # Map conversion keywords to Mr
        conversion_map = {
            "M1+M2": "m1_m2",
            "M1": "m1",
            "M2": "m2",
            "glimepiride": "gli",
        }
        for fig_id in ["Fig1", "Fig2", "Tab3A"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                if "urine" in label:
                    factor = 1 / self.Mr.m1_m2
                else:
                    factor = None
                    for key, attr in conversion_map.items():
                        if key in label:
                            factor = 1 / getattr(self.Mr, attr)
                            break
                if factor is not None:
                    dset.unit_conversion("mean", factor)
                dsets[label] = dset
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        po_sim = TimecourseSim([
            Timecourse(
                start=0,
                end=170 * 60,
                steps=2000,
                changes={
                    **self.default_changes(),
                    "PODOSE_gli": Q_(1.0, "mg"),
                },
            )
        ])

        iv_sim = TimecourseSim([
            Timecourse(
                start=0,
                end=1,
                steps=20,
                changes={
                    **self.default_changes(),
                    "Ri_gli": Q_(1.0, "mg/min"),
                },
            ),
            Timecourse(
                start=0,
                end=170.0 * 60.0 - 1,
                steps=2000,
                changes={
                    "Ri_gli": Q_(0.0, "mg/min"),
                },
            )
        ])
        return {"GLI1_po": po_sim, "GLI1_iv": iv_sim}


    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            if "urine" in name:
                tissue = Tissue.URINE
            elif "feces" in name:
                tissue = Tissue.FECES
            else:
                tissue = Tissue.SERUM

            for route in self.colors:
                if route == "iv" and tissue == Tissue.FECES:
                    continue

                dataset_label = f"{name}_GLI1{route}"

                mappings[f"fm_{route}_{name}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=dataset_label,
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd" if tissue == Tissue.FECES else None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_GLI1_{route}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=GlimepirideMappingMetaData(
                        tissue=tissue,
                        route=Route(route),
                        application_form=ApplicationForm.TABLET if route == "po" else ApplicationForm.SOLUTION,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                    ),
                )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.fig1(),
            **self.fig2(),
            **self.tab3a(),
        }

    def fig1(self) -> Dict[str, Figure]:
        # Fig1: Concentrations
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=1,
            num_cols=4,
            name=f"{self.__class__.__name__} (Plasma, Healthy)"
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True
        )
        plots[0].set_yaxis(self.label_gli_plasma, unit=self.unit_gli)
        plots[1].set_yaxis(self.label_m1_plasma, unit=self.unit_m1)
        plots[2].set_yaxis(self.label_m2_plasma, unit=self.unit_m2)
        plots[3].set_yaxis(self.label_m1_m2_plasma, unit=self.unit_m1_m2)

        for k in range(len(plots)):
            plots[k].xaxis.min = -0.5
            plots[k].xaxis.max = 25

        for route, color in self.colors.items():
            for k, substance in enumerate(["gli", "m1", "m2", "m1_m2"]):
                sid = f"[Cve_{substance}]"
                name = self.info[sid]
                plots[k].add_data(
                    task=f"task_GLI1_{route}",
                    xid="time",
                    yid=sid,
                    label=f"Sim 1mg {route.upper()}",
                    color=color,
                )
                plots[k].add_data(
                    dataset=f"{name}_GLI1{route}",
                    xid="time",
                    yid="mean",
                    count="count",
                    label=f"Data 1mg {route.upper()}",
                    color=color,
                )

        return {fig.sid: fig}

    def fig2(self) -> Dict[str, Figure]:
        # Fig2: Urine M1+M2
        fig = Figure(
            experiment=self,
            sid="Fig2",
            name=f"{self.__class__.__name__} (Urine, Healthy)"
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            yaxis=Axis(self.label_m1_m2_urine, unit=self.unit_m1_m2_urine),
            legend=True
        )

        plots[0].xaxis.max = 51
        plots[0].xaxis.min = -0.5

        for route, color in self.colors.items():
            plots[0].add_data(
                task=f"task_GLI1_{route}",
                xid="time",
                yid="Aurine_m1_m2",
                label=f"Sim 1mg {route.upper()}",
                color=color,
            )
            plots[0].add_data(
                dataset=f"urine_M1+M2_GLI1{route}",
                xid="time",
                yid="mean",
                yid_sd=None,
                count="count",
                label=f"Data 1mg {route.upper()}",
                color=color,
            )
        return {fig.sid: fig}

    def tab3a(self) -> Dict[str, Figure]:

        # Fig4: Feces M1, M2, M1+M2, and glimepiride
        fig = Figure(
            experiment=self,
            sid="Tab3",
            num_rows=1,
            name=f"{self.__class__.__name__} (Feces, Healthy)"
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis("Excretion Feces", unit=self.unit_gli_feces)
        colors = {"gli": "black", "m1": "#018571", "m2": "#dfc27d", "m1_m2": "tab:brown"}

        plots[0].xaxis.min = -0.5

        route = "po"

        for k, substance in enumerate(["gli", "m1", "m2", "m1_m2"]):
            color = colors[substance]
            sid = f"Afeces_{substance}"
            name = self.info[sid]

            plots[0].add_data(
                task=f"task_GLI1_{route}",
                xid="time",
                yid=sid,
                label=f"Sim 1mg {route.upper()} - M1+M2" if substance == "m1_m2" else f"Data 1mg {route.upper()} - {substance.upper()}",
                color=color,
            )
            plots[0].add_data(
                dataset=f"{name}_GLI1{route}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"Data 1mg {route.upper()} - M1+M2" if substance == "m1_m2" else f"Data 1mg {route.upper()} - {substance.upper()}",
                color=color,
                linestyle="",
            )

        return {fig.sid: fig}


if __name__ == "__main__":
    out = glimepiride.RESULTS_PATH_SIMULATION / FDA1995.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(FDA1995, output_dir=out)