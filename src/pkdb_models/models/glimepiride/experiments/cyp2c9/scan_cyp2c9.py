"""CYP2C9 scan for glimepiride pharmacokinetics.
Computationally intensive simulation that only produces TSV files, used by other CYP2C9 analysis scripts."""

from typing import Dict
import numpy as np
import pandas as pd
import xarray as xr
from sbmlsim.model import RoadrunnerSBMLModel
from sbmlsim.result import XResult
from sbmlsim.simulation import Timecourse, TimecourseSim, ScanSim, Dimension
from pkdb_models.models.glimepiride.experiments.cyp2c9.analyze_cyp2c9_activity import get_lognorm_pars
from pkdb_models.models.glimepiride.experiments.base_experiment import GlimepirideSimulationExperiment
from pkdb_models.models.glimepiride.helpers import run_experiments
from pkdb_models.models.glimepiride.experiments.cyp2c9.population_sampling import sample_diplotype
from pkdb_models.models.glimepiride import DATA_PATH_GLIMEPIRIDE, RESULTS_PATH
from pkdb_models.models.glimepiride.glimepiride_pk import calculate_glimepiride_pk

# Known CYP2C9 allele activities
cyp2c9_allele_activity = GlimepirideSimulationExperiment.cyp2c9_allele_activity

def _filter_diplotype(diplotype_str, allele_dict):
    """Check if both alleles exist in allele-activity dict."""
    a1, a2 = diplotype_str.split("/")
    return (a1 in allele_dict) and (a2 in allele_dict)

def filename(s: str) -> str:
    """Replace characters that may cause invalid filenames."""
    return s.replace(" ", "_").replace("/", "_").replace("\\", "_")


class GlimepirideCYP2C9Scan(GlimepirideSimulationExperiment):
    """Genotype-based PK simulations for glimepiride."""

    lognorm_pars = get_lognorm_pars()
    shape_all = lognorm_pars["All"]["shape"]
    genotypes_of_interest = ["*1/*1", "*1/*2", "*1/*3", "*3/*3"]
    n_samples_all = 1000000

    genotypes_all_list = []
    activities_all_list = []

    for genotype in genotypes_of_interest:
        allele1, allele2 = genotype.split("/")
        genotype_samples = sample_diplotype(allele1, allele2, shape_all, n=n_samples_all)
        genotypes_all_list.append([genotype]*n_samples_all)
        activities_all_list.append(genotype_samples)

    genotypes_all = np.concatenate(genotypes_all_list)
    activities_all = np.concatenate(activities_all_list)

    # Ethnicity approach
    excel_path = DATA_PATH_GLIMEPIRIDE / "cyp2c9_data.xlsx"
    df_diplotype_original = pd.read_excel(excel_path, sheet_name="diplotype_frequencies", comment="#")

    population_cols = [
        "African American/Afro-Caribbean", "American", "Central/South Asian",
        "East Asian", "European", "Latino", "Near Eastern", "Oceanian", "Sub-Saharan African"
    ]

    # Filter diplotypes with non-zero frequencies and valid alleles
    df_filled = df_diplotype_original[population_cols].fillna(0.0)
    df_diplotype = df_diplotype_original[df_filled.sum(axis=1) > 0].copy()
    df_diplotype["mean_freq"] = df_diplotype[population_cols].fillna(0).mean(axis=1)
    mask_valid = df_diplotype["CYP2C9_allele"].apply(lambda x: _filter_diplotype(x, cyp2c9_allele_activity))
    df_diplotype = df_diplotype[mask_valid].sort_values(by="mean_freq", ascending=False)

    n_samples_ethnicity = 1000000
    samples_by_ethnicity = {}

    for eth in population_cols:
        population_data = []
        for _, row in df_diplotype.iterrows():
            genotype = row["CYP2C9_allele"]
            freq = row.get(eth, 0.0)
            if pd.isna(freq):
                freq = 0.0

            count = int(n_samples_ethnicity * freq)
            if count < 1:
                continue

            a1, a2 = genotype.split("/")
            genotype_samples = sample_diplotype(a1, a2, shape_all, n=count)
            population_data.extend(genotype_samples)

        samples_by_ethnicity[eth] = np.array(population_data)

    task_eth_map = {}

    def simulations(self) -> Dict[str, ScanSim]:
        """Create ScanSims for all approaches."""
        Q_ = self.Q_
        sim_dict = {}

        sim_dict["all"] = ScanSim(
            simulation=TimecourseSim(
                Timecourse(
                    start=0,
                    end=24 * 60,
                    steps=2500,
                    changes={
                        **self.default_changes(),
                        "PODOSE_gli": Q_(4, "mg"),
                        "KI__f_renal_function": Q_(1.0, "dimensionless"),
                        "f_cirrhosis": Q_(0.0, "dimensionless"),
                    },
                )
            ),
            dimensions=[
                Dimension(
                    "dim_scan",
                    changes={"LI__f_cyp2c9": Q_(self.activities_all, "dimensionless")},
                )
            ],
        )

        # Ethnicity-specific approaches
        for eth, arr in self.samples_by_ethnicity.items():
            safe_eth = filename(eth)
            self.task_eth_map[safe_eth] = eth
            sim_dict[safe_eth] = ScanSim(
                simulation=TimecourseSim(
                    Timecourse(
                        start=0,
                        end=24 * 60,
                        steps=2500,
                        changes={
                        **self.default_changes(),
                        "PODOSE_gli": Q_(4, "mg"),
                        "KI__f_renal_function": Q_(1.0, "dimensionless"),
                        "f_cirrhosis": Q_(0.0, "dimensionless"),
                        },
                    )
                ),
                dimensions=[
                    Dimension(
                        "dim_scan",
                        changes={"LI__f_cyp2c9": Q_(arr, "dimensionless")},
                    )
                ],
            )
        return sim_dict

    def data(self) -> Dict:
        """Specify model variables to record."""
        self.add_selections_data(
            selections=[
                "time",
                "PODOSE_gli",
                # venous plasma
                "[Cve_gli]",
                "[Cve_m1]",
                "[Cve_m2]",
                # urine
                "Aurine_m1",
                "Aurine_m2",
                # feces
                "Afeces_gli",
                "Afeces_m1",
                "Afeces_m2",
                'LI__f_cyp2c9',
            ]
        )
        return {}

    def calculate_pharmacokinetics_parameters(self, output_dir: str):
        """Compute PK parameters and save TSV files."""
        task_keys = ["all"] + list(self.task_eth_map.keys())
        df_ethnicities = []

        for task_key in task_keys:
            nc_file = RESULTS_PATH / output_dir / output_dir / f"GlimepirideCYP2C9Scan_task_{task_key}.nc"
            dset = xr.open_dataset(nc_file)

            abstract_model = self.models()["model"]
            model = RoadrunnerSBMLModel.from_abstract_model(
                abstract_model=abstract_model, ureg=self.ureg
            )
            xres = XResult(xdataset=dset, uinfo=model.uinfo)
            df_pk = calculate_glimepiride_pk(experiment=self, xres=xres)

            if task_key == "all":
                df_pk.insert(0, "genotype", np.concatenate([self.genotypes_all] * 3))
                df_pk.insert(1, "f_cyp2c9", np.concatenate([self.activities_all] * 3))
                df_pk.insert(0, "task", "all_genotypes")

                # Save 'all' approach results
                df_pk.to_csv(
                    RESULTS_PATH / output_dir / output_dir / "GlimepirideCYP2C9Scan_cyp2c9_pharmacokinetics.tsv",
                    sep="\t", index=False
                )
            else:
                df_pk.insert(0, "ethnicity", self.task_eth_map[task_key])
                df_ethnicities.append(df_pk)

        # Save ethnicity results
        if df_ethnicities:
            pd.concat(df_ethnicities, ignore_index=True).to_csv(
                RESULTS_PATH / output_dir / output_dir / "GlimepirideCYP2C9Scan_ethnicities_pharmacokinetics.tsv",
                sep="\t", index=False
            )


if __name__ == "__main__":
    output_dir = "GlimepirideCYP2C9Scan"
    run_experiments(GlimepirideCYP2C9Scan, output_dir=output_dir, save_results=True)
    GlimepirideCYP2C9Scan().calculate_pharmacokinetics_parameters(output_dir)