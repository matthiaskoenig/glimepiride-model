"""
Reusable functionality for multiple simulation experiments.
"""

from collections import namedtuple
from typing import Dict
import pandas as pd
from pkdb_models.models.glimepiride import MODEL_PATH
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.model import AbstractModel
from sbmlsim.task import Task

# Constants for conversion
MolecularWeights = namedtuple("MolecularWeights", "gli m1 m2 m1_m2")


class GlimepirideSimulationExperiment(SimulationExperiment):
    """Base class for all SimulationExperiments."""

    font = {"weight": "bold", "size": 22}
    scan_font = {"weight": "bold", "size": 15}
    tick_font_size = 12
    legend_font_size = 9
    suptitle_font_size = 25

    # labels
    label_time = "Time"
    label_gli = "Glimepiride"
    label_m1 = "M1"
    label_m2 = "M2"
    label_m1_m2 = "M1+M2"

    # labels plasma
    label_gli_plasma = label_gli + " Plasma"
    label_m1_plasma = label_m1 + " Plasma"
    label_m2_plasma = label_m2 + " Plasma"
    label_m1_m2_plasma = label_m1_m2 + " Plasma"

    # labels serum
    label_gli_serum = label_gli + " Serum"
    label_m1_serum = label_m1 + " Serum"
    label_m2_serum = label_m2 + " Serum"
    label_m1_m2_serum = label_m1_m2 + " Serum"

    # labels urine
    label_m1_urine = label_m1 + " Urine"
    label_m2_urine = label_m2 + " Urine"
    label_m1_m2_urine = label_m1_m2 + " Urine"

    # labels feces
    label_gli_feces = label_gli + " Feces"
    label_m1_feces = label_m1 + " Feces"
    label_m2_feces = label_m2 + " Feces"
    label_m1_m2_feces = label_m1_m2 + " Feces"

    labels: Dict[str, str] = {
        "time": "Time",
        "[Cve_gli]": label_gli + " Plasma",
        "[Cve_m1]": label_m1 + " Plasma",
        "[Cve_m2]": label_m2 + " Plasma",
        "[Cve_m1_m2]": label_m1_m2 + " Plasma",
        "Aurine_m1": label_m1_urine,
        "Aurine_m2": label_m2_urine,
        "Aurine_m1_m2": label_m1_m2_urine,
        "Afeces_gli": label_gli_feces,
        "Afeces_m1": label_m1_feces,
        "Afeces_m2": label_m2_feces,
        "Afeces_m1_m2": label_m1_m2_feces,

        "[LI__gli]": label_gli + " Liver",
        "[LI__m1]": label_m1 + " Liver",
        "[LI__m2]": label_m2 + " Liver",

        "KI__crcl": "reatinine clearance",
    }

    # units
    unit_time = "hr"
    unit_metabolite = "µM"
    unit_metabolite_urine = "µmole"
    unit_metabolite_feces = "µmole"

    unit_gli = unit_metabolite
    unit_m1 = unit_metabolite
    unit_m2 = unit_metabolite
    unit_m1_m2 = unit_metabolite

    unit_m1_urine = unit_metabolite_urine
    unit_m2_urine = unit_metabolite_urine
    unit_m1_m2_urine = unit_metabolite_urine
    unit_gli_feces = unit_metabolite_feces
    unit_m1_feces = unit_metabolite_feces
    unit_m2_feces = unit_metabolite_feces
    unit_m1_m2_feces = unit_metabolite_feces

    units: Dict[str, str] = {
        "time": unit_time,
        "[Cve_gli]": unit_gli,
        "[Cve_m1]": unit_m1,
        "[Cve_m2]": unit_m2,
        "[Cve_m1_m2]": unit_m1_m2,
        "Aurine_m1": unit_m1_urine,
        "Aurine_m2": unit_m2_urine,
        "Aurine_m1_m2": unit_m1_m2_urine,
        "Afeces_gli": unit_gli_feces,
        "Afeces_m1": unit_m1_feces,
        "Afeces_m2": unit_m2_feces,
        "Afeces_m1_m2": unit_m1_m2_feces,

        "[LI__gli]": unit_gli,
        "[LI__m1]": unit_m1,
        "[LI__m2]": unit_m2,

        "KI__crcl": "ml/min",
    }

    # ----------- Renal map --------------
    renal_map = {
        "Normal renal function": 101.0 / 101.0,  # 1.0,
        "Mild renal impairment": 50.0 / 101.0,  # 0.5
        "Moderate renal impairment": 35.0 / 101.0,  # 0.35
        "Severe renal impairment": 20.0 / 101.0,  # 0.20
        # "End stage renal disease": 10.5 / 101.0,  # 0.1
    }
    renal_colors = {
        "Normal renal function": "black",
        "Mild renal impairment": "#66c2a4",
        "Moderate renal impairment": "#2ca25f",
        "Severe renal impairment": "#006d2c",
        # "End stage renal disease": "#006d5e"
    }

    # ----------- Cirrhosis map --------------
    cirrhosis_map = {
        "Control": 0,
        "Mild cirrhosis": 0.3994897959183674,  # CPT A
        "Moderate cirrhosis": 0.6979591836734694,  # CPT B
        "Severe cirrhosis": 0.8127551020408164,  # CPT C
    }
    cirrhosis_colors = {
        "Control": "black",
        "Mild cirrhosis": "#74a9cf",  # CPT A
        "Moderate cirrhosis": "#2b8cbe",  # CPT B
        "Severe cirrhosis": "#045a8d",  # CPT C
    }

    # ----------- Genotypes --------------
    # see https://www.pharmgkb.org/page/cyp2c9RefMaterials

    cyp2c9_allele_activity = {
        "*1": 1.00,  # [Reference allele]
        "*2": 0.68,  # [Dai2014, Yang2018]
        "*3": 0.23,  # [Dai2014, Maekawa2009, Suzuki2006, Yang2018]
        "*8": 0.09,  # [Dai2014]
        "*11": 0.61,  # [Dai2014]
        "*13": 0.02,  # [Dai2014]
        "*14": 0.06,  # [Dai2014]
        "*16": 0.04,  # [Dai2014]
        "*19": 0.01,  # [Dai2014]
        "*23": 0.07,  # [Dai2014]
        "*26": 0.10,  # [Dai2014]
        "*27": 0.15,  # [Dai2014]
        "*28": 0.44,  # [Dai2014]
        "*29": 0.36,  # [Dai2014]
        "*30": 0.25,  # [Dai2014]
        "*31": 0.14,  # [Dai2014]
        "*33": 0.04,  # [Dai2014]
        "*34": 0.41,  # [Dai2014]
        "*36": 1.46,  # [Dai2014]
        "*37": 0.75,  # [Dai2014]
        "*38": 0.64,  # [Dai2014]
        "*39": 0.10,  # [Dai2014]
        "*40": 1.03,  # [Dai2014]
        "*41": 0.75,  # [Dai2014]
        "*42": 0.03,  # [Dai2014]
        "*43": 0.07,  # [Dai2014]
        "*44": 0.15,  # [Dai2014]
        "*45": 0.08,  # [Dai2014]
        "*46": 0.23,  # [Dai2014]
        "*47": 1.15,  # [Dai2014]
        "*48": 0.67,  # [Dai2014]
        "*49": 0.51,  # [Dai2014]
        "*50": 0.29,  # [Dai2014]
        "*51": 0.91,  # [Dai2014]
        "*52": 0.03,  # [Dai2014]
        "*53": 0.87,  # [Dai2014]
        "*54": 0.95,  # [Dai2014]
        "*55": 0.11,  # [Dai2014]
        "*56": 0.81,  # [Dai2014]

        "no": 0.0,  # [No activity]
    }

    cyp2c9_activity = {
        "*1/*1": (cyp2c9_allele_activity["*1"] + cyp2c9_allele_activity["*1"])/2,
        "*1/*2": (cyp2c9_allele_activity["*1"] + cyp2c9_allele_activity["*2"])/2,
        "*1/*3": (cyp2c9_allele_activity["*1"] + cyp2c9_allele_activity["*3"])/2,
        "*3/*3": (cyp2c9_allele_activity["*3"] + cyp2c9_allele_activity["*3"])/2,

        #"no": (cyp2c9_allele_activity["no"] + cyp2c9_allele_activity["no"]) / 2,
    }

    # cyp2c9 colors for time course plots
    cyp2c9_colors = {
        "*1/*1": "black",
        "*1/*2": "#b589d6",
        "*1/*3": "#804fb3",
        "*3/*3": "#552586",
        #"no": "grey",
    }

    # cyp2c9 colors for all other analyses
    cyp2c9_colors_analyses = {
        "*1/*1": "lightgrey",
        "*1/*2": "#b589d6",
        "*1/*3": "#804fb3",
        "*3/*3": "#552586",
    }

    # dose_values = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    # dose_palette = ["#FFF0E0", "#FFDAB9", "#FFBC7D", "#FFA54F", "#FF8C00", "#F06000", "#D04A00", "#A83800", "#8B2500"]
    # dose_colors: Dict[float, str] = dict(zip(dose_values, dose_palette))

    dose_values = [0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4, 5, 6, 7, 8]
    dose_palette = ["#FFFCF0", "#FFEDCC", "#FFDB99", "#FFC966", "#FFB533", "#FFA500", "#FF9933", "#FF8C00", "#FF7F00",
                    "#FF6600", "#E55100", "#CC4400", "#B33300"]
    dose_colors: Dict[float, str] = dict(zip(dose_values, dose_palette))

    bodyweights = [45, 70, 95, 120, 145, 170]
    bodyweight_palette = ["#f5c6d0", "#de93b0", "#c76091", "#a93b72", "#842655", "#601038"]
    bodyweight_colors: Dict[float, str] = dict(zip(bodyweights, bodyweight_palette))

    def models(self) -> Dict[str, AbstractModel]:
        Q_ = self.Q_
        return {
            "model": AbstractModel(
                source=MODEL_PATH,
                language_type=AbstractModel.LanguageType.SBML,
                changes={},
            )
        }

    @staticmethod
    def _default_changes(Q_):
        """Default changes to simulations."""

        changes = {
            # 20250410_143850__6baa2/GLIMEPIRIDE_LSQ_CYP2C9
            # 'LI__GLIIM_k': Q_(100, '1/min'),  # [0.001 - 100]
            # 	>>> !Optimal parameter 'Kp_gli' within 5% of lower bound! <<<
            # 	>>> !Optimal parameter 'LI__M2EX_k' within 5% of upper bound! <<<
            # 	'ftissue_gli': Q_(0.00070864179236918, 'l/min'),  # [0.0001 - 1]
            # 	'Kp_gli': Q_(10.020596565282165, 'dimensionless'),  # [10 - 1000]
            # 	'GU__GLIABS_k': Q_(0.015900835518327966, '1/min'),  # [0.0001 - 1]
            # 	'GU__MREABS_k': Q_(0.015920183788739137, '1/min'),  # [0.0001 - 1]
            # 	'GU__MEXC_k': Q_(0.0001719300784461187, '1/min'),  # [0.0001 - 1]
            # 	'LI__GLI2M1_Vmax': Q_(5.3620361192856574e-05, 'mmole/min/l'),  # [1e-07 - 0.01]
            # 	'LI__M1EX_k': Q_(0.07774345977190338, '1/min'),  # [0.001 - 100]
            # 	'LI__M12M2_k': Q_(0.014851987496532965, '1/min'),  # [0.001 - 100]
            # 	'LI__M2EX_k': Q_(99.99817447866164, '1/min'),  # [0.001 - 100]
            # 	'KI__M1EX_k': Q_(0.14800981624897125, '1/min'),  # [0.0001 - 1]
            # 	'KI__M2EX_k': Q_(0.09849263263810941, '1/min'),  # [0.0001 - 1]

        }

        return changes

    def default_changes(self: SimulationExperiment) -> Dict:
        """Default changes to simulations."""
        return GlimepirideSimulationExperiment._default_changes(Q_=self.Q_)

    def tasks(self) -> Dict[str, Task]:
        if self.simulations():
            return {
                f"task_{key}": Task(model="model", simulation=key)
                for key in self.simulations()
            }
        return {}

    def data(self) -> Dict:
        self.add_selections_data(
            selections=[
                "time",

                # dosing
                "IVDOSE_gli",
                "PODOSE_gli",
                "IVDOSE_m1",

                # venous plasma
                "[Cve_gli]",
                "[Cve_m1]",
                "[Cve_m2]",
                "[Cve_m1_m2]",

                # urine
                "Aurine_m1",
                "Aurine_m2",
                "Aurine_m1_m2",

                # feces
                "Afeces_gli",
                "Afeces_m1",
                "Afeces_m2",
                "Afeces_m1_m2",

                # cases
                'KI__f_renal_function',
                'f_cirrhosis',
                'LI__f_cyp2c9',

                # liver
                "[LI__gli]",
                "[LI__m1]",
                "[LI__m2]",

                # kidney
                "KI__crcl",  # creatinine clearance
            ]
        )
        return {}

    @property
    def Mr(self):
        return MolecularWeights(
            gli=self.Q_(490.616, "g/mole"),
            m1=self.Q_(506.62, "g/mole"),
            m2=self.Q_(520.6, "g/mole"),
            m1_m2=self.Q_(513.6, "g/mole") # mean metabolites weight
        )