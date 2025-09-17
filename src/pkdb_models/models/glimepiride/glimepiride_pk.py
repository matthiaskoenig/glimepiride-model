import pandas as pd
from pkdb_analysis.pk.pharmacokinetics import TimecoursePK
from sbmlsim.result import XResult


def process_substance_pk(experiment, xres, scandim, dose_index, dose_value, substance, keys):
    """Process PK calculations for a substance."""
    Q_ = experiment.Q_

    # Get time and concentration vectors
    t_vec = Q_(xres.dim_mean("time").magnitude, xres.uinfo["time"])
    conc_key = keys["conc_key"]
    c_vec = Q_(
        xres[conc_key].sel({scandim: dose_index}).values,
        xres.uinfo[conc_key]
    )

    # Set up TimecoursePK parameters
    tpkw = {
        "time": t_vec,
        "concentration": c_vec,
        "substance": substance,
        "ureg": experiment.ureg,
        "dose": dose_value if keys.get("dose_used", False) else None
    }

    if keys.get("dose_used", False) and "min_threshold" in keys:
        tpkw["min_treshold"] = keys["min_threshold"]

    # Calculate basic PK parameters
    tcpk = TimecoursePK(**tpkw)
    pk_dict = tcpk.pk.to_dict()
    pk_dict["substance"] = substance

    # For metabolites: calculate additional clearance parameters
    if substance in ["M1", "M2"]:
        # Renal clearance
        aurine_key = keys["aurine_key"]
        aurine_vec = Q_(
            xres[aurine_key].sel({scandim: dose_index}).values,
            xres.uinfo[aurine_key]
        )
        pk_dict[aurine_key] = aurine_vec.magnitude[-1]
        pk_dict[f"{aurine_key}_unit"] = aurine_vec.units
        cl_renal = aurine_vec[-1] / tcpk.pk.auc
        pk_dict["cl_renal"] = cl_renal.magnitude
        pk_dict["cl_renal_unit"] = cl_renal.units

        # Fecal clearance
        afeces_key = keys["afeces_key"]
        afeces_vec = Q_(
            xres[afeces_key].sel({scandim: dose_index}).values,
            xres.uinfo[afeces_key]
        )
        pk_dict[afeces_key] = afeces_vec.magnitude[-1]
        pk_dict[f"{afeces_key}_unit"] = afeces_vec.units
        cl_fecal = afeces_vec[-1] / tcpk.pk.auc
        pk_dict["cl_fecal"] = cl_fecal.magnitude
        pk_dict["cl_fecal_unit"] = cl_fecal.units

        # Total clearance (sum of renal and fecal)
        pk_dict["cl_total"] = pk_dict["cl_renal"] + pk_dict["cl_fecal"]
        pk_dict["cl_total_unit"] = pk_dict["cl_renal_unit"]

        # Apparent clearance (dose/AUC)
        if dose_value is not None:
            dose_mmol = dose_value / experiment.Mr.gli
            # Correction factor for metabolites
            if substance == "M1":
                correction_factor = experiment.Mr.gli / experiment.Mr.m1
            elif substance == "M2":
                correction_factor = experiment.Mr.gli / experiment.Mr.m2
            else:
                correction_factor = 1.0
            cl_apparent = dose_mmol / tcpk.pk.auc * correction_factor
            pk_dict["cl"] = cl_apparent.magnitude
            pk_dict["cl_unit"] = cl_apparent.units
        else:
            pk_dict["cl"] = float('nan')
            pk_dict["cl_unit"] = tcpk.pk.auc.units.replace('*min', '')

        # Add AUC and volume of distribution
        pk_dict["auc"] = tcpk.pk.auc.magnitude
        pk_dict["auc_unit"] = tcpk.pk.auc.units
        pk_dict["vd"] = pk_dict['cl'] / tcpk.pk.kel.magnitude
        pk_dict["vd_unit"] = f"{pk_dict['cl_unit']}/{tcpk.pk.kel.units}"

    return pk_dict


def calculate_glimepiride_pk(experiment, xres):
    """Calculate PK parameters for glimepiride, M1, and M2."""
    # Get scanned dimension and dose vector
    scandim = xres._redop_dims()[0]
    dose_vec = experiment.Q_(xres["PODOSE_gli"].values[0], xres.uinfo["PODOSE_gli"])

    # Define substance info dictionaries
    substance_info = {
        "glimepiride": {
            "conc_key": "[Cve_gli]",
            "dose_used": True,
            "min_threshold": 100,
        },
        "M1": {
            "conc_key": "[Cve_m1]",
            "aurine_key": "Aurine_m1",
            "afeces_key": "Afeces_m1",
            "dose_used": False,
        },
        "M2": {
            "conc_key": "[Cve_m2]",
            "aurine_key": "Aurine_m2",
            "afeces_key": "Afeces_m2",
            "dose_used": False,
        },
    }

    # Process each substance for each dose
    pk_dicts = []
    for substance, keys in substance_info.items():
        for idx, dose in enumerate(dose_vec):
            pk_dict = process_substance_pk(experiment, xres, scandim, idx, dose, substance, keys)
            pk_dicts.append(pk_dict)

    return pd.DataFrame(pk_dicts)