"""Kidney model for glimepiride."""

import pandas as pd
from sbmlutils.factory import create_model, FactoryResult
from sbmlutils.converters import odefac
from sbmlutils import cytoscape as cyviz
import numpy as np
from sbmlutils.converters import odefac
from sbmlutils.cytoscape import visualize_sbml
from sbmlutils.factory import *
from sbmlutils.metadata import *

from pkdb_models.models.glimepiride.models import annotations
from pkdb_models.models.glimepiride.models import templates
from pkdb_models.models.glimepiride import MODEL_BASE_PATH


class U(templates.U):
    """UnitDefinitions"""

    mg_per_g = UnitDefinition("mg_per_g", "mg/g")
    ml_per_l = UnitDefinition("ml_per_l", "ml/l")
    ml_per_min = UnitDefinition("ml_per_min", "ml/min")
    ml_per_min_bsa = UnitDefinition("ml_per_min_bsa", "ml/min/(1.73 * m^2)")

mid = "glimepiride_kidney"
version = 1

_m = Model(
    sid=mid,
    name="Model for renal M1 and M2 excretion.",
    notes=f"""
    Model for renal M1 and M2 excretion.
    
    - Glimepiride is fully metabolized in the liver, and no unchanged drug is excreted [FDA label].
    - ~60 % of the total radioactivity is excreted in urine (M1: predominant), and ~40% is 
      excreted in feces (M2: predominant) [FDA label].
    
    **version** {version}
    
    ## Changelog
    **version 1**
    - initial model for glimepiride
    """ + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
)

_m.compartments = [
    Compartment(
        "Vext",
        value=1.5,
        unit=U.liter,
        name="plasma",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["plasma"],
        port=True
    ),
    Compartment(
        "Vki",
        value=0.3,
        unit=U.liter,
        name="kidney",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ki"],
        port=True
    ),
    Compartment(
        "Vurine",
        1.0,
        name="urine",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["urine"],
    ),
]

_m.species = [
    Species(
        "m1_ext",
        name="M1 (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m1"],
        port=True
    ),
    Species(
        "m1_urine",
        name="M1 (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m1"],
        port=True
    ),
    Species(
        "m2_ext",
        name="M2 (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m2"],
        port=True
    ),
    Species(
        "m2_urine",
        name="M2 (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m2"],
        port=True
    ),
]


# glomerular filtration rate, creatinine clearance & kidney function
_m.parameters.extend([
    Parameter(
        "BSA",
        1.73,
        U.m2,
        constant=False,
        name="body surface area [m^2]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        port=True
    ),
    Parameter(
        "f_renal_function",
        name="parameter for renal function",
        value=1.0,
        unit=U.dimensionless,
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""scaling factor for renal function. 1.0: normal renal function; 
        <1.0: reduced renal function
        """
    ),
    Parameter(
        "egfr",
        np.nan,
        unit=U.ml_per_min_bsa,
        constant=False,
        name="estimated GFR [ml/min/(1.73 m^2)]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations = [
            (BQB.IS, "ncit/C110935"),
        ],
        port=True,
        notes="""A laboratory test that estimates kidney function. It is calculated using an individual's 
        serum creatinine measurement, age, gender, and race."""
    ),
    Parameter(
        "crcl",
        np.nan,
        U.ml_per_min,
        constant=False,
        name="creatinine clearance [ml/min]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS, "cmo/CMO:0000765"),
        ],
        port=True,
        notes="""The clearance rate of creatinine, that is, the volume of plasma that is cleared of creatinine 
        by the kidneys per unit time. Creatinine clearance is calculated using the level of creatinine in a 
        sample of urine, usually one collected over a period of 24 hours, the corresponding plasma creatinine 
        level, and the volume of urine excreted. It is used as an approximation of the glomerular 
        filtration rate (GFR)."""
    ),
    Parameter(
        sid="egfr_healthy",
        name="estimated glomerular filtration (eGFR) rate healthy",
        value=100.0,
        unit=U.ml_per_min_bsa,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
        ],
        notes="""eGFR is often estimated via creatinine clearance, creatinine urinary amount,
        or creatinine plasma amounts.
        
        CKD-EPI Creatinine Equation (2021): This is the most recent and recommended equation for 
        estimating GFR35. It is more accurate than previous formulas and does not include 
        race as a variable.
        
        MDRD (Modification of Diet in Renal Disease) Study Equation
    """,
    )
])
_m.rules.extend([
    AssignmentRule(
        "egfr", "f_renal_function * egfr_healthy", unit=U.ml_per_min_bsa,
        name="estimated eGFR"
    ),
    AssignmentRule(
        "crcl", "egfr * BSA/(1.73 dimensionless) * 1.1 dimensionless", unit=U.ml_per_min,
        name="creatinine clearance",
        notes="CrCl typically overestimates GFR by 10-20% due to the active secretion"
              "of creatinine in the proximal tubules."
              "Using a 10 % overestimation in the model."
    ),
])

_m.reactions = [
    Reaction(
        sid="M1EX",
        name="M1 excretion (M1EX)",
        equation="m1_ext -> m1_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "M1EX_k",
                0.14800981624897125,
                U.per_min,
                name="rate of M1 urinary excretion",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=(
            "f_renal_function * Vki * M1EX_k * m1_ext"
        )
    ),
    Reaction(
        sid="M2EX",
        name="M2 excretion (M2EX)",
        equation="m2_ext -> m2_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "M2EX_k",
                0.09849263263810941,
                U.per_min,
                name="rate of M2 urinary excretion",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=(
            "f_renal_function * Vki * M2EX_k * m2_ext"
        )
    ),
]

model_kidney = _m


def kidney_layout(dx=200, dy=200) -> pd.DataFrame:
    """Kidney layout for Cytoscape visualization."""
    delta_x = 0.8 * dx
    delta_y = 0.4 * dy

    positions = [
        ["m1_ext", 0, 0],
        ["M1EX",   0, 1 * delta_y],
        ["m1_urine", 0, 1.75 * delta_y],

        ["m2_ext", 1 * delta_x, 0],
        ["M2EX",   1 * delta_x, 1 * delta_y],
        ["m2_urine", 1 * delta_x, 1.75 * delta_y],
    ]

    df = pd.DataFrame(positions, columns=["id", "x", "y"])
    df.set_index("id", inplace=True)
    return df

def kidney_annotations(dx=200, dy=200) -> list:
    """Bounding boxes for 'plasma' and 'urine'."""
    from sbmlutils import cytoscape as cyviz

    delta_x = 0.8 * dx
    delta_y = 0.4 * dy

    kwargs = {
        "type": cyviz.AnnotationShapeType.ROUND_RECTANGLE,
        "opacity": 20,
        "border_color": "#000000",
        "border_thickness": 2,
    }

    # plasma
    top_box = cyviz.AnnotationShape(
        x_pos=-0.5 * delta_x,
        y_pos=-0.5 * delta_y,
        width=2.0 * delta_x,
        height=1.5 * delta_y,
        fill_color="#FF0000",  # pale pink
        **kwargs
    )

    # urine
    bottom_box = cyviz.AnnotationShape(
        x_pos=-0.5 * delta_x,
        y_pos=top_box.y_pos + top_box.height,
        width=2.0 * delta_x,
        height=1.25 * delta_y,
        fill_color="#FFFACD",  # pale yellow
        **kwargs
    )

    return [top_box, bottom_box]



if __name__ == "__main__":
    results: FactoryResult = create_model(
        model=model_kidney,
        filepath=MODEL_BASE_PATH / f"{model_kidney.sid}.xml",
        sbml_level=3, sbml_version=2,
    )

    # ODE equations
    ode_factory = odefac.SBML2ODE.from_file(sbml_file=results.sbml_path)
    ode_factory.to_markdown(md_file=results.sbml_path.parent / f"{results.sbml_path.stem}.md")

    # visualization in Cytoscape
    cyviz.visualize_sbml(sbml_path=results.sbml_path, delete_session=False)
    cyviz.apply_layout(layout=kidney_layout())
    cyviz.add_annotations(annotations=kidney_annotations())

    # PNG output
    # cyviz.export_image(
    #    MODEL_BASE_PATH / "figures" / f"{model_kidney.sid}.png",
    #    fit_content=True,
    # )

