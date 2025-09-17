"""Liver model for glimepiride."""

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

    pass


mid = "glimepiride_liver"
version = 1

_m = Model(
    sid=mid,
    name="Model for hepatic glimepiride metabolism.",
    notes=f"""
    Model for hepatic glimepiride metabolism.
    Glimepiride is completely biotransformed by hepatic oxidative metabolism. The CYP2C9 enzyme transforms glimepiride 
    to the cyclohexylhydroxymethyl derivative (M1), which is further metabolized to form the carboxy derivative (M2) 
    by cytosolic enzymes. [FDA label]
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
        "Vli",
        value=1.5,
        unit=U.liter,
        name="liver",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["li"],
        port=True
    ),
]

_m.species = [
    Species(
        "gli_ext",
        name="glimepiride (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["gli"],
        port=True
    ),
    Species(
        "gli",
        name="glimepiride (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["gli"],
        port=True
    ),
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
        "m1",
        name="M1 (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
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
        "m2",
        name="M2 (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m2"],
        port=True
    ),
]


_m.reactions = [
    Reaction(
        sid="GLIIM",
        name="glimepiride import",
        equation="gli_ext <-> gli",
        compartment="Vext",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "GLIIM_k",
                100,
                U.per_min,
                name="rate glimepiride import",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
        ],
        formula=(
            "GLIIM_k * Vli * (gli_ext - gli)"
        ),
    ),
    Reaction(
        sid="GLI2M1",
        name="glimepiride → M1 (CYP2C9)",
        equation="gli -> m1",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "GLI2M1_Vmax",
                5.3620361192856574e-05,
                U.mmole_per_min_l,
                name="rate glimepiride conversion",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
            Parameter(
                "GLI2M1_Km_gli",
                0.00225,
                U.mM,
                name="Km gli",
                sboTerm=SBO.MICHAELIS_CONSTANT,
                notes="""
                0.18 +- 0.03 µM [Maekawa2009]
                2.25 +- 0.80 µM [Zhang2023]
                4.1 +- 0.4 µM [Suzuki2006]
                """
            ),
            Parameter(
                "f_cyp2c9",
                1,
                U.dimensionless,
                name="scaling factor CYP2C9 activity",
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                notes="""Scaling factor to vary CYP2C9 activity.
                1.0: unchanged activity; < 1.0 decreased activity; >1.0 increased activity.
                """
            )
        ],
        formula=(
            "f_cyp2c9 * GLI2M1_Vmax * Vli * gli/(gli + GLI2M1_Km_gli)"
        )
    ),
    Reaction(
        sid="M1EX",
        name="M1 export",
        equation="m1 <-> m1_ext",
        compartment="Vli",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "M1EX_k",
                0.07774345977190338,
                unit=U.per_min,
                name="rate M1 export",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
        ],
        formula=(
            "M1EX_k * Vli * (m1 - m1_ext)"
        )
    ),

    Reaction(
        sid="M12M2",
        name="M1 → M2",
        equation="m1 -> m2",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "M12M2_k",
                0.014851987496532965,
                U.per_min,
                name="rate m1 -> m2 conversion",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
        ],
        formula=(
            "M12M2_k * Vli * m1"
        )
    ),
    Reaction(
        sid="M2EX",
        name="M2 export",
        equation="m2 <-> m2_ext",
        compartment="Vli",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "M2EX_k",
                99.99817447866164,
                unit=U.per_min,
                name="rate M2 export",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            )
        ],
        formula=(
            "M2EX_k * Vli * (m2 - m2_ext)"
        )
    ),
]

model_liver = _m


def liver_layout(dx=200, dy=200) -> pd.DataFrame:
    """Liver layout for Cytoscape visualization."""

    delta_x = 0.8 * dx
    delta_y = 0.4 * dy

    positions = [
        ["gli_ext", 0, 0],
        ["GLIIM", 0, 1 * delta_y],
        ["gli", 0, 2 * delta_y],

        ["m1_ext", 1 * delta_x, 0],
        ["M1EX", 1 * delta_x, 1 * delta_y],
        ["m1", 1 * delta_x, 2 * delta_y],

        ["m2_ext", 2 * delta_x, 0],
        ["M2EX", 2 * delta_x, 1 * delta_y],
        ["m2", 2 * delta_x, 2 * delta_y],

        ["GLI2M1", 0.5 * delta_x, 3 * delta_y],
        ["M12M2", 1.5 * delta_x, 3 * delta_y],
    ]

    df = pd.DataFrame(positions, columns=["id", "x", "y"])
    df.set_index("id", inplace=True)
    return df

def liver_annotations(dx=200, dy=200) -> list:
    """Bounding boxes for 'plasma' and 'liver'."""

    delta_x = 0.8 * dx
    delta_y = 0.4 * dy

    kwargs = {
        "type": cyviz.AnnotationShapeType.ROUND_RECTANGLE,
        "opacity": 20,
        "border_color": "#000000",
        "border_thickness": 2,
    }
    annotations = [

        #plasma
        cyviz.AnnotationShape(
            x_pos=-0.5 *delta_x,
            y_pos=-0.5 * delta_y,
            width=2.8 * delta_x,
            height=1.5 * delta_y,
            fill_color="#FF0000",  # pale pink
            **kwargs
        ),

        # liver
        cyviz.AnnotationShape(
            x_pos=-0.5 *delta_x,
            y_pos=delta_y,
            width=2.8 * delta_x,
            height=2.25 * delta_y,
            fill_color="#FFFFFF",
            **kwargs
        ),
    ]
    return annotations


if __name__ == "__main__":
    results: FactoryResult = create_model(
        model=model_liver,
        filepath=MODEL_BASE_PATH / f"{model_liver.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_factory = odefac.SBML2ODE.from_file(sbml_file=results.sbml_path)
    ode_factory.to_markdown(md_file=results.sbml_path.parent / f"{results.sbml_path.stem}.md")

    # visualization in Cytoscape
    cyviz.visualize_sbml(results.sbml_path, delete_session=True)
    cyviz.apply_layout(layout=liver_layout())
    cyviz.add_annotations(annotations=liver_annotations())

    # PNG output
    # cyviz.export_image(
    #     MODEL_BASE_PATH / "figures" / f"{model_liver.sid}.png",
    #     fit_content=True,
    # )