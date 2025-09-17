"""Glimepiride intestine model."""

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

    per_hr = UnitDefinition("per_hr", "1/hr")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")


_m = Model(
    "glimepiride_intestine",
    name="Model for glimepiride absorption in the small intestine",
    notes="""
    # Model for glimepiride absorption
    - 100% fraction absorbed
    - reabsorption of glimepiride metabolites into the intestinal lumen (not biliary)
    - 40% of the dose as metabolites in feces
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
)

_m.compartments = [
    Compartment(
        "Vext",
        1.0,
        name="plasma",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["plasma"],
    ),
    Compartment(
        "Vlumen",
        value=1.0,
        name="intestinal lumen (inner part of intestine)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        constant=False,
        port=True,
        annotations=annotations.compartments["gu_lumen"],
    ),
    Compartment(
        "Vfeces",
        metaId="meta_Vfeces",
        value=1.0,
        unit=U.liter,
        constant=True,
        name="feces",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        port=True,
        annotations=annotations.compartments["feces"],
    ),
    Compartment(
        "Vstomach",
        metaId="meta_Vstomach",
        value=1,
        unit=U.liter,
        constant=True,
        name="stomach",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        port=True,
        annotations=annotations.compartments["stomach"],
    ),
]

_m.species = [
    Species(
        f"gli_stomach",
        metaId=f"meta_gli_stomach",
        initialConcentration=0.0,
        compartment="Vstomach",
        substanceUnit=U.mmole,
        name=f"glimepiride (stomach)",
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["gli"],
        boundaryCondition=True,
    ),
    Species(
        "gli_lumen",
        initialConcentration=0.0,
        name="glimepiride (intestines)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["gli"],
        port=True,
    ),
    Species(
        "gli_ext",
        initialConcentration=0.0,
        name="glimepiride (plasma)",
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["gli"],
        port=True,
    ),
    Species(
        "m1_ext",
        initialConcentration=0.0,
        name="M1 (plasma)",
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m1"],
        port=True,
    ),
    Species(
        "m2_ext",
        initialConcentration=0.0,
        name="M2 (plasma)",
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m2"],
        port=True,
    ),
    Species(
        "m1_lumen",
        initialConcentration=0.0,
        name="M1 (intestines)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m1"],
        port=True,
    ),
    Species(
        "m2_lumen",
        initialConcentration=0.0,
        name="M2 (intestines)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m2"],
        port=True,
    ),
    Species(
        "m1_feces",
        initialConcentration=0.0,
        name="M1 (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m1"],
        port=True,
    ),
    Species(
        "m2_feces",
        initialConcentration=0.0,
        name="M2 (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m2"],
        port=True,
    ),
]

_m.parameters = [
    Parameter(
        "GLIABS_k",
        0.015900835518327966,
        unit=U.per_min,
        name="rate of glimepiride absorption",
        sboTerm=SBO.KINETIC_CONSTANT,
    ),
    Parameter(
        "f_absorption",
        1,
        unit=U.dimensionless,
        name="scaling factor for absorption rate",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""1.0: normal absorption corresponding to tablet under fasting conditions;
        food decreases the absorption rate, i.e. < 1.0.
        """
    ),

    Parameter(
        "MREABS_k",
        0.015920183788739137,
        unit=U.per_min,
        name="rate constant for M1 and M2 reabsorption from plasma to intestine",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""
        Assumption: same rate for M1 and M2.
        """
    ),
    Parameter(
        "MEXC_k",
        0.0001719300784461187,
        unit=U.per_min,
        name="rate feces excretion of M1 and M2",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""
        Assumption: same rate for M1 and M2.
        """
    ),
]

_m.reactions = [
    Reaction(
        "GLIABS",
        name="absorption glimepiride",
        equation="gli_lumen -> gli_ext",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vlumen",
        formula="f_absorption * GLIABS_k * Vlumen * gli_lumen",
        notes="""Due to 100% fraction absorbed no process for excretion of gli in feces is added."""
    ),
    Reaction(
        sid="M1REABS",
        name="reabsorption M1",
        compartment="Vlumen",
        equation="m1_ext -> m1_lumen",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[],
        formula=(
            f"MREABS_k * m1_ext * Vlumen"
        ),
    ),
    Reaction(
        sid="M2REABS",
        name="reabsorption M2",
        compartment="Vlumen",
        equation="m2_ext -> m2_lumen",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[],
        formula=(
            f"MREABS_k * m2_ext * Vlumen"
        ),
    ),
    Reaction(
        sid="M1EXC",
        name=f"excretion M1",
        compartment="Vlumen",
        equation=f"m1_lumen -> m1_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[],
        formula=(
            f"MEXC_k * m1_lumen * Vlumen"
        ),
        notes="""Assumption: same rate for M1 and M2."""
    ),
    Reaction(
        sid="M2EXC",
        name=f"excretion M2",
        compartment="Vlumen",
        equation="m2_lumen -> m2_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[],
        formula=(
            f"MEXC_k * m2_lumen * Vlumen"
        ),
        notes="""Assumption: same rate for M1 and M2."""
    ),
]


_m.parameters.extend([
    Parameter(
        f"PODOSE_gli",
        0,
        unit=U.mg,
        constant=False,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"oral dose glimepiride [mg]",
        port=True,
    ),
    Parameter(
        f"Ka_dis_gli",
        2.0,
        unit=U.per_hr,
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"Ka_dis [1/hr] dissolution glimepiride",
        port=True
    ),
    Parameter(
        f"Mr_gli",
        490.616,
        unit=U.g_per_mole,
        constant=True,
        name=f"Molecular weight glimepiride [g/mole]",
        sboTerm=SBO.MOLECULAR_MASS,
        port=True,
    ),
])

# -------------------------------------------------------------------------------------------------
# Dissolution of tablet/dose in stomach
# -------------------------------------------------------------------------------------------------
_m.reactions.extend(
    [
        # fraction dose available for absorption from stomach
        Reaction(
            sid=f"dissolution_gli",
            name=f"dissolution glimepiride",
            formula=(
                f"Ka_dis_gli/60 min_per_hr * PODOSE_gli/Mr_gli",
                U.mmole_per_min
            ),
            equation=f"gli_stomach -> gli_lumen",
            compartment="Vlumen",
            sboTerm=SBO.BIOCHEMICAL_REACTION,
            notes="""Swallowing, dissolution of tablet, and transport into intestine.
            Overall process describing the rates of this processes.
            """
        ),
    ]
)
_m.rate_rules.append(
    RateRule(
        f"PODOSE_gli",
        f"-dissolution_gli * Mr_gli",
        unit=U.mg_per_min,
        name="rate of dose dissolution"
    ),
)
# -------------------------------------------------------------------------------------------------
# Feces amounts
# -------------------------------------------------------------------------------------------------
_m.parameters.append(
    Parameter(
        "mtot_feces",
        value=np.nan,
        unit=U.mmole,
        name="total fecal metabolites (M1 + M2)",
        constant=False,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    )
)
_m.rules.append(
    AssignmentRule(
        "mtot_feces",
        "m1_feces + m2_feces",
        unit=U.mmole,
        name="total fecal metabolites (M1 + M2)",
    )
)

model_intestine = _m


def intestine_layout(dx=200, dy=200) -> pd.DataFrame:
    """Intestine layout for cytoscape visualization."""

    delta_x = 0.8 * dx
    delta_y = 0.4 * dy

    positions = [

        ["gli_stomach",        0 * delta_x, 0],
        ["dissolution_gli",    0 * delta_x, 1 * delta_y],
        ["gli_lumen",          0 * delta_x, 2 * delta_y],
        ["GLIABS",             0 * delta_x, 3 * delta_y],
        ["gli_ext",            0 * delta_x, 3.75 * delta_y],

        ["m2_lumen",           1 * delta_x, 2 * delta_y],
        ["M2REABS",            1 * delta_x, 3 * delta_y],
        ["m2_ext",             1 * delta_x, 3.75 * delta_y],
        ["M2EXC",              1 * delta_x, 1 * delta_y],
        ["m2_feces",           1 * delta_x, 0 * delta_y],

        ["m1_lumen",           2 * delta_x, 2 * delta_y],
        ["M1REABS",            2 * delta_x, 3 * delta_y],
        ["m1_ext",             2 * delta_x, 3.75 * delta_y],
        ["M1EXC",              2 * delta_x, 1 * delta_y],
        ["m1_feces",           2 * delta_x, 0 * delta_y],
    ]

    df = pd.DataFrame(positions, columns=["id", "x", "y"])
    df.set_index("id", inplace=True)
    return df


def intestine_annotations(dx=200, dy=200) -> list:
    """Bounding boxes for 'stomach', 'intestines', 'plasma' and 'feces'."""
    from sbmlutils.cytoscape import AnnotationShape, AnnotationShapeType

    delta_x = 0.8 * dx
    delta_y = 0.4 * dy

    kwargs = dict(
        type=AnnotationShapeType.ROUND_RECTANGLE,
        opacity=20,
        border_color="#000000",
        border_thickness=2,
    )

    # Stomach region
    stomach_box = AnnotationShape(
        x_pos=-0.6 * delta_x,
        y_pos=-0.5 * delta_y,
        width=1.2 * delta_x,
        height=1.5 * delta_y,
        fill_color="#FFC0CB",      # pinkish
        **kwargs
    )

    # Intestinal lumen
    lumen_box = AnnotationShape(
        x_pos=-0.6 * delta_x,
        y_pos=stomach_box.y_pos + stomach_box.height,
        width=3.1 * delta_x,
        height=2 * delta_y,
        fill_color="#ffe3e3",      # very light pink
        **kwargs
    )

    # Plasma
    plasma_box = AnnotationShape(
        x_pos=-0.6 * delta_x,
        y_pos=lumen_box.y_pos + lumen_box.height,
        width=3.1 * delta_x,
        height=1.25 * delta_y,
        fill_color="#FF0000",  # pale pink
        **kwargs
    )

    # Feces
    feces_box = AnnotationShape(
        x_pos=0.6 * delta_x,
        y_pos=-0.5 * delta_y,
        width=1.9 * delta_x,
        height=1.5 * delta_y,
        fill_color="#d7c2af",  # beige
        **kwargs
    )

    return [stomach_box, lumen_box, plasma_box, feces_box]



if __name__ == "__main__":
    results: FactoryResult = create_model(
        model=model_intestine,
        filepath=MODEL_BASE_PATH / f"{model_intestine.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_factory = odefac.SBML2ODE.from_file(sbml_file=results.sbml_path)
    ode_factory.to_markdown(md_file=results.sbml_path.parent / f"{results.sbml_path.stem}.md")

    # visualization in Cytoscape
    cyviz.visualize_sbml(results.sbml_path, delete_session=False)
    cyviz.apply_layout(layout=intestine_layout())
    cyviz.add_annotations(annotations=intestine_annotations())

    # PNG output
    # cyviz.export_image(
    #    MODEL_BASE_PATH / "figures" / f"{model_intestine.sid}.png",
    #    fit_content=True,
    # )