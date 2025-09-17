"""Glimepiride model factory."""
from pathlib import Path
from typing import Dict, Any

from sbmlutils.converters import odefac


from pkdb_models.models.glimepiride.models.model_kidney import model_kidney
from pkdb_models.models.glimepiride.models.model_liver import model_liver
from pkdb_models.models.glimepiride.models.model_intestine import model_intestine
from pkdb_models.models.glimepiride.models.model_body import model_body


from sbmlutils.comp import flatten_sbml
from sbmlutils.console import console
from sbmlutils.cytoscape import visualize_sbml
from sbmlutils.factory import create_model
from pymetadata.omex import *


def create_models(
    model_output_dir: Path, create_tissues: bool = True
) -> Dict[str, Path]:
    """Creates tissue and whole-body model."""
    results: Dict[str, Dict[str, Any]] = {
        "README": {
            "path": MODEL_BASE_PATH.parent / "README.md",
            "entry": ManifestEntry(
                location=f"./README.md",
                format=EntryFormat.MARKDOWN,
                master=False,
            ),
        },
        "cc-by": {
            "path": MODEL_BASE_PATH.parent / "cc-by-sa-4.0.txt",
            "entry": ManifestEntry(
                location=f"./cc-by-sa-4.0.txt",
                format=EntryFormat.TXT,
                master=False,
            ),
        },
        "mit": {
            "path": MODEL_BASE_PATH.parent / "mit.txt",
            "entry": ManifestEntry(
                location=f"./mit.txt",
                format=EntryFormat.TXT,
                master=False,
            ),
        },
        "model_figure": {
            "path": MODEL_BASE_PATH.parent / "figures" / "glimepiride_model.png",
            "entry": ManifestEntry(
                location=f"./figures/glimepiride_model.png",
                format=EntryFormat.PNG,
                master=False,
            ),
        },
    }
    if create_tissues:
        for model in [
            model_kidney,
            model_liver,
            model_intestine,
            model_body,
        ]:
            factory_results = create_model(
                model=model,
                filepath=model_output_dir / f"{model.sid}.xml", sbml_level=3, sbml_version=2
            )
            sbml_path = factory_results.sbml_path
            results[model.sid] = {
                "path": sbml_path,
                "entry": ManifestEntry(
                    location=f"./models/{sbml_path.name}",
                    format=EntryFormat.SBML_L3V2,
                    master=False,
                ),
            }

            # create differential equations
            md_path = model_output_dir / f"{model.sid}.md"
            ode_factory = odefac.SBML2ODE.from_file(sbml_file=sbml_path)
            ode_factory.to_markdown(md_file=md_path)
            results[f"{model.sid}_md"] = {
                "path": md_path,
                "entry": ManifestEntry(
                    location=f"./models/{md_path.name}",
                    format=EntryFormat.MARKDOWN,
                    master=False,
                ),
            }

    # create whole-body model
    sbml_path = results["glimepiride_body"]["path"]
    sbml_path_flat = model_output_dir / f"{model_body.sid}_flat.xml"
    flatten_sbml(sbml_path, sbml_flat_path=sbml_path_flat)

    results["glimepiride_body_flat"] = {
        "path": sbml_path_flat,
        "entry": ManifestEntry(
            location=f"./models/{sbml_path_flat.name}",
            format=EntryFormat.SBML_L3V2,
            master=True,
        ),
    }

    # create differential equations
    md_path = model_output_dir / f"{model_body.sid}_flat.md"
    ode_factory = odefac.SBML2ODE.from_file(sbml_file=sbml_path_flat)
    ode_factory.to_markdown(md_file=md_path)
    results[f"{model_body.sid}_flat_md"] = {
        "path": md_path,
        "entry": ManifestEntry(
            location=f"./models/{md_path.name}",
            format=EntryFormat.MARKDOWN,
            master=False,
        ),
    }

    # create omex
    omex = Omex()
    for info in results.values():
        omex.add_entry(entry_path=info["path"], entry=info["entry"])
    omex.to_omex(omex_path=model_output_dir / "glimepiride_model.omex")

    console.print(omex.manifest.model_dump())

    return results


if __name__ == "__main__":
    from pkdb_models.models.glimepiride import MODEL_BASE_PATH
    results = create_models(model_output_dir=MODEL_BASE_PATH, create_tissues=True)

    # visualize_sbml(sbml_path=fac_result.sbml_path)
    for k, key in enumerate(results):
        info = results[key]
        path: Path = info["path"]
        if path.suffix != ".xml":
            continue

        console.print(path)
        visualize_sbml(sbml_path=path, delete_session=(k == 0))
