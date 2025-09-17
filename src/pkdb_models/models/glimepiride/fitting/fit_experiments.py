"""Parameter fit problems for glimepiride."""
from typing import Dict, List
from sbmlsim.fit.helpers import f_fitexp, filter_empty
from sbmlutils.console import console
from sbmlutils.log import get_logger

from sbmlsim.fit import FitExperiment, FitMapping

from pkdb_models.models.glimepiride import GLIMEPIRIDE_PATH, DATA_PATHS
from pkdb_models.models.glimepiride.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health,
    Fasting, GlimepirideMappingMetaData, Coadministration, Genotype
)
from pkdb_models.models.glimepiride.experiments.studies import *


logger = get_logger(__name__)


# --- Filters ---
def filter_baseline(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Return baseline experiments/mappings for reference model."""

    metadata: GlimepirideMappingMetaData = fit_mapping.metadata

    # only PO and IV (no SL, MU, RE)
    if metadata.route not in {Route.PO, Route.IV}:
        return False

    # filter coadminstration
    if metadata.coadministration != Coadministration.NONE:
        return False

    # filter health (no renal, cardiac impairment, ...)
    if metadata.health not in {Health.HEALTHY, Health.T2DM, Health.HYPERTENSION}:
        return False

    # filter multiple dosing (only single dosing)
    if metadata.dosing == Dosing.MULTIPLE:
        return False

    # only fasted subjects
    if metadata.fasting not in {Fasting.FASTED, Fasting.NR}:
        return False

    # remove outliers
    if metadata.outlier is True:
        return False

    return True


def filter_cyp2c9_genotype(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Return control experiments/mappings."""

    metadata: GlimepirideMappingMetaData = fit_mapping.metadata

    # filter genotypes
    return metadata.genotype in {Genotype.CYP2C9_1_1, Genotype.NR}


# --- Fit experiments ---
f_fitexp_kwargs = dict(
    experiment_classes  = [
        Ahmed2016,
        Badian1994,
        Badian1996,
        Choi2014,
        FDA1995,
        Helmy2013,
        Kasichayanula2011c,
        Kim2017,
        Lee2012,
        Lehr1990,
        Liu2010,
        Malerczyk1994,
        Matsuki2007,
        Niemi2002,
        Ratheiser1993,
        Rosenkranz1996a,
        Shukla2004,
        Suzuki2006,
        Wang2005,
        Yoo2011,
    ],
    base_path=GLIMEPIRIDE_PATH,
    data_path=DATA_PATHS,
)
# --- Experiment classes ---

def f_fitexp_all():
    """All data."""
    return f_fitexp(metadata_filters=filter_empty, **f_fitexp_kwargs)

def f_fitexp_cyp2c9() -> Dict[str, List[FitExperiment]]:
    """Control data with all CYP2C9 data."""
    return f_fitexp(metadata_filters=[filter_baseline], **f_fitexp_kwargs)

def f_fitexp_control() -> Dict[str, List[FitExperiment]]:
    """Control data with CYP2C9 wildtype."""
    return f_fitexp(metadata_filters=[filter_baseline, filter_cyp2c9_genotype], **f_fitexp_kwargs)




if __name__ == "__main__":
    """Test construction of FitExperiments."""

    for f in [
        f_fitexp_all,
        f_fitexp_cyp2c9,
        f_fitexp_control,

    ]:
        console.rule(style="white")
        console.print(f"{f.__name__}")
        fitexp = f()
