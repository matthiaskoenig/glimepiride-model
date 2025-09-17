"""Runs all standalone CYP2C9 analyses scripts and collects their figures."""

import subprocess
import sys
from pkdb_models.models.glimepiride import RESULTS_PATH_ANALYSES, GLIMEPIRIDE_PATH

CYP2C9_PATH = GLIMEPIRIDE_PATH / "experiments" / "cyp2c9"

cyp2c9_scripts = [
    CYP2C9_PATH / "analyze_cyp2c9_activity.py",
    CYP2C9_PATH / "plot_cyp2c9_pk.py",
    CYP2C9_PATH / "population_sampling.py",
]


if __name__ == "__main__":
    RESULTS_PATH_ANALYSES.mkdir(parents=True, exist_ok=True)
    for script_path in cyp2c9_scripts:
        subprocess.run(
            [sys.executable, script_path.name],
            cwd=script_path.parent,
            check=True,
        )