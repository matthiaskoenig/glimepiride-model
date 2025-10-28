"""Run all simulation experiments."""
import shutil
from typing import List, Dict
from sbmlutils import log
from sbmlutils.console import console
from sbmlsim.plot import Figure
import pkdb_models.models.glimepiride as glimepiride
from pkdb_models.models.glimepiride.helpers import run_experiments
from pkdb_models.models.glimepiride.experiments.studies import *
from pkdb_models.models.glimepiride.experiments.misc import *
from pkdb_models.models.glimepiride.experiments.scans import *
from matplotlib import pyplot as plt


Figure.legend_fontsize = 11
Figure.axes_labelsize = 16
Figure.xtick_labelsize = 12
Figure.ytick_labelsize = 12
# Figure.fig_dpi = 600

logger = log.get_logger(__name__)

def run_simulation_experiments(selected: str = None, specific_experiments: List[str] = None, list_only: bool = False) -> None:
    """Run simulation experiments."""

    experiments = {
        "studies": [
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
        "dose_dependency": [
            Malerczyk1994,
            Helmy2013,
            Ratheiser1993,
        ],
        "bodyweight": [
            Shukla2004,
        ],
        "cyp2c9": [
            Lee2012,
            Niemi2002,
            Suzuki2006,
            Wang2005,
            Yoo2011,
        ],
        "hepatic_impairment": [
        ],
        "renal_impairment": [
            Rosenkranz1996a,
        ],
        "misc": [
            BodyweightExperiment,
            DoseDependencyExperiment,
            GeneticVariantExperiment,
            HepaticRenalImpairmentExperiment,
        ],
        "scans": [
            GlimepirideBodyweightScan,
            GlimepirideCirrhosisScan,
            GlimepirideCrClScan,
            GlimepirideDoseScan,
            GlimepirideParameterScan,
        ]

    }
    experiments["all"] = experiments["studies"] + experiments["misc"] + experiments["scans"]

    # experiment name collector for run_glimepiride.py
    if list_only:
        console.rule("[bold cyan]Available Simulation Experiments[/bold cyan]", style="cyan")
        console.print("\n[bold]You can use these group names:[/bold]")
        console.print(f"  {', '.join([g for g in experiments.keys()])}")
        console.print("\n[bold]Or these individual experiment names:[/bold]")
        for group_name in ["studies", "misc", "scans"]:
            if group_name in experiments and experiments[group_name]:
                console.print(f"\n[yellow]{group_name}:[/yellow]")
                for exp in experiments[group_name]:
                    console.print(f"{exp.__name__}")
        console.print("\n[dim]Use '--experiments' with comma-separated names to run specific experiments.[/dim]")
        console.print('[dim]Example: run_glimepiride --action simulate --experiments "cyp2c9,Ahmed2016"[/dim]')
        console.print('[dim]Or use "all" to run all experiments: run_glimepiride --action simulate --experiments all[/dim]\n')
        return

    # Determine which experiments to run
    if specific_experiments:
        # dictionary of all available experiments
        all_available_exp = {}
        for category, exp_list in experiments.items():
            if category != "all":
                for exp in exp_list:
                    all_available_exp[exp.__name__] = exp

        # Validate and collect requested experiments
        experiment_classes = []
        not_found = []

        for exp_name in specific_experiments:
            # group name
            if exp_name in experiments:
                experiment_classes.extend(experiments[exp_name])
            elif exp_name in all_available_exp:
                # individual experiment
                experiment_classes.append(all_available_exp[exp_name])
            else:
                not_found.append(exp_name)

        # Report any experiments that weren't found
        if not_found:
            console.rule(style="red bold")
            console.print(f"[red]Warning: The following experiments were not found: {', '.join(not_found)}[/red]")
            console.rule(style="red bold")
        if not experiment_classes:
            console.rule(style="red bold")
            console.print("[red]Error: No valid experiments to run![/red]")
            console.rule(style="red bold")
            return

        # output directory for custom selection
        output_dir = glimepiride.RESULTS_PATH_SIMULATION / "custom_selection"

    else:
        if not selected:
            console.print("\n[red bold]Error: No experiments specified![/red bold]")
            console.print("[yellow]For example, use selected='all' or selected='studies' or specific_experiments=[...][/yellow]\n")
            return

        if selected not in experiments:
            console.rule(style="red bold")
            console.print(
                f"[red]Error: Unknown group '{selected}'. Valid groups: {', '.join(experiments.keys())}[/red]")
            console.rule(style="red bold")
            return

        experiment_classes = experiments[selected]
        output_dir = glimepiride.RESULTS_PATH_SIMULATION / selected

    run_experiments(
        experiment_classes=experiment_classes,
        output_dir=output_dir,
    )
    # collect figures
    plt.close('all')
    import gc
    gc.collect()  # Ensures file objects are fully garbage collected
    figures_dir = output_dir / "_figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    for f in output_dir.glob("**/*.png"):
        if f.parent == figures_dir:
            continue
        try:
            shutil.copy2(f, figures_dir / f.name)
            # print(f"file {f.name} in {f.parent} copied to {figures_dir / f.name}")
        except Exception as err:
            print(f"file {f.name} in {f.parent} fails, skipping. Error: {err}")
    console.print(f"Figures copied to: file://{figures_dir}", style="info")


if __name__ == "__main__":
    """
    # Run predefined groups
    run_simulation_experiments(selected="studies")      # "all", "studies", "misc", etc.

    # Run specific experiments
    run_simulation_experiments(specific_experiments=["Ahmed2016", "Badian1994"])

    # Mix groups and individual experiments
    run_simulation_experiments(specific_experiments=["cyp2c9", "GlimepirideDoseScan", "Choi2014"])
    """

    run_simulation_experiments(selected="all")
    # run_simulation_experiments(specific_experiments=["Ahmed2016", "Badian1994", "misc"])
