"""Tool to run the glimepiride model factory, simulation experiments, or standalone CYP2C9 analysis scripts."""

import sys
import os
import subprocess
from enum import Enum
import optparse
from pathlib import Path
from pkdb_models.models.glimepiride import GLIMEPIRIDE_PATH
from pkdb_models.models.glimepiride.simulations import run_simulation_experiments
from pkdb_models.models.glimepiride.run_cyp2c9_analyses import cyp2c9_scripts
from sbmlutils.console import console

FACTORY_SCRIPT_PATH = GLIMEPIRIDE_PATH / "models" / "factory.py"

class Action(str, Enum):
    # simulations
    SIMULATE = "simulate"
    LIST_EXPERIMENTS = "list_experiments"
    # CYP2C9 scripts
    ANALYZE_CYP2C9 = "analyze_cyp2c9"
    LIST_CYP2C9 = "list_cyp2c9"
    # factory
    FACTORY = "factory"
    # run all
    ALL = "all"

# Available CYP2C9 scripts
CYP2C9_SCRIPT_MAP = {
    "activity": "analyze_cyp2c9_activity.py",
    "pk": "plot_cyp2c9_pk.py",
    "population": "population_sampling.py",
}

def _setup_custom_results_paths(results_dir: str):
    """Override default results paths with custom directory for custom figure output directory."""
    import pkdb_models.models.glimepiride as glimepiride
    custom_path = Path(results_dir).resolve()
    # Create directory if it doesn't exist
    custom_path.mkdir(parents=True, exist_ok=True)
    # Override the module paths
    glimepiride.RESULTS_PATH = custom_path
    glimepiride.RESULTS_PATH_SIMULATION = custom_path / "simulation"
    glimepiride.RESULTS_PATH_ANALYSES = custom_path / "cyp2c9_analyses"

    console.print(f"Figure output directory set to: [cyan]{custom_path}[/cyan]")
    return custom_path

def _get_current_results_path():
    """Get the current results path (default or custom)."""
    from pkdb_models.models.glimepiride import RESULTS_PATH
    return RESULTS_PATH

def _run_factory():
    """Executes the model factory script."""
    console.rule("[bold cyan]Running Model Factory[/bold cyan]", style="cyan")
    subprocess.run(
        [sys.executable, str(FACTORY_SCRIPT_PATH)],
        cwd=FACTORY_SCRIPT_PATH.parent,
        check=True,
    )
    console.print("[bold green]Factory finished.[/bold green]")

def _list_cyp2c9_scripts():
    """Lists available CYP2C9 analysis scripts."""
    console.rule("[bold cyan]Available CYP2C9 Analysis Scripts[/bold cyan]", style="cyan")
    console.print("\n[bold]Available scripts:[/bold]")
    for key, script_name in CYP2C9_SCRIPT_MAP.items():
        console.print(f"  • [cyan]{key}[/cyan]: {script_name}")
    console.print("\n[dim]Use '--scripts' or '-s' with comma-separated names to run specific scripts.[/dim]")
    console.print('[dim]Example: run_glimepiride --action analyze_cyp2c9 --scripts "activity,pk"[/dim]')
    console.print('[dim]Or use "all" to run all scripts: run_glimepiride --action analyze_cyp2c9 --scripts all[/dim]\n')


def _run_cyp2c9_analysis(selected_scripts=None):
    """Runs the standalone CYP2C9 analysis scripts.
    selected_scripts: List of script keys to run, or None to run all
    """
    results_path = _get_current_results_path()
    console.rule("[bold cyan]Running CYP2C9 Analysis[/bold cyan]", style="cyan")
    console.print(f"[cyan]Results will be saved to: {results_path / 'cyp2c9_analyses'}[/cyan]")

    if selected_scripts == ["all"]:
        # Run all scripts
        scripts_to_run = cyp2c9_scripts
        console.print("[cyan]Running all CYP2C9 analysis scripts...[/cyan]\n")
    elif selected_scripts is None:
        # Run all scripts
        scripts_to_run = cyp2c9_scripts
        console.print("[cyan]Default: Running all CYP2C9 analysis scripts...[/cyan]\n")
    else:
        # Run selected scripts
        scripts_to_run = []
        for script_key in selected_scripts:
            script_key = script_key.strip().lower()
            if script_key in CYP2C9_SCRIPT_MAP:
                script_name = CYP2C9_SCRIPT_MAP[script_key]
                for script_path in cyp2c9_scripts:
                    if script_path.name == script_name:
                        scripts_to_run.append(script_path)
                        break
            else:
                console.print(f"[yellow]Warning: Unknown script key '{script_key}'. Skipping.[/yellow]")
                console.print(f"[dim]Available keys: {', '.join(CYP2C9_SCRIPT_MAP.keys())}[/dim]")

        if not scripts_to_run:
            console.print("[red]No valid scripts selected. Use '--action list_cyp2c9' to see available options.[/red]")
            return

    env = {**os.environ, "GLIMEPIRIDE_RESULTS_PATH": str(_get_current_results_path())}
    for script_path in scripts_to_run:
        console.print(f" -> Running: [bold cyan]{script_path.name}[/bold cyan]")
        subprocess.run(
            [sys.executable, script_path.name],
            cwd=script_path.parent,
            check=True,
            env=env,
        )

    console.print(f"\n[bold green]CYP2C9 analysis finished - ({len(scripts_to_run)} script(s) completed).[/bold green]")
    console.print(f"[bold green]Results saved to: {results_path / 'cyp2c9_analyses'}[/bold green]")


def main() -> None:

    parser = optparse.OptionParser()
    parser.add_option(
        "-a", "--action",
        dest="action",
        help=f"The main action to perform. Choices: {[a.value for a in Action]} (required)",
    )
    parser.add_option(
        "-r", "--results-dir",
        dest="results_dir",
        help="Optional: Custom directory to save all results/figures (default: ./results)",
    )
    parser.add_option(
        "-s", "--scripts",
        dest="scripts",
        help="Comma-separated list of CYP2C9 analysis scripts to run (for '--action analyze_cyp2c9'). "
             "Use '--action list_cyp2c9' to see available options.",
    )
    parser.add_option(
        "-e", "--experiments",
        dest="experiments",
        help="Comma-separated list of simulation experiments and/or groups (for '--action simulate'). "
             "Use '--action list_experiments' to see all available options.",
    )

    console.rule(style="white")
    console.print("GLIMEPIRIDE PBPK MODEL")
    console.rule(style="white")

    options, args = parser.parse_args()

    def _parser_message(text: str) -> None:
        console.print(f"[bold red]Error: {text}[/bold red]")
        parser.print_help()
        console.rule(style="red")
        sys.exit(1)

    if not options.action:
        _parser_message("Required argument '--action' is missing.")

    try:
        action = Action(options.action.lower())
    except ValueError:
        _parser_message(f"Invalid action '{options.action}'. Please choose from {[a.value for a in Action]}.")

    # Setup custom results directory if provided
    if options.results_dir:
        _setup_custom_results_paths(options.results_dir)
    else:
        import pkdb_models.models.glimepiride as glimepiride
        default_path = _get_current_results_path()
        if options.action and options.action.lower() == "factory":
            console.print(f"[cyan]Factory will use model base path: {glimepiride.MODEL_BASE_PATH}[/cyan]")
        else:
            console.print(f"[cyan]Using figure results directory: {default_path}[/cyan]")

    if action == Action.FACTORY:
        _run_factory()

    elif action == Action.LIST_EXPERIMENTS:
        run_simulation_experiments(list_only=True)

    elif action == Action.SIMULATE:

        if not options.experiments:
            _parser_message("For '--action simulate', the '--experiments' argument is required.")
        exp_list = [e.strip() for e in options.experiments.split(",")]
        results_path = _get_current_results_path()
        console.rule("[bold cyan]Running Simulations[/bold cyan]", style="cyan")
        run_simulation_experiments(specific_experiments=exp_list)
        console.print("[bold green]Simulations finished.[/bold green]")
        console.print(f"[bold green]Results saved to: {results_path / 'simulation'}[/bold green]")

    elif action == Action.LIST_CYP2C9:
        _list_cyp2c9_scripts()

    elif action == Action.ANALYZE_CYP2C9:
        if options.scripts:
            script_list = [s.strip() for s in options.scripts.split(",")]
            _run_cyp2c9_analysis(selected_scripts=script_list)
        else:
            _run_cyp2c9_analysis()

    elif action == Action.ALL:
        console.print("[bold magenta]Running: Factory, all Simulations, and all CYP2C9 Analysis[/bold magenta]")
        _run_factory()                              # Run factory
        run_simulation_experiments(selected="all")  # Run all simulations
        _run_cyp2c9_analysis()                      # Run the standalone CYP2C9 analyses
        console.print("\n[bold green]All scripts completed successfully![/bold green]")

    console.rule(style="white")


if __name__ == "__main__":
    """
    This script is intended to be run from the command line as 'run_glimepiride'.

    -------------------------------------------------
    Usage Examples (to be run in your terminal):
    -------------------------------------------------

    1. Help:
       Shows all available options and actions.
       $ run_glimepiride --help


    2. Run the Model Factory:
       Generates all SBML model files.
       $ run_glimepiride --action factory


    3. Run Simulations:
       List available experiments:
       $ run_glimepiride --action list_experiments

       Run specific experiments:
       $ run_glimepiride --action simulate --experiments "Ahmed2016,Badian1994"

       Run all experiments:
       $ run_glimepiride --action simulate --experiments all


    4. Run CYP2C9 Analysis:
       List available analysis scripts:
       $ run_glimepiride --action list_cyp2c9

       Run all CYP2C9 analyses:
       $ run_glimepiride --action analyze_cyp2c9


    5. Run Everything:
       Runs factory, all simulations, and all analyses in sequence.
       $ run_glimepiride --action all
       
       With custom results directory for figures:
       $ run_glimepiride --action all --results-dir '/path/to/my/results'
    """
    main()