import click
import os
from breakdown.auto_detect_isomerization import AutoDetectIsomerization
from breakdown.generate_charts_stats import GenerateChartsStats, inferMoleculeFromLogs
from breakdown.molecule_transition_graph import MoleculeTransitionGraph
from breakdown.optimize_isomers import IsomerOptimizer
from breakdown.energies_table import EnergiesTable
from breakdown.functional_group_comparison import generateComparison
from breakdown.elemental_stability_chart import generateElementalComparison
from breakdown.oxygen_fate_pie_chart import generateOxygenFates
from breakdown.run_summary import loadBlock
from breakdown.event_survival_chart import generateEventComparison
from breakdown.pipeline import (CONFIG_FILE, ConfigError, analyseFiles, defaultConfig, formatConfig, listMolFiles,
                                 loadConfig, logPathFor, runPipeline)

def packageVersion():
    """Installed version of the package, or unknown when it is run from a source folder."""
    try:
        from importlib.metadata import version
        return version('breakdown')
    except Exception:
        return 'unknown'

@click.group()
@click.version_option(version=packageVersion(), prog_name='breakdown')
def main_group():
    pass

def reportPlot(path):
    """Tell the user where a plot was saved, or that there was nothing to plot."""
    if path:
        click.echo(f"Saved {path} and a png of it next to it.")
    else:
        raise click.ClickException("Nothing to plot, none of the folders has the needed data.")

def analysisOptions(function):
    """Options shared by the commands that analyse .mol files."""
    options = [
        click.option("--molecule_size", "-s", default=None, type=int,
                     help="Number of atoms in starting molecule. Read from the file if not given."),
        click.option("--elements", "-e", multiple=True, default=(), type=str,
                     help="Elements present in the starting molecule, passed separately, like "
                          "-e C -e H -e O. Read from the file if not given."),
        click.option("--formation_distance", "-fd", default=1.63, type=float,
                     help="Distance in angstroms where a new bond will form."),
        click.option("--break_distance", "-bd", default=2.5, type=float,
                     help="Distance in angstroms where a bond will dissociate."),
        click.option("--min_lifetime", "-ml", default=1, type=int,
                     help="Number of steps a new molecule state must last to be reported, shorter "
                          "lived states are ignored as noise. 1 reports everything."),
    ]
    for option in reversed(options):
        function = option(function)
    return function

@main_group.command()
@click.argument("folder")
@click.argument("file_name")
@analysisOptions
def analyseDemon(folder, file_name, molecule_size, elements,
                 formation_distance, break_distance, min_lifetime):
    """Convert one deMon .mol trajectory file into a _log.txt file of SMILES strings.

    Args:
        folder (str): Folder containing deMon.mol files.
        file_name (str): Name of the file to be analyzed.
        molecule_size (int): Number of atoms in starting molecule.
        elements (list): List of elements present in the starting molecule.
        formation_distance (float): Distance in angstroms where a new bond will form.
        break_distance (float): Distance in angstroms where a bond will dissociate.
        min_lifetime (int): Number of steps a new molecule state must last to be reported.
    """
    autoDetectIsomerization = AutoDetectIsomerization(
        folder, file_name, molecule_size, elements, formation_distance, break_distance,
        min_lifetime)
    try:
        autoDetectIsomerization.autoAnalyseFile()
    except ValueError as error:
        raise click.ClickException(str(error))

@main_group.command()
@click.argument("folder")
@click.option("--jobs", "-j", default=os.cpu_count() or 1, type=int,
              help="Number of files analysed at the same time.")
@click.option("--pattern", "-p", default="*.mol", type=str,
              help="Pattern for the names of the files to analyse.")
@click.option("--skip_existing/--no-skip_existing", default=False,
              help="Skip files that already have a _log.txt file.")
@analysisOptions
def analyseDemonFolder(folder, jobs, pattern, skip_existing, molecule_size, elements,
                       formation_distance, break_distance, min_lifetime):
    """Convert every deMon .mol trajectory file in a folder into a _log.txt file,
    analysing several files at once.

    Args:
        folder (str): Folder containing deMon.mol files.
        jobs (int): Number of files analysed at the same time.
        pattern (str): Pattern for the names of the files to analyse.
        skip_existing (bool): Skip files that already have a log.
    """
    fileNames = listMolFiles(folder, pattern)
    if skip_existing:
        fileNames = [name for name in fileNames if not os.path.exists(logPathFor(folder, name))]
    if not fileNames:
        click.echo(f"No files matching {pattern} to analyse in {folder}.")
        return
    options = (molecule_size, list(elements), formation_distance, break_distance, min_lifetime)
    failures = analyseFiles(folder, fileNames, options, jobs)
    if failures:
        click.echo(f"{len(failures)} files failed: {', '.join(failures)}")
        raise SystemExit(1)

@main_group.command()
@click.argument("folder")
@click.argument("chart_name")
@click.argument("molecule_name")
@click.option("--num_carbons", "-s", default=None, type=int,
              help="Number of carbon atoms in starting molecule. Read from the logs if not given.")
@click.option("--elements", "-e", default=None,
              help="Elements that make up the molecule, comma-separated. Read from the logs if not given.")
@click.option("--func_group", "-f", default=None, type=str,
              help="Original functional group to track, e.g. ketone, hydroxyl, thiol")
@click.option("--legend_carbons", "-lc", default=None, type=str,
              help="Carbon counts to show in legend, comma-separated (e.g. 14,12,10). Default is all.")
@click.option("--time_per_step_fs", "-t", default=5.0, type=float,
              help="Simulated time between steps in the logs, in femtoseconds. The default is "
                   "0.5 fs time steps with output every 10 steps.")
@click.option('--draw/--no-draw', default=False, help="Draw each isomer.")
@click.option('--coordinates/--no-coordinates', default=False, help="Make coordinate files for each isomer.")
def makeChartsAndStats(folder, chart_name, molecule_name, num_carbons, elements,
                       func_group, legend_carbons, time_per_step_fs, draw, coordinates):
    """Combine all _log.txt files in a folder into full_log.txt with statistics, and make plots."""
    legend_list = None
    if legend_carbons:
        try:
            legend_list = [int(c.strip()) for c in legend_carbons.split(",")]
        except ValueError:
            click.echo("Error: --legend_carbons must be a comma-separated list of integers.")
            return
    element_list = elements.split(",") if elements else None
    if num_carbons is None or element_list is None:
        foundCarbons, foundElements = inferMoleculeFromLogs(folder)
        num_carbons = foundCarbons if num_carbons is None else num_carbons
        element_list = foundElements if element_list is None else element_list
        click.echo(f"Using {num_carbons} carbons and elements {','.join(element_list)}.")
    generateChartsStats = GenerateChartsStats(
        folder, chart_name, molecule_name, num_carbons, element_list,
        func_group, draw, coordinates, legendCarbons=legend_list,
        timePerStepFs=time_per_step_fs)
    generateChartsStats.autoAnalyseInFolder()

@main_group.command()
@click.argument("folder")
@click.argument("start_state")
@click.option("--num_top_states", "-n", default=10, type=int,
              help="Number of top states to include in the flow chart.")
@click.option("--num_scattered", "-ns", default=5, type=int,
                help="Number of scattered states to include in the flow chart.")
@click.option("--max_depth_best", "-mb", default=5, type=int,
              help="Depth to which the best path will be looked for.")
@click.option("--max_depth_possible", "-mp", default=20, type=int,
              help="Depth to which the any path will be looked for.")
@click.option('--direct/--no-direct', default=False,
              help="Only do the start state and directly connected states.")
@click.option('--states', '-s', multiple=True, default=[],
              help="List of states that must be visited.")
def makeFlowCharts(folder, start_state, num_top_states, num_scattered, max_depth_best, max_depth_possible,
                   direct, states):
    """Make flow charts of the transitions between molecule states, from full_log.txt.

    Args:
        folder (str): Folder containing full_log.txt.
        start_state (str): Name of start state isomer.
        num_top_states (int): Number of top states to include in the flow chart.
        num_scattered (int): Number of scattered states to include in the flow chart.
        max_depth_best (int): Depth to which the best path will be looked for.
        max_depth_possible (int): Depth to which the any path will be looked for.
        direct (bool): Only do the start state and directly connected states.
        states (list): List of states that must be visited.
    """
    sections = [loadBlock(folder, key) for key in
                ('isomerNames', 'isomerStats', 'isomerizationStats', 'netIsomerizationStats')]
    if None in sections:
        raise click.ClickException(f'{folder} has no full_log.txt with the needed statistics, '
                                   'run makechartsandstats first.')
    isomerNames, isomerStats, isomerizationStats, netTransitions = sections
    if start_state not in isomerNames.values():
        raise click.ClickException(f'{start_state} is not a state name in {folder}/full_log.txt, '
                                   'the names are listed under "Automatically generated molecule names."')
    moleculeTransitionGraph = MoleculeTransitionGraph(folder, isomerNames, isomerizationStats, isomerStats,
                                                      netTransitions, start_state, num_top_states, num_scattered,
                                                      max_depth_best, max_depth_possible, states)
    if direct:
        moleculeTransitionGraph.makeDirectTransitionGraph()
    else:
        moleculeTransitionGraph.makeTransitionGraph()

@main_group.command()
@click.argument("folder")
@click.argument("isomer_name")
def optimizeIsomer(folder, isomer_name):
    """Main function to run isomer optimizer on all xyz files in a folder.

    Args:
        folder (str): Folder containing .xyz files.
        isomer_name (str): Name of the isomer to optimize.
    """
    isomerOptimizer = IsomerOptimizer(folder, isomer_name)
    isomerOptimizer.optimizeIsomerGeometry()

@main_group.command()
@click.argument("folder")
def makeEnergiesTable(folder):
    """Main function to create energies table from all .out files in a folder.

    Args:
        folder (str): Folder containing .out files.
    """
    energiesTable = EnergiesTable(folder)
    energiesTable.processOutFiles()
    energiesTable.saveToCsv()

@main_group.command()
@click.option("--directory", "-d", required=True, type=str,
              help="Base directory containing the molecule folders and output folder.")
@click.option("--folders", "-f", required=True, type=str,
              help="Semicolon-separated names of molecule folders (e.g., mol1; mol2)")
@click.option("--names", "-n", required=True, type=str,
              help="Semicolon-separated abbreviations for the molecules (e.g., ANT; PHEN)")
@click.option("--groups", "-g", required=True, type=str,
              help="Semicolon-separated functional groups to track (e.g., ketone; hydroxyl)")
@click.option("--title", "-t", default="Group_Comparison", type=str,
              help="Title of the output chart and filename.")
@click.option("--output_dir", "-o", default="comparisons", type=str,
              help="Folder name inside the base directory where the plot will be saved.")
@click.option("--time_per_step_fs", default=None, type=float,
              help="Simulated time between steps in the logs, in femtoseconds. Read from run_summary.json if not given.")
def compareFunctionalGroups(directory, folders, names, groups, title, output_dir, time_per_step_fs):
    """Compares functional group survival across different molecule simulations.

    Uses semicolons (;) as separators. Folders and output_dir are relative to --directory.
    """
    folderList = [os.path.join(directory, f.strip()) for f in folders.split(";")]
    nameList = [n.strip() for n in names.split(";")]
    groupList = [g.strip() for g in groups.split(";")]
    
    finalOutputDir = os.path.join(directory, output_dir)

    if not (len(folderList) == len(nameList) == len(groupList)):
        click.echo(f"Error: Mismatch in counts.")
        click.echo(f"Folders found: {len(folderList)}")
        click.echo(f"Names found: {len(nameList)}")
        click.echo(f"Groups found: {len(groupList)}")
        return
    click.echo(f"Base Directory: {directory}")
    click.echo(f"Generating comparison chart: {title}...")
    path = generateComparison(folderList, nameList, groupList, title, finalOutputDir, time_per_step_fs)
    reportPlot(path)

@main_group.command()
@click.option("--directory", "-d", required=True, type=str,
              help="Base directory containing the molecule folders and output folder.")
@click.option("--folders", "-f", required=True, type=str,
              help="Semicolon-separated names of molecule folders (e.g., mol1; mol2)")
@click.option("--names", "-n", required=True, type=str,
              help="Semicolon-separated abbreviations for the molecules (e.g., ANT; PHEN)")
@click.option("--title", "-t", default="Elemental_Stability", type=str,
              help="Title of the output chart and filename.")
@click.option("--output_dir", "-o", default="comparisons", type=str,
              help="Folder name inside the base directory where the plot will be saved.")
@click.option("--heteroatom", default="O", type=click.Choice(["O", "S"]),
              help="Heteroatom whose retention is compared with carbon.")
@click.option("--time_per_step_fs", default=None, type=float,
              help="Simulated time between steps in the logs, in femtoseconds. Read from run_summary.json if not given.")
def compareElementalStability(directory, folders, names, title, output_dir, heteroatom, time_per_step_fs):
    """Compares Carbon and Oxygen (or Sulfur) retention across different molecule simulations.

    Creates a dual-subplot figure showing runs that maintain original atom counts.
    Uses semicolons (;) as separators. Folders and output_dir are relative to --directory.
    """
    folderList = [os.path.join(directory, f.strip()) for f in folders.split(";")]
    nameList = [n.strip() for n in names.split(";")]
    finalOutputDir = os.path.join(directory, output_dir)

    if not (len(folderList) == len(nameList)):
        click.echo(f"Error: Mismatch in counts.")
        click.echo(f"Folders found: {len(folderList)}")
        click.echo(f"Names found: {len(nameList)}")
        return
    
    click.echo(f"Base Directory: {directory}")
    click.echo(f"Generating elemental stability chart: {title}...")
    
    path = generateElementalComparison(folderList, nameList, title, finalOutputDir, heteroatom, time_per_step_fs)
    
    reportPlot(path)


@main_group.command()
@click.option("--directory", "-d", required=True, type=str,
              help="Base directory containing the molecule folders and output folder.")
@click.option("--folders", "-f", required=True, type=str,
              help="Semicolon-separated names of molecule folders (e.g., mol1; mol2)")
@click.option("--names", "-n", required=True, type=str,
              help="Semicolon-separated abbreviations for the molecules (e.g., ANT; PHEN)")
@click.option("--event", "-e", default="fragmentation",
              type=click.Choice(["anyChange", "isomerization", "fragmentation", "recombination"]),
              help="Which first event to compare.")
@click.option("--title", "-t", default="First_Event", type=str,
              help="Title of the output chart and filename.")
@click.option("--output_dir", "-o", default="comparisons", type=str,
              help="Folder name inside the base directory where the plot will be saved.")
@click.option("--time_per_step_fs", default=None, type=float,
              help="Simulated time between steps in the logs, in femtoseconds. Read from run_summary.json if not given.")
def compareFirstEvents(directory, folders, names, event, title, output_dir, time_per_step_fs):
    """Compares how long molecules survive until their first isomerization or fragmentation.

    Draws survival curves (fraction of runs without the event yet) with 95% bands, using the
    firstEvents made by makechartsandstats. Uses semicolons (;) as separators.
    """
    folderList = [os.path.join(directory, f.strip()) for f in folders.split(";")]
    nameList = [n.strip() for n in names.split(";")]
    if len(folderList) != len(nameList):
        raise click.ClickException("The number of folders and names differ.")
    try:
        path = generateEventComparison(folderList, nameList, event, title,
                                       os.path.join(directory, output_dir), time_per_step_fs)
    except ValueError as error:
        raise click.ClickException(str(error))
    reportPlot(path)


@main_group.command()
@click.option("--directory", "-d", required=True, type=str,
              help="Base directory containing the molecule folders.")
@click.option("--folders", "-f", required=True, type=str,
              help="Semicolon-separated names of molecule folders.")
@click.option("--names", "-n", required=True, type=str,
              help="Semicolon-separated abbreviations for the molecules.")
@click.option("--groups", "-g", required=True, type=str,
              help="Semicolon-separated starting oxygen groups for each folder.")
@click.option("--title", "-t", default="Oxygen_Fate_Analysis", type=str,
              help="Title of the output chart and filename.")
@click.option("--output_dir", "-o", default="comparisons", type=str,
              help="Folder name inside the base directory for the output.")
def plotOxygenFates(directory, folders, names, groups, title, output_dir):
    """Generates pie charts showing the final state of oxygen for changed runs.
    
    Filters out the starting group to show where oxygen went (fragments or new groups).
    """
    folderList = [os.path.join(directory, f.strip()) for f in folders.split(";")]
    nameList = [n.strip() for n in names.split(";")]
    groupList = [g.strip() for g in groups.split(";")]
    finalOutputDir = os.path.join(directory, output_dir)

    if not (len(folderList) == len(nameList) == len(groupList)):
        click.echo("Error: Mismatch in folder, name, or group counts.")
        return

    click.echo(f"Analyzing oxygen fate for {len(folderList)} molecules...")
    reportPlot(generateOxygenFates(folderList, nameList, groupList, title, finalOutputDir))


@main_group.command()
@click.argument("folder")
@click.option("--config", "-c", "config_path", default=None, type=str,
              help=f"TOML file with the options. Default: FOLDER/{CONFIG_FILE} if it exists, else the defaults.")
@click.option("--write_config", is_flag=True,
              help=f"Only write FOLDER/{CONFIG_FILE} with every option, its default and what it means, then stop.")
def run(folder, config_path, write_config):
    """Do the whole analysis of a folder: .mol files -> logs -> full_log.txt and plots -> flow charts -> report.html.

    Logs that were already made with the same settings are kept. The options used are written to
    FOLDER/breakdown_used.toml so the analysis can be repeated with --config.
    """
    if not os.path.isdir(folder):
        raise click.ClickException(f"{folder} is not a folder.")
    defaultPath = os.path.join(folder, CONFIG_FILE)
    if write_config:
        if os.path.exists(defaultPath):
            raise click.ClickException(f"{defaultPath} already exists, edit it or delete it first.")
        with open(defaultPath, 'w') as file:
            file.write(f"# Options of 'breakdown run'. Empty strings, 0 and [] mean: read from the files / not used.\n\n"
                       + formatConfig(defaultConfig()))
        click.echo(f"Wrote {defaultPath}, edit it and run: breakdown run {folder}")
        return
    if config_path is None and os.path.exists(defaultPath):
        config_path = defaultPath
    try:
        config = loadConfig(config_path)
        click.echo(f"Options from {config_path}." if config_path else "No configuration file, using the defaults.")
        runPipeline(folder, config)
    except (ConfigError, ValueError) as error:
        raise click.ClickException(str(error))


@main_group.command()
@click.option("--directory", "-d", required=True, type=str,
              help="Base directory containing the folders.")
@click.option("--folders", "-f", required=True, type=str,
              help="Semicolon-separated names of the folders, each with full_log.txt.")
@click.option("--names", "-n", default=None, type=str,
              help="Semicolon-separated short labels for the folders, the folder names if not given.")
@click.option("--largest_fragment/--whole_state", default=False,
              help="Match states by their largest fragment only (C12H8 + H and C12H8 + H2 count as the same).")
@click.option("--top", "-t", default=30, type=int,
              help="Number of rows to print, the csv file has all of them.")
@click.option("--output_dir", "-o", default="comparisons", type=str,
              help="Folder name inside the base directory for the csv file.")
def matchStates(directory, folders, names, largest_fragment, top, output_dir):
    """Show which state names in different folders are the same structure (same connectivity).

    State names are numbered by how common a state is in its own folder, so PHE_C5 at 15 eV and PHE_C5 at 20 eV
    are usually different structures. This lines them up by connectivity key and writes DIR/OUTPUT_DIR/state_matches.csv.
    """
    import csv
    from breakdown.state_keys import matchStates as match
    folderNames = [f.strip() for f in folders.split(";")]
    labels = [n.strip() for n in names.split(";")] if names else folderNames
    if len(labels) != len(folderNames):
        raise click.ClickException(f"{len(folderNames)} folders but {len(labels)} names.")
    try:
        rows = match([os.path.join(directory, f) for f in folderNames], labels, largest_fragment)
    except ValueError as error:
        raise click.ClickException(str(error))

    def cell(states):
        return ', '.join(f'{name} ({runs} runs)' for name, runs, _ in states) or '-'
    outputFolder = os.path.join(directory, output_dir)
    os.makedirs(outputFolder, exist_ok=True)
    csvPath = os.path.join(outputFolder, 'state_matches.csv')
    with open(csvPath, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['key', 'formula'] + [f'{label} {column}' for label in labels
                                              for column in ('states', 'runs', 'steps')])
        for row in rows:
            values = [row['key'], row['formula']]
            for label in labels:
                states = row['states'][label]
                values += [' '.join(name for name, _, _ in states), sum(s[1] for s in states), sum(s[2] for s in states)]
            writer.writerow(values)
    shared = sum(1 for row in rows if all(row['states'][label] for label in labels))
    click.echo(f"{len(rows)} distinct {'largest fragments' if largest_fragment else 'states'}, "
               f"{shared} found in every folder. Most common first:")
    header = ['key', 'formula'] + labels
    table = [header] + [[row['key'], row['formula']] + [cell(row['states'][l]) for l in labels] for row in rows[:top]]
    widths = [max(len(r[i]) for r in table) for i in range(len(header))]
    for line in table:
        click.echo('  '.join(value.ljust(width) for value, width in zip(line, widths)).rstrip())
    click.echo(f"Saved {csvPath}.")


@main_group.command()
@click.argument("folder")
@click.argument("state")
@click.option("--from_state", "-f", default=None, type=str,
              help="Write the transitions from this state into STATE (a multi-frame xyz plus reactant and product) "
                   "instead of one frame of STATE.")
@click.option("--before", "-b", default=20, type=int,
              help="Steps before the transition to write, at most to the start of the visit of FROM_STATE.")
@click.option("--after", "-a", default=20, type=int,
              help="Steps after the transition to write, at most to the end of the visit of STATE.")
@click.option("--max_runs", "-n", default=5, type=int, help="Number of runs to take frames from.")
@click.option("--largest_fragment", is_flag=True,
              help="Only write the atoms of the largest fragment (for a transition, of its first or last frame).")
@click.option("--output_dir", "-o", default=None, type=str, help="Where to write, FOLDER/frames by default.")
def extractFrames(folder, state, from_state, before, after, max_runs, largest_fragment, output_dir):
    """Write real MD geometries from the .mol files as .xyz, as starting points for DFT.

    Without --from_state: for each run with STATE, the frame in the middle of its longest visit of STATE.
    With --from_state A: for the first A -> STATE transition of each run, the frames around it (a first guess for a
    NEB: the first frame is in A, the last in STATE). Needs the .mol files and their _log.txt files in FOLDER.
    """
    from breakdown.frame_extractor import FrameExtractor
    try:
        paths = FrameExtractor(folder, state, from_state, before, after, max_runs, largest_fragment,
                               output_dir).extract()
    except ValueError as error:
        raise click.ClickException(str(error))
    if not paths:
        what = f"a {from_state} -> {state} transition" if from_state else state
        raise click.ClickException(f"No run with {what} and its .mol file found in {folder}.")
    for path in paths:
        click.echo(f"Wrote {path}")


@main_group.command()
@click.argument("folder")
@click.option("--max_states", "-n", default=100, type=int,
              help="Number of states in the table, the most common by number of runs.")
def makeReport(folder, max_states):
    """Write FOLDER/report.html: a sortable table of the states with drawings, the transitions, final fragments,
    settings and plots, in one file. Needs makechartsandstats to have been run."""
    from breakdown.html_report import writeReport
    try:
        click.echo(f"Wrote {writeReport(folder, maxStates=max_states)}")
    except ValueError as error:
        raise click.ClickException(str(error))


if __name__ == "__main__":
    main_group()