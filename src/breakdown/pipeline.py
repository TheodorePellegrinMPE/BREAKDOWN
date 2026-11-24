"""The whole analysis of a folder in one command, with its options in a TOML file.

`breakdown run FOLDER` analyses the .mol files, makes full_log.txt with the statistics and plots, the flow charts and
the HTML report. The options come from FOLDER/breakdown.toml (or the file given with --config), anything not in the
file has the same default as the separate commands. The options that were used, with the names and start state that
were picked filled in, are written to FOLDER/breakdown_used.toml, so
`breakdown run FOLDER --config FOLDER/breakdown_used.toml` repeats the analysis.

A run log that already exists is kept when it is complete and was made with the same bond settings, otherwise the
.mol file is analysed again, so changing a distance in the file never mixes old and new logs.
"""
import copy
import fnmatch
import os
import shutil
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib

from breakdown.auto_detect_isomerization import AutoDetectIsomerization
from breakdown.run_log import LOG_SUFFIX, isComplete
from breakdown.run_settings import parseSettings
from breakdown.run_summary import loadBlock

CONFIG_FILE = 'breakdown.toml'
USED_CONFIG_FILE = 'breakdown_used.toml'

# Every option, with its default and what it means. An empty string, 0 or an empty list means "not set": read from
# the files, or not used. TOML has no null.
CONFIG_OPTIONS = {
    'analysis': {
        'pattern': ('*.mol', 'Pattern for the names of the .mol files to analyse (files with "keep" in the name are skipped).'),
        'jobs': (0, 'Number of files analysed at the same time, 0 for one per processor.'),
        'reanalyse': (False, 'Analyse every .mol file again, even when its log was made with the same settings.'),
        'formation_distance': (1.63, 'Distance in angstroms under which two atoms become bonded.'),
        'break_distance': (2.5, 'Distance in angstroms over which a bond breaks. The most important setting.'),
        'min_lifetime': (1, 'Number of steps a new state must last to be reported, 1 reports everything.'),
        'molecule_size': (0, 'Number of atoms, 0 to read it from the files.'),
        'elements': ([], 'Elements of the molecule, like ["C", "H", "O"], empty to read them from the files.'),
    },
    'statistics': {
        'chart_name': ('', 'Name used in the plot file names, the folder name if empty.'),
        'molecule_name': ('', 'Start of the state names (PHE gives PHE_A1, ...), the folder name if empty.'),
        'func_group': ('', 'Functional group to track, like "hydroxyl" or "ketone_2", empty for none.'),
        'time_per_step_fs': (5.0, 'Simulated time between output steps in fs (0.5 fs time step, output every 10).'),
        'legend_carbons': ([], 'Carbon counts to show in the dissociation legend, empty for all.'),
        'draw': (False, 'Draw every state (needed for the drawings in the flow charts, slow with thousands of states).'),
        'coordinates': (False, 'Make force-field .xyz files for every state.'),
    },
    'flowcharts': {
        'enabled': (True, 'Make the flow charts (skipped with a message when Graphviz "dot" is not installed).'),
        'start_state': ('', 'Name of the start state, empty for the state most runs start in.'),
        'num_top_states': (10, 'Number of most common states to draw a flow chart to.'),
        'num_scattered': (5, 'Number of most common fragmented states to draw a flow chart to.'),
        'max_depth_best': (5, 'Depth to which the best path is looked for.'),
        'max_depth_possible': (20, 'Depth to which any path is looked for.'),
    },
    'report': {
        'enabled': (True, 'Write report.html.'),
        'max_states': (100, 'Number of states in the report table, the most common by number of runs.'),
    },
}


class ConfigError(ValueError):
    """A configuration file that cannot be used."""


def defaultConfig():
    """All options with their defaults.

    Returns:
        dict: Section to {option: value}.
    """
    return {section: {name: copy.deepcopy(option[0]) for name, option in options.items()}
            for section, options in CONFIG_OPTIONS.items()}


def loadConfig(path=None):
    """Read a configuration file and fill in the defaults.

    Args:
        path (str, optional): TOML file, None for only the defaults.

    Returns:
        dict: Section to {option: value}.

    Raises:
        ConfigError: Raised if the file has an unknown section or option, or a value of the wrong type.
    """
    config = defaultConfig()
    if path is None:
        return config
    try:
        with open(path, 'rb') as file:
            given = tomllib.load(file)
    except tomllib.TOMLDecodeError as error:
        raise ConfigError(f'{path} is not valid TOML: {error}')
    for section, options in given.items():
        if section not in CONFIG_OPTIONS:
            raise ConfigError(f'Unknown section [{section}] in {path}, the sections are {", ".join(CONFIG_OPTIONS)}.')
        if not isinstance(options, dict):
            raise ConfigError(f'{section} in {path} must be a [section].')
        for name, value in options.items():
            if name not in CONFIG_OPTIONS[section]:
                raise ConfigError(f'Unknown option {name} in [{section}] of {path}, the options are '
                                  f'{", ".join(CONFIG_OPTIONS[section])}.')
            default = CONFIG_OPTIONS[section][name][0]
            if isinstance(default, float) and isinstance(value, int) and not isinstance(value, bool):
                value = float(value)
            if type(value) is not type(default):
                raise ConfigError(f'{name} in [{section}] of {path} must be a {type(default).__name__}, '
                                  f'not {value!r}.')
            config[section][name] = value
    return config


def formatTomlValue(value):
    """Write a value as TOML (str, bool, int, float or a list of those)."""
    if isinstance(value, bool):
        return 'true' if value else 'false'
    if isinstance(value, str):
        return '"' + value.replace('\\', '\\\\').replace('"', '\\"') + '"'
    if isinstance(value, list):
        return '[' + ', '.join(formatTomlValue(item) for item in value) + ']'
    return repr(value)


def formatConfig(config, comments=True):
    """Write a configuration as TOML text.

    Args:
        config (dict): Section to {option: value}.
        comments (bool): Put the meaning of every option above it.

    Returns:
        str: The TOML text.
    """
    lines = []
    for section, options in config.items():
        lines.append(f'[{section}]')
        for name, value in options.items():
            if comments:
                lines.append(f'# {CONFIG_OPTIONS[section][name][1]}')
            lines.append(f'{name} = {formatTomlValue(value)}')
        lines.append('')
    return '\n'.join(lines)


def analyseOneFile(arguments):
    """Analyse one file, for use by worker processes, returning the name of the file and an error message."""
    folder, fileName, options = arguments
    try:
        AutoDetectIsomerization(folder, fileName, *options).autoAnalyseFile()
        return fileName, None
    except Exception as error:
        return fileName, f'{type(error).__name__}: {error}'


def listMolFiles(folder, pattern):
    """The .mol files of a folder that match a pattern, without the ones with keep in their name, sorted."""
    return sorted(name for name in os.listdir(folder) if fnmatch.fnmatch(name, pattern) and 'keep' not in name)


def logPathFor(folder, molName):
    return os.path.join(folder, f'{os.path.splitext(molName)[0]}{LOG_SUFFIX}')


def analyseFiles(folder, fileNames, options, jobs):
    """Analyse .mol files several at a time, printing the progress.

    Args:
        folder (str): Folder with the files.
        fileNames (list): Names of the .mol files.
        options (tuple): molecule size, elements, formation distance, break distance, min lifetime.
        jobs (int): Number of files analysed at the same time.

    Returns:
        list: Names of the files that failed.
    """
    failures = []
    with ProcessPoolExecutor(max_workers=max(1, min(jobs, len(fileNames)))) as executor:
        tasks = [(folder, name, options) for name in fileNames]
        for done, (name, error) in enumerate(executor.map(analyseOneFile, tasks), start=1):
            if error:
                failures.append(name)
                print(f'[{done}/{len(fileNames)}] {name} FAILED: {error}')
            else:
                print(f'[{done}/{len(fileNames)}] {name}')
    return failures


def logMatchesSettings(logPath, analysis):
    """True if a run log is complete and was made with the bond settings of the configuration.

    Args:
        logPath (str): The run log.
        analysis (dict): The [analysis] section.

    Returns:
        bool: Whether the log can be kept.
    """
    if not os.path.exists(logPath):
        return False
    with open(logPath) as file:
        text = file.read()
    settings = parseSettings(text)
    if not isComplete(text) or settings is None:
        return False
    for name in ('formation_distance', 'break_distance', 'min_lifetime'):
        if settings.get(name) != analysis[name]:
            return False
    if analysis['molecule_size'] and settings.get('molecule_size') != analysis['molecule_size']:
        return False
    if analysis['elements'] and settings.get('elements') != analysis['elements']:
        return False
    return True


def findStartState(folder):
    """The state most runs start in, from the runSequences block.

    Args:
        folder (str): Folder with full_log.txt.

    Returns:
        str or None: Its name, None if there is no runSequences block.
    """
    sequences = loadBlock(folder, 'runSequences')
    starts = Counter(sequence['0'] for sequence in (sequences or {}).values() if '0' in sequence)
    return starts.most_common(1)[0][0] if starts else None


def runPipeline(folder, config):
    """Analyse a folder from .mol files to the report.

    Args:
        folder (str): Folder with the .mol files.
        config (dict): The configuration, see loadConfig.

    Returns:
        dict: The configuration that was used, with the chart and molecule names and the start state filled in.

    Raises:
        ConfigError: Raised if there is nothing to analyse or files failed to analyse.
    """
    # Imported here so that reading a configuration does not load matplotlib and graphviz
    from breakdown.generate_charts_stats import GenerateChartsStats, inferMoleculeFromLogs
    from breakdown.html_report import writeReport
    from breakdown.molecule_transition_graph import MoleculeTransitionGraph

    config = copy.deepcopy(config)
    analysis, statistics = config['analysis'], config['statistics']
    folderName = os.path.basename(os.path.normpath(os.path.abspath(folder)))

    print('== Analysing trajectories')
    molFiles = listMolFiles(folder, analysis['pattern'])
    if not molFiles:
        raise ConfigError(f'No files matching {analysis["pattern"]} in {folder}.')
    toAnalyse = [name for name in molFiles
                 if analysis['reanalyse'] or not logMatchesSettings(logPathFor(folder, name), analysis)]
    print(f'{len(molFiles) - len(toAnalyse)} of {len(molFiles)} files already have a log made with these settings.')
    if toAnalyse:
        options = (analysis['molecule_size'] or None, list(analysis['elements']), analysis['formation_distance'],
                   analysis['break_distance'], analysis['min_lifetime'])
        failures = analyseFiles(folder, toAnalyse, options, analysis['jobs'] or os.cpu_count() or 1)
        if failures:
            raise ConfigError(f'{len(failures)} files failed to analyse: {", ".join(failures)}. Fix or move them, '
                              'the statistics would be missing these runs.')

    print('== Statistics and plots')
    numCarbons, elements = inferMoleculeFromLogs(folder)
    statistics['chart_name'] = statistics['chart_name'] or folderName
    statistics['molecule_name'] = statistics['molecule_name'] or folderName
    if analysis['elements']:
        elements = list(analysis['elements'])
    GenerateChartsStats(folder, statistics['chart_name'], statistics['molecule_name'], numCarbons, elements,
                        statistics['func_group'] or None, statistics['draw'], statistics['coordinates'],
                        legendCarbons=statistics['legend_carbons'] or None,
                        timePerStepFs=statistics['time_per_step_fs']).autoAnalyseInFolder()
    print(f'Wrote {os.path.join(folder, "full_log.txt")} and the plots.')

    flowcharts = config['flowcharts']
    if flowcharts['enabled']:
        print('== Flow charts')
        startState = flowcharts['start_state'] or findStartState(folder)
        flowcharts['start_state'] = startState or ''
        if shutil.which('dot') is None:
            print('Skipped: Graphviz "dot" is not installed.')
        elif startState is None:
            print('Skipped: no start state found.')
        else:
            sections = [loadBlock(folder, key) for key in
                        ('isomerNames', 'isomerStats', 'isomerizationStats', 'netIsomerizationStats')]
            try:
                MoleculeTransitionGraph(folder, *sections[:3], sections[3], startState,
                                        flowcharts['num_top_states'], flowcharts['num_scattered'],
                                        flowcharts['max_depth_best'], flowcharts['max_depth_possible'],
                                        []).makeTransitionGraph()
                print(f'Wrote the flow charts from {startState} to {os.path.join(folder, "transitionGraphs")}.')
            except Exception as error:
                print(f'Flow charts failed ({type(error).__name__}: {error}), the rest of the analysis is kept.')

    if config['report']['enabled']:
        print('== Report')
        print(f'Wrote {writeReport(folder, maxStates=config["report"]["max_states"])}.')

    usedPath = os.path.join(folder, USED_CONFIG_FILE)
    with open(usedPath, 'w') as file:
        file.write('# Options used by the last "breakdown run" of this folder.\n'
                   f'# Repeat it with: breakdown run {folder} --config {usedPath}\n\n' + formatConfig(config))
    print(f'Wrote {usedPath}.')
    return config
