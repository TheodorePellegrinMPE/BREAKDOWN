"""The statistics blocks of a folder of runs, and the small json summary that makechartsandstats writes next to
full_log.txt.

full_log.txt is the run logs followed by statistics blocks, each a heading line and a json dictionary. The same
blocks are kept in run_summary.json under `sections`. This module is the only place that knows the headings:
the code that makes the statistics and the code that reads them use the keys of BLOCKS.

full_log.txt is often 100 MB. The plotting commands that compare molecules only need the statistics blocks
at its end and a few per-step counts, so they read this summary instead when it exists, and fall back to
streaming through full_log.txt when it does not (for logs made by older versions).
"""
import json
import math
import os
from breakdown.run_log import FULL_LOG_FILE

SUMMARY_FILE = 'run_summary.json'
SUMMARY_VERSION = 1

# Heading of every statistics block, in the order they are written. Files made by older versions are found by
# these headings, so never change their text. {group} is the tracked functional group, {heteroatom} oxygen or sulfur.
BLOCKS = {
    'isomerNames': 'Automatically generated molecule names.',
    'isomerStats': 'Total number of steps at which an isomer is present.',
    'isomerRuns': 'Total number of runs in which an isomer appeared.',
    'endFragmentStats': 'Number of runs at which each chemical fragment was present at the end.',
    'isomerizationStats': 'Total number of transitions to an isomer from an isomer.',
    'netIsomerizationStats': 'Net number of transitions to an isomer from an isomer.',
    'fragmentationStats': 'Number of transitions only including fragmentation.',
    'isomerChemicalNames': 'Chemical names of fragments of this isomerization state.',
    'isomerMasses': 'Molecular masses of different observed isomers.',
    'groupSteps': 'Runs containing original functional group {group} at step.',
    'endingGroups': 'Runs with each {heteroatom} in each final group and/or fragment.',
    'firstEvents': 'First step of each event type in each run.',
    'stateKeys': 'Connectivity key of each molecule state.',
    'runSequences': 'States visited in each run, by the step they were entered.',
    'settings': 'Settings used for these results.',
}
HETEROATOM_NAMES = {'O': 'oxygen', 'S': 'sulfur'}
# The first line after the run logs in full_log.txt
STATISTICS_START = BLOCKS['isomerNames']
FIRST_EVENTS_HEADING = BLOCKS['firstEvents']


def getHeadings(key, group=None, heteroatom=None):
    """The heading lines a block can have.

    Args:
        key (str): Key of the block in BLOCKS.
        group (str, optional): Functional group, needed for groupSteps.
        heteroatom (str, optional): O or S for endingGroups, both are tried if not given.

    Returns:
        list: Headings, the one that is written first.

    Raises:
        ValueError: Raised if the key is unknown or the block needs a group that was not given.
    """
    if key not in BLOCKS:
        raise ValueError(f'Unknown statistics block {key}.')
    if '{group}' in BLOCKS[key] and group is None:
        raise ValueError(f'The {key} block needs a functional group.')
    elements = [heteroatom] if heteroatom else list(HETEROATOM_NAMES)
    return list(dict.fromkeys(BLOCKS[key].format(group=group, heteroatom=HETEROATOM_NAMES[element])
                              for element in elements))


def formatBlocks(blocks, group=None, heteroatom='O'):
    """Write statistics as the blocks at the end of full_log.txt.

    Args:
        blocks (dict): Key of BLOCKS to its dictionary, blocks that are not given are left out.
        group (str, optional): Functional group, for the heading of groupSteps.
        heteroatom (str): O or S, for the heading of endingGroups.

    Returns:
        Tuple[str, dict]: The text, and the sections of the summary (heading to the dictionary as read back from json).
    """
    text = ''
    sections = {}
    for key in BLOCKS:
        if key in blocks:
            heading = getHeadings(key, group, heteroatom)[0]
            blockText = json.dumps(blocks[key], indent=4) + '\n'
            text += heading + '\n' + blockText
            sections[heading] = json.loads(blockText)
    return text, sections


def getSummaryPath(folder):
    return os.path.join(folder, SUMMARY_FILE)


def saveSummary(folder, summary):
    """Save the summary of a folder of runs.

    Args:
        folder (str): Folder with full_log.txt.
        summary (dict): Dictionary with the keys version, molecule, heteroatom, timePerStepFs,
            numRuns, sections (heading of a full_log.txt block to its dictionary),
            firstEvents (see IsomerDicts.createFirstEventStats) and elementalIntact.
    """
    summary = dict(summary, version=SUMMARY_VERSION)
    with open(getSummaryPath(folder), 'w') as file:
        json.dump(summary, file)


def loadSummary(folder):
    """Load the summary of a folder, if it exists and is not older than full_log.txt.

    Args:
        folder (str): Folder with full_log.txt.

    Returns:
        dict or None: The summary, None if there is none or it is out of date.
    """
    path = getSummaryPath(folder)
    fullLogPath = os.path.join(folder, FULL_LOG_FILE)
    if not os.path.exists(path):
        return None
    if os.path.exists(fullLogPath) and os.path.getmtime(path) < os.path.getmtime(fullLogPath):
        return None
    with open(path) as file:
        summary = json.load(file)
    return summary if summary.get('version') == SUMMARY_VERSION else None


def readSectionFromFullLog(folder, heading):
    """Read one statistics block from full_log.txt without loading the whole file.

    Args:
        folder (str): Folder with full_log.txt.
        heading (str): Heading line of the block, without the newline.

    Returns:
        dict or None: The block, None if the file or the block does not exist.
    """
    path = os.path.join(folder, FULL_LOG_FILE)
    if not os.path.exists(path):
        return None
    with open(path, encoding='utf-8') as file:
        for line in file:
            if line.rstrip('\n') == heading:
                lines = []
                for blockLine in file:
                    lines.append(blockLine)
                    if blockLine.rstrip('\n') in ('{}', '}'):
                        break
                return json.loads(''.join(lines))
    return None


def loadSection(folder, *headings):
    """Get a statistics block of a folder by its heading, from the summary if possible.

    Args:
        folder (str): Folder with full_log.txt.
        *headings (str): One or more headings the block may have, the first one that exists is used.

    Returns:
        dict or None: The block, None if it is not there.
    """
    summary = loadSummary(folder)
    for heading in headings:
        if summary is not None and heading in summary['sections']:
            return summary['sections'][heading]
        if summary is None:
            section = readSectionFromFullLog(folder, heading)
            if section is not None:
                return section
    return None


def loadBlock(folder, key, group=None, heteroatom=None):
    """Get a statistics block of a folder, from the summary if possible.

    Args:
        folder (str): Folder with full_log.txt.
        key (str): Key of the block in BLOCKS.
        group (str, optional): Functional group, needed for groupSteps.
        heteroatom (str, optional): O or S for endingGroups, both are tried if not given.

    Returns:
        dict or None: The block, None if it is not there.
    """
    return loadSection(folder, *getHeadings(key, group, heteroatom))


def kaplanMeier(times, observed, confidenceZ=1.96):
    """Kaplan-Meier estimate of the fraction of runs that have not had an event yet.

    Runs that ended without the event are censored, they stay in the count up to their end.
    The confidence interval uses Greenwood's variance.

    Args:
        times (list): Time of the event, or of the end of the run if there was no event, for every run.
        observed (list): True for every run where the event happened.
        confidenceZ (float): Number of standard errors of the interval, 1.96 is 95%.

    Returns:
        dict: Lists times, survival, lower and upper, each starting at time 0 with survival 1,
        and median, the first time survival is 0.5 or lower (None if it never gets there).
    """
    order = sorted(range(len(times)), key=lambda i: times[i])
    atRisk = len(times)
    survival = 1.0
    greenwood = 0.0
    result = {'times': [0.0], 'survival': [1.0], 'lower': [1.0], 'upper': [1.0], 'median': None}
    index = 0
    while index < len(order):
        time = times[order[index]]
        events = 0
        leaving = 0
        while index < len(order) and times[order[index]] == time:
            events += observed[order[index]]
            leaving += 1
            index += 1
        if events:
            survival *= 1 - events / atRisk
            if atRisk > events:
                greenwood += events / (atRisk * (atRisk - events))
            error = confidenceZ * survival * math.sqrt(greenwood) if survival > 0 else 0.0
            result['times'].append(time)
            result['survival'].append(survival)
            result['lower'].append(max(0.0, survival - error))
            result['upper'].append(min(1.0, survival + error))
            if result['median'] is None and survival <= 0.5:
                result['median'] = time
        atRisk -= leaving
    return result


def getSurvivalCurve(firstEvents, event, timePerStepFs):
    """Kaplan-Meier curve of the time to the first isomerization, fragmentation or any change.

    Args:
        firstEvents (dict): Per run, the dictionary made by IsomerDicts.createFirstEventStats.
        event (str): One of isomerization, fragmentation, recombination or anyChange.
        timePerStepFs (float): Simulated time between steps, in femtoseconds.

    Returns:
        dict: The result of kaplanMeier with times in picoseconds, plus numRuns and numEvents.
    """
    stepToPs = timePerStepFs / 1000
    times, observed = [], []
    for record in firstEvents.values():
        step = record[event]
        observed.append(step is not None)
        times.append((step if step is not None else record['lastStep']) * stepToPs)
    curve = kaplanMeier(times, observed)
    curve['numRuns'] = len(times)
    curve['numEvents'] = sum(observed)
    return curve


def getElementalIntactSeries(runs, elements, checkFunction):
    """Count, at every step, the runs in which each element is not in a fragment away from the parent.

    Args:
        runs (list): Runs as read by run_log.parseRuns.
        elements (list): Element symbols to track, apart from carbon which is always tracked.
        checkFunction (callable): Function(smiles, element) returning
            (carbon in fragment, element in fragment), see checkElementalDissociation.

    Returns:
        dict: Element symbol to {step: number of runs where it is intact}, and 'runsAtStep'.
    """
    series = {element: {} for element in ['C'] + list(elements)}
    runsAtStep = {}
    for run in runs:
        for step, _, smiles, _ in run.steps:
            runsAtStep[step] = runsAtStep.get(step, 0) + 1
            carbonLost = None
            for element in elements or ['O']:
                carbonLost, elementLost = checkFunction(smiles, element)
                if element in elements:
                    series[element][step] = series[element].get(step, 0) + int(not elementLost)
            series['C'][step] = series['C'].get(step, 0) + int(not carbonLost)
    series['runsAtStep'] = runsAtStep
    return series


def resolveTimePerStep(folders, timePerStepFs=None, default=5.0):
    """Pick the time between steps: the given one, else the one stored in the first summary, else the default.

    Args:
        folders (list): Folders with full_log.txt.
        timePerStepFs (float or None): Value given by the user.
        default (float): Femtoseconds, 0.5 fs time steps with output every 10 steps.

    Returns:
        float: Time between steps in femtoseconds.
    """
    if timePerStepFs:
        return timePerStepFs
    for folder in folders:
        summary = loadSummary(folder)
        if summary:
            return summary['timePerStepFs']
    return default
