"""The run log: the text file `<name>_log.txt` that analysedemon writes for one trajectory.

This module is the only place that knows the format, for writing and for reading it::

    In file: <name>.mol
    Settings: {"break_distance": 2.5, ...}
    At step 0:               <smiles>
    At step 1: Isomerization <smiles>
    <analysis of the new state>
    ...
    Final molecules present: <formula>, <formula>

full_log.txt starts with the run logs of a folder written one after the other, so the same functions read it.
Everything downstream works on `Run` objects instead of parsing lines itself.
"""
import os
from dataclasses import dataclass, field
from typing import List, NamedTuple, Optional
from breakdown.run_settings import SETTINGS_PREFIX, formatSettingsLine, parseSettings

EVENT_LABELS = ('Isomerization', 'Fragmentation', 'Recombination')
# Written instead of a label when nothing changed, the same width as the labels
BLANK_EVENT = ' ' * 13
STEP_PREFIX = 'At step'
HEADER_PREFIX = 'In file'
FINAL_PREFIX = 'Final molecules present: '
LOG_SUFFIX = '_log.txt'
FULL_LOG_FILE = 'full_log.txt'


class Step(NamedTuple):
    """One step of a run: its number, the event label ('' if nothing changed), the SMILES string and the
    analysis line written after an event ('' if there is none)."""
    step: int
    event: str
    smiles: str
    analysis: str


@dataclass
class Run:
    """One run read from a run log.

    Attributes:
        name (str): The .mol file name of the run, or run_<n> for logs without an In file line.
        steps (list): The Steps of the run.
        final (list or None): Chemical formulas present at the end, None if the line is missing.
        settings (dict or None): Settings the run was analysed with, None for logs of older versions.
    """
    name: str
    steps: List[Step] = field(default_factory=list)
    final: Optional[List[str]] = None
    settings: Optional[dict] = None

    @property
    def isomerized(self):
        """True if the run had any isomerization."""
        return any(step.event == 'Isomerization' for step in self.steps)

    @property
    def everFragmented(self):
        """True if the run fragmented at some point, even if it recombined later."""
        return any(step.event == 'Fragmentation' for step in self.steps)

    @property
    def fragmentedAtEnd(self):
        """True if the last step is in more than one piece. A run that loses a hydrogen for a moment and takes
        it back (fragmentation followed by recombination) is not fragmented."""
        return bool(self.steps) and '.' in self.steps[-1].smiles


def formatStepLine(step, event, smiles):
    """The line of one step.

    Args:
        step (int): Step number.
        event (str): One of EVENT_LABELS, or BLANK_EVENT.
        smiles (str): SMILES string of the state at this step.

    Returns:
        str: The line, without a newline.
    """
    return f'{STEP_PREFIX} {step}: {event} {smiles}'


def formatRunLog(fileName, settings, stepsText, finalFormulas):
    """The whole text of a run log.

    Args:
        fileName (str): Name of the .mol file.
        settings (dict): Settings of the analysis, written as the Settings line.
        stepsText (str): Step lines and analysis lines, ending with a newline.
        finalFormulas (list): Chemical formulas present at the last step, largest first.

    Returns:
        str: Text of the log file.
    """
    return (f'{HEADER_PREFIX}: {fileName} \n' + formatSettingsLine(settings) + '\n' + stepsText +
            FINAL_PREFIX + ', '.join(finalFormulas))


def parseStepLine(line):
    """Read an At step line.

    Args:
        line (str): A line that starts with 'At step'.

    Returns:
        Tuple[int, str, str]: Step number, event label ('' if there was none) and SMILES string.
    """
    parts = line.split()
    event = parts[3] if len(parts) == 5 and parts[3] in EVENT_LABELS else ''
    return int(parts[2][:-1]), event, parts[-1]


def parseRuns(text):
    """Read every run in the text of one or more run logs (or of full_log.txt).

    Args:
        text (str): Text of run logs.

    Returns:
        list: The Runs, in the order they appear.
    """
    runs = []
    run = None
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if line.startswith(HEADER_PREFIX):
            run = Run(line.split(': ', 1)[1].strip())
            runs.append(run)
        elif line.startswith(STEP_PREFIX):
            if run is None or run.final is not None:
                run = Run(f'run_{len(runs) + 1}')
                runs.append(run)
            step, event, smiles = parseStepLine(line)
            analysis = lines[index + 1] if event and index + 1 < len(lines) else ''
            run.steps.append(Step(step, event, smiles, analysis))
        elif line.startswith(FINAL_PREFIX) and run is not None:
            run.final = line.split(': ', 1)[1].split(', ')
        elif line.startswith(SETTINGS_PREFIX) and run is not None and not run.steps:
            run.settings = parseSettings(line)
    return runs


def parseRun(text):
    """Read the text of one run log.

    Args:
        text (str): Text of the log.

    Returns:
        Run: The run, with no steps if the text has none.
    """
    runs = parseRuns(text)
    return runs[0] if runs else Run('run_1')


def isComplete(text):
    """True if the text of a run log has steps and its final line, logs of unfinished analyses do not."""
    return FINAL_PREFIX in text and STEP_PREFIX in text


def listRunLogs(folder):
    """Names of the run logs in a folder, sorted, without full_log.txt.

    Args:
        folder (str): Folder with _log.txt files.

    Returns:
        list: File names.
    """
    return sorted(name for name in os.listdir(folder)
                  if name.endswith(LOG_SUFFIX) and not name.endswith(FULL_LOG_FILE))


def iterSteps(lines, stopAt=None):
    """Read the steps from lines one at a time, for files too big to load, such as full_log.txt.

    Args:
        lines (iterable): Lines, like an open file.
        stopAt (str, optional): Stop at the first line that starts with this.

    Yields:
        Tuple[int, str, str]: Step number, event label ('' if there was none) and SMILES string.
    """
    for line in lines:
        if stopAt is not None and line.startswith(stopAt):
            return
        if line.startswith(STEP_PREFIX):
            yield parseStepLine(line)
