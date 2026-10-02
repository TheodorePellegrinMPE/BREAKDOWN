"""Settings of an analysis, recorded in the output so that a result can be traced back to how it was made.

Every run log gets a line `Settings: {...}` under its `In file:` line with the distances and options used to find
the bonds. `makechartsandstats` collects them, warns if the runs of a folder were analysed with different settings, and
writes them with the version of the program into full_log.txt and run_summary.json.
"""
import json
import os
import platform
import subprocess
from datetime import datetime, timezone

SETTINGS_PREFIX = 'Settings: '
SETTINGS_HEADING = 'Settings used for these results.'


def getVersionInfo():
    """Version of the program, the git commit it was run from if there is one, and of the main libraries.

    Returns:
        dict: Version information, values are None if they could not be found.
    """
    try:
        from importlib.metadata import version
        package = version('breakdown')
    except Exception:
        package = None
    commit = None
    try:
        folder = os.path.dirname(os.path.abspath(__file__))
        commit = subprocess.run(['git', 'describe', '--always', '--dirty', '--tags'], cwd=folder, capture_output=True,
                                text=True, timeout=5).stdout.strip() or None
    except Exception:
        pass
    try:
        from rdkit import rdBase
        rdkit = rdBase.rdkitVersion
    except Exception:
        rdkit = None
    return {'breakdown': package, 'git': commit, 'rdkit': rdkit, 'python': platform.python_version()}


def formatSettingsLine(settings):
    """Turn the settings of one run into the line that is written under its In file line."""
    return SETTINGS_PREFIX + json.dumps(settings, sort_keys=True)


def parseSettings(logString):
    """Get the settings from the log of one run.

    Args:
        logString (str): Log of one run.

    Returns:
        dict or None: The settings, None for logs made by versions that did not record them.
    """
    for line in logString.split('\n', 5)[:5]:
        if line.startswith(SETTINGS_PREFIX):
            return json.loads(line[len(SETTINGS_PREFIX):])
    return None


def collectSettings(runSettings, statisticsOptions):
    """Combine the settings of all runs and of the statistics step.

    Args:
        runSettings (list): The result of parseSettings for every run (None where there were none).
        statisticsOptions (dict): The options that makechartsandstats was run with.

    Returns:
        dict: version information, the distinct settings the runs were analysed with, whether all runs used the same,
        the options of the statistics step and the date.
    """
    known = [s for s in runSettings if s is not None]
    distinct = []
    for settings in known:
        if settings not in distinct:
            distinct.append(settings)
    return {
        'program': getVersionInfo(),
        'analysis': distinct,
        'runsWithoutRecordedSettings': len(runSettings) - len(known),
        'analysisConsistent': len(distinct) <= 1 and len(known) == len(runSettings),
        'statistics': statisticsOptions,
        'created': datetime.now(timezone.utc).isoformat(timespec='seconds'),
    }
