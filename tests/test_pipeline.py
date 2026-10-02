import json
import os
import numpy as np
import pytest
from breakdown.auto_detect_isomerization import AutoDetectIsomerization
from breakdown.generate_charts_stats import GenerateChartsStats
from breakdown.isomer_dicts import IsomerDicts
from breakdown.auto_detect_utils import processSmileToFragment, extractToDict
from conftest import ATOMS, writeMol


def makeRun(tmp_path, name, leaves):
    """Write a .mol file where the H on atom 3 leaves for good at step `leaves`."""
    from conftest import POSITIONS
    coords = np.tile(POSITIONS, (30, 1, 1))
    if leaves is not None:
        coords[leaves:, 3] = [12.0, 0.0, 0.0]
    path = os.path.join(tmp_path, name)
    writeMol(path, ATOMS, coords)
    return name


def test_analyse_file_infers_size_and_elements(tmp_path):
    name = makeRun(tmp_path, 'run_1.mol', 10)
    analyser = AutoDetectIsomerization(str(tmp_path), name)
    logPath = analyser.autoAnalyseFile()
    text = open(logPath).read()
    assert text.startswith('In file: run_1.mol')
    assert analyser.moleculeSize == 5
    assert analyser.moleculeElements == ['C', 'H', 'S']
    assert 'At step 10: Fragmentation' in text
    assert text.endswith('Final molecules present: C2HS, H')
    assert '[HH]' not in text


def test_dotted_file_names_and_wrong_size(tmp_path):
    name = makeRun(tmp_path, '2.5-thing.mol', None)
    logPath = AutoDetectIsomerization(str(tmp_path), name).autoAnalyseFile()
    assert os.path.basename(logPath) == '2.5-thing_log.txt'
    with pytest.raises(ValueError):
        AutoDetectIsomerization(str(tmp_path), name, moleculeSize=4).autoAnalyseFile()


def test_charts_and_stats_ignore_full_log(tmp_path):
    for i, leaves in enumerate([10, None, 20]):
        AutoDetectIsomerization(str(tmp_path), makeRun(tmp_path, f'run_{i}.mol', leaves)).autoAnalyseFile()
    args = (str(tmp_path), 'toy', 'TOY', 2, ['C', 'H', 'S'], 'thiol', False, False)
    GenerateChartsStats(*args).autoAnalyseInFolder()
    def withoutDate(path):
        return [line for line in open(path).read().splitlines() if '"created"' not in line]
    firstFull = open(tmp_path / 'full_log.txt').read()
    # A second run must not read full_log.txt as if it were a run
    GenerateChartsStats(*args).autoAnalyseInFolder()
    assert withoutDate(tmp_path / 'full_log.txt') == [line for line in firstFull.splitlines() if '"created"' not in line]
    assert firstFull.count('In file:') == 3
    assert (tmp_path / 'plots' / 'toy_dissociation.pdf').exists()

    names = extractToDict(firstFull, 'Automatically generated molecule names.',
                          'Total number of steps at which an isomer is present.')
    assert len(names) == 2
    runs = extractToDict(firstFull, 'Total number of runs in which an isomer appeared.',
                         'Number of runs at which each chemical fragment was present at the end.')
    assert sorted(runs.values()) == [2, 3]


def test_isomer_dicts_from_log_text():
    log = ('In file: a.mol \n'
           'At step 0:               [H]CC\n'
           'At step 1: Fragmentation [H].CC\n'
           'Number of 6-carbon rings: 1. Thiol group count: 1. \n'
           'At step 2:               [H].CC\n'
           'Final molecules present: C2, H\n\n'
           'In file: b.mol \n'
           'At step 0:               [H]CC\n'
           'Final molecules present: C2H\n')
    dicts = IsomerDicts(log, 'M', '.', 'thioketone', heteroatom='S')
    dicts.createIsomerDicts()
    assert dicts.isomerStats == {'[H]CC': 2, '[H].CC': 2}
    assert dicts.isomerRuns == {'[H]CC': 2, '[H].CC': 1}
    assert dicts.isomerizationStats == {'[H].CC': {'[H]CC': 1}}
    # The thioketone is gone after the fragmentation in the first run
    assert dicts.groupSteps == {0: 2, 1: 0, 2: 0}
    # The first run ended with a thiol, the second never changed so it kept its group
    assert dicts.endingGroups == {'thiol': 1, 'thioketone': 1}


def test_formula_ignores_bracket_hydrogens():
    assert processSmileToFragment('[H]C1[SH]C1') == 'C2HS'


def test_settings_are_recorded_and_checked(tmp_path):
    import json
    import warnings
    from breakdown.run_settings import parseSettings
    for i in range(2):
        AutoDetectIsomerization(str(tmp_path), makeRun(tmp_path, f'run_{i}.mol', 10),
                                breakDistance=2.5 + 0.5 * i).autoAnalyseFile()
    text = open(tmp_path / 'run_0_log.txt').read()
    assert text.splitlines()[1].startswith('Settings: ')
    settings = parseSettings(text)
    assert settings['break_distance'] == 2.5 and settings['formation_distance'] == 1.63
    assert settings['elements'] == ['C', 'H', 'S'] and settings['min_lifetime'] == 1
    args = (str(tmp_path), 'toy', 'TOY', 2, ['C', 'H', 'S'], None, False, False)
    with pytest.warns(UserWarning, match='different sets of settings'):
        GenerateChartsStats(*args).autoAnalyseInFolder()
    summary = json.load(open(tmp_path / 'run_summary.json'))
    assert summary['settings']['analysisConsistent'] is False
    assert len(summary['settings']['analysis']) == 2
    assert summary['settings']['statistics']['time_per_step_fs'] == 5.0
    assert 'Settings used for these results.' in open(tmp_path / 'full_log.txt').read()
