import json
import pytest
from breakdown.run_log import (BLANK_EVENT, Step, formatRunLog, formatStepLine, isComplete, iterSteps, parseRun,
                               parseRuns)
from breakdown.run_summary import (BLOCKS, STATISTICS_START, formatBlocks, getHeadings, loadBlock, saveSummary)

SETTINGS = {'break_distance': 2.5, 'formation_distance': 1.63}


def makeLog(name='a.mol'):
    steps = '\n'.join([formatStepLine(0, BLANK_EVENT, 'CC'),
                       formatStepLine(1, 'Fragmentation', 'C.C'),
                       'Number of 3-carbon rings: 1. ',
                       formatStepLine(2, BLANK_EVENT, 'C.C')]) + '\n'
    return formatRunLog(name, SETTINGS, steps, ['C', 'C'])


def test_run_log_round_trip():
    text = makeLog()
    assert text.splitlines()[0] == 'In file: a.mol '
    assert text.splitlines()[2] == 'At step 0:               CC'
    assert isComplete(text)
    run = parseRun(text)
    assert run.name == 'a.mol' and run.settings == SETTINGS and run.final == ['C', 'C']
    assert run.steps == [Step(0, '', 'CC', ''), Step(1, 'Fragmentation', 'C.C', 'Number of 3-carbon rings: 1. '),
                         Step(2, '', 'C.C', '')]
    assert (run.isomerized, run.everFragmented, run.fragmentedAtEnd) == (False, True, True)


def test_logs_without_header_get_numbered_names():
    text = 'At step 0:               CC\nFinal molecules present: C2\n\nAt step 0:               CC\nFinal molecules present: C2\n'
    assert [run.name for run in parseRuns(text)] == ['run_1', 'run_2']
    assert not isComplete('In file: a.mol \nAt step 0:               CC\n')


def test_streamed_steps_stop_at_the_statistics():
    lines = (makeLog() + '\n\n' + STATISTICS_START + '\n{}\nAt step 9:               CC\n').splitlines(True)
    assert [step for step, _, _ in iterSteps(lines, stopAt=STATISTICS_START)] == [0, 1, 2]


def test_block_headings_never_change():
    # Files made by older versions are found by these headings
    assert list(BLOCKS.values()) == [
        'Automatically generated molecule names.',
        'Total number of steps at which an isomer is present.',
        'Total number of runs in which an isomer appeared.',
        'Number of runs at which each chemical fragment was present at the end.',
        'Total number of transitions to an isomer from an isomer.',
        'Net number of transitions to an isomer from an isomer.',
        'Number of transitions only including fragmentation.',
        'Chemical names of fragments of this isomerization state.',
        'Molecular masses of different observed isomers.',
        'Runs containing original functional group {group} at step.',
        'Runs with each {heteroatom} in each final group and/or fragment.',
        'First step of each event type in each run.',
        'Connectivity key of each molecule state.',
        'States visited in each run, by the step they were entered.',
        'Settings used for these results.',
    ]
    assert getHeadings('endingGroups') == ['Runs with each oxygen in each final group and/or fragment.',
                                           'Runs with each sulfur in each final group and/or fragment.']
    with pytest.raises(ValueError):
        getHeadings('groupSteps')


@pytest.fixture(params=['summary', 'full_log'])
def writtenFolder(tmp_path, request):
    blocks = {'isomerNames': {'CC': 'M_A1'}, 'groupSteps': {0: 2, 1: 1}, 'endingGroups': {'thiol': 2},
              'firstEvents': {'a.mol': {'fragmentation': 1}}}
    text, sections = formatBlocks(blocks, group='thiol', heteroatom='S')
    with open(tmp_path / 'full_log.txt', 'w') as file:
        file.write(makeLog() + '\n\n' + text)
    if request.param == 'summary':
        saveSummary(str(tmp_path), {'sections': sections})
    return str(tmp_path)


def test_blocks_are_read_back_by_key(writtenFolder):
    assert loadBlock(writtenFolder, 'isomerNames') == {'CC': 'M_A1'}
    assert loadBlock(writtenFolder, 'groupSteps', group='thiol') == {'0': 2, '1': 1}
    assert loadBlock(writtenFolder, 'groupSteps', group='hydroxyl') is None
    assert loadBlock(writtenFolder, 'endingGroups') == {'thiol': 2}
    assert loadBlock(writtenFolder, 'firstEvents') == {'a.mol': {'fragmentation': 1}}
    assert loadBlock(writtenFolder, 'isomerStats') is None


def test_blocks_are_written_in_catalogue_order():
    text, _ = formatBlocks({'settings': {}, 'isomerNames': {}})
    assert text.splitlines()[0] == BLOCKS['isomerNames'] and json.loads(text.splitlines()[1]) == {}
    assert text.splitlines()[2] == BLOCKS['settings']
