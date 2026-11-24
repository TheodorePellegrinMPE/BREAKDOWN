"""Connectivity keys, run sequences, MD frame extraction, the HTML report and `breakdown run`."""
import json
import os
import numpy as np
import pytest
from breakdown.auto_detect_isomerization import AutoDetectIsomerization
from breakdown.frame_extractor import FrameExtractor, getVisits
from breakdown.generate_charts_stats import GenerateChartsStats
from breakdown.html_report import writeReport
from breakdown.pipeline import ConfigError, defaultConfig, formatConfig, loadConfig, runPipeline
from breakdown.run_log import parseRun
from breakdown.run_summary import loadBlock
from breakdown.state_keys import largestFragment, matchStates, stateKey
from conftest import ATOMS, POSITIONS, writeMol

SHIFTED_H = [1.5, 1.1, 0.0]   # H 3 next to carbon 1 instead of carbon 0
GONE_H = [12.0, 0.0, 0.0]


def writeShiftRun(folder, name, shiftAt=10, leaveAt=20, nSteps=30):
    """A run where H 3 moves from carbon 0 to carbon 1 at shiftAt and leaves at leaveAt."""
    coords = np.tile(POSITIONS, (nSteps, 1, 1))
    coords[shiftAt:, 3] = SHIFTED_H
    if leaveAt is not None:
        coords[leaveAt:, 3] = GONE_H
    # Small, different displacements so the frames can be told apart
    coords += np.arange(nSteps)[:, None, None] * 1e-4
    writeMol(os.path.join(folder, name), ATOMS, coords)
    return coords


def analyseFolder(folder, molecule='TOY'):
    for name in sorted(os.listdir(folder)):
        if name.endswith('.mol'):
            AutoDetectIsomerization(str(folder), name).autoAnalyseFile()
    GenerateChartsStats(str(folder), 'toy', molecule, 2, ['C', 'H', 'S'], None, False, False).autoAnalyseInFolder()


@pytest.fixture
def shiftFolder(tmp_path):
    coords = writeShiftRun(tmp_path, 'run_1.mol')
    writeShiftRun(tmp_path, 'run_2.mol', shiftAt=5, leaveAt=None)
    analyseFolder(tmp_path)
    return tmp_path, coords


def namesBySmiles(folder):
    return loadBlock(str(folder), 'isomerNames')


# ---------- connectivity keys ----------

def test_key_ignores_bond_orders_atom_order_fragment_order_and_bracket_hydrogens():
    assert stateKey('[H]C([H])=O') == stateKey('[H]C([H])O') == stateKey('O[C]([H])[H]')
    assert stateKey('[H]OC([H])([H])[H]') == stateKey('[H]C([H])([H])O[H]')
    assert stateKey('[H].[H]CC') == stateKey('[H]CC.[H]')
    assert stateKey('[SH]C') == stateKey('SC')
    # Ethanol and dimethyl ether have the same formula but not the same connectivity
    ethanol = '[H]OC([H])([H])C([H])([H])[H]'
    ether = '[H]C([H])([H])OC([H])([H])[H]'
    assert stateKey(ethanol) != stateKey(ether)


def test_largest_fragment_counts_heavy_atoms_first():
    assert largestFragment('[H][H].[H][H].CC') == 'CC'
    assert largestFragment('[H].[H]C([H])C') == '[H]C([H])C'


def test_state_keys_block_and_matching_across_folders_with_different_names(tmp_path):
    # In folder a the shifted state lasts longest, in folder b the start state, so their names swap
    first, second = tmp_path / 'a', tmp_path / 'b'
    first.mkdir(), second.mkdir()
    writeShiftRun(first, 'run_1.mol', shiftAt=5, leaveAt=None)
    writeShiftRun(second, 'run_1.mol', shiftAt=25, leaveAt=None)
    analyseFolder(first)
    analyseFolder(second)
    keysA = loadBlock(str(first), 'stateKeys')
    keysB = loadBlock(str(second), 'stateKeys')
    assert set(keysA) == set(keysB) == {'TOY_A1', 'TOY_A2'}
    assert keysA['TOY_A1']['key'] == keysB['TOY_A2']['key']
    assert keysA['TOY_A1']['formula'] == 'C2H2S'

    rows = matchStates([str(first), str(second)], ['a', 'b'])
    assert len(rows) == 2
    for row in rows:
        assert len(row['states']['a']) == len(row['states']['b']) == 1
        assert row['states']['a'][0][0] != row['states']['b'][0][0]


def test_matching_by_largest_fragment(shiftFolder):
    folder, _ = shiftFolder
    whole = matchStates([str(folder)], ['x'])
    largest = matchStates([str(folder)], ['x'], byLargestFragment=True)
    assert len(whole) == 3
    # The shifted state and the state after the H left have the same carbon skeleton but not the same H count,
    # so all three largest fragments differ, and the fragmented one is C2HS
    assert len(largest) == 3
    assert {row['formula'] for row in largest} == {'C2H2S', 'C2HS'}


# ---------- run sequences ----------

def test_run_sequences_block(shiftFolder):
    folder, _ = shiftFolder
    names = namesBySmiles(folder)
    sequences = loadBlock(str(folder), 'runSequences')
    assert list(sequences) == ['run_1.mol', 'run_2.mol']
    assert list(sequences['run_1.mol']) == ['0', '10', '20']
    assert list(sequences['run_2.mol']) == ['0', '5']
    assert set(sequences['run_1.mol'].values()) <= set(names.values())
    assert sequences['run_1.mol']['0'] == sequences['run_2.mol']['0']
    assert sequences['run_1.mol']['10'] == sequences['run_2.mol']['5']


# ---------- frames ----------

def test_visits():
    run = parseRun('In file: a.mol \nAt step 0:               A\nAt step 1:               A\n'
                   'At step 2: Isomerization B\n\nAt step 3:               B\nAt step 4: Isomerization A\n\n'
                   'Final molecules present: X')
    assert getVisits(run) == [('A', 0, 1), ('B', 2, 3), ('A', 4, 4)]


def readXyz(path):
    """Frames of an xyz file as (comment, elements, coordinates)."""
    lines = open(path).read().splitlines()
    frames = []
    index = 0
    while index < len(lines):
        count = int(lines[index])
        atoms = [line.split() for line in lines[index + 2:index + 2 + count]]
        frames.append((lines[index + 1], [a[0] for a in atoms], np.array([[float(v) for v in a[1:]] for a in atoms])))
        index += 2 + count
    return frames


def test_state_frame_is_the_middle_of_the_longest_visit(shiftFolder):
    folder, coords = shiftFolder
    sequences = loadBlock(str(folder), 'runSequences')
    shifted = sequences['run_1.mol']['10']
    paths = FrameExtractor(str(folder), shifted, maxRuns=1).extract()
    assert [os.path.basename(p) for p in paths] == [f'{shifted}__run_1_step14.xyz']
    [(comment, elements, xyz)] = readXyz(paths[0])
    assert comment == f'run_1.mol step 14 t=0.070 ps state {shifted}'
    assert elements == ATOMS
    assert np.allclose(xyz, coords[14], atol=1e-6)


def test_transition_frames_stay_inside_both_states(shiftFolder):
    folder, coords = shiftFolder
    sequence = loadBlock(str(folder), 'runSequences')['run_1.mol']
    start, shifted = sequence['0'], sequence['10']
    movie, reactant, product = FrameExtractor(str(folder), shifted, fromState=start, before=3, after=50,
                                              maxRuns=1).extract()
    frames = readXyz(movie)
    steps = [int(comment.split()[2]) for comment, _, _ in frames]
    assert steps == list(range(7, 20))   # after is cut at the end of the visit (step 19)
    assert frames[0][0].endswith(f'state {start}') and frames[-1][0].endswith(f'state {shifted}')
    assert np.allclose(readXyz(reactant)[0][2], coords[7], atol=1e-6)
    assert np.allclose(readXyz(product)[0][2], coords[19], atol=1e-6)


def test_largest_fragment_frames_drop_the_lost_hydrogen(shiftFolder):
    folder, _ = shiftFolder
    sequence = loadBlock(str(folder), 'runSequences')['run_1.mol']
    lost = sequence['20']
    [path] = FrameExtractor(str(folder), lost, largestFragment=True).extract()
    assert readXyz(path)[0][1] == ['C', 'C', 'S', 'H']
    # Fragmentation transition: the union of the first and last largest fragments keeps the leaving H
    movie = FrameExtractor(str(folder), lost, fromState=sequence['10'], largestFragment=True).extract()[0]
    assert all(len(elements) == 5 for _, elements, _ in readXyz(movie))


def test_frames_without_the_mol_file_or_with_unknown_state(shiftFolder):
    folder, _ = shiftFolder
    with pytest.raises(ValueError, match='not a state name'):
        FrameExtractor(str(folder), 'TOY_Q9').extract()
    os.remove(folder / 'run_1.mol')
    lost = loadBlock(str(folder), 'runSequences')['run_1.mol']['20']
    with pytest.warns(UserWarning, match='missing'):
        assert FrameExtractor(str(folder), lost).extract() == []


def test_frames_warn_when_settings_do_not_match_the_log(shiftFolder):
    folder, _ = shiftFolder
    logPath = folder / 'run_1_log.txt'
    text = open(logPath).read().replace('"break_distance": 2.5', '"break_distance": 2.5, "x": 1')
    # Pretend the run was analysed with a formation distance so short that no bond forms
    text = text.replace('"formation_distance": 1.63', '"formation_distance": 0.5')
    open(logPath, 'w').write(text)
    shifted = loadBlock(str(folder), 'runSequences')['run_1.mol']['10']
    with pytest.warns(UserWarning, match='bonds of the logged state'):
        FrameExtractor(str(folder), shifted, maxRuns=1).extract()


# ---------- report ----------

def test_report_lists_states_and_embeds_plots(shiftFolder):
    folder, _ = shiftFolder
    path = writeReport(str(folder))
    text = open(path, encoding='utf-8').read()
    for name in namesBySmiles(folder).values():
        assert f'<b>{name}</b>' in text
    assert text.count('<svg') == 3
    assert text.count('data:image/png;base64,') == len([f for f in os.listdir(folder / 'plots') if f.endswith('.png')])
    assert '50%' in text   # the fragment was present at the end of one of the two runs
    # Without run_summary.json the blocks are read from full_log.txt
    os.remove(folder / 'run_summary.json')
    assert writeReport(str(folder), maxStates=1, fileName='small.html')
    assert open(folder / 'small.html', encoding='utf-8').read().count('<svg') == 1


def test_report_needs_statistics(tmp_path):
    with pytest.raises(ValueError, match='makechartsandstats'):
        writeReport(str(tmp_path))


# ---------- breakdown run ----------

def test_config_round_trip_and_errors(tmp_path):
    config = defaultConfig()
    config['analysis']['elements'] = ['C', 'H', 'S']
    config['statistics']['molecule_name'] = 'a "quoted" name'
    path = tmp_path / 'c.toml'
    path.write_text(formatConfig(config))
    assert loadConfig(str(path)) == config
    path.write_text('[analysis]\nbreak_distance = 3\n')
    assert loadConfig(str(path))['analysis']['break_distance'] == 3.0
    for text, message in [('[analysis]\nbreakdistance = 3\n', 'Unknown option'),
                          ('[analyse]\n', 'Unknown section'),
                          ('[analysis]\nmin_lifetime = "2"\n', 'must be a int'),
                          ('[analysis\n', 'not valid TOML')]:
        path.write_text(text)
        with pytest.raises(ConfigError, match=message):
            loadConfig(str(path))


def test_run_pipeline_end_to_end_and_reuse_of_logs(tmp_path):
    writeShiftRun(tmp_path, 'run_1.mol')
    writeShiftRun(tmp_path, 'run_2.mol', shiftAt=5, leaveAt=None)
    config = defaultConfig()
    config['analysis']['jobs'] = 1
    config['statistics']['molecule_name'] = 'TOY'
    used = runPipeline(str(tmp_path), config)
    for name in ('full_log.txt', 'run_summary.json', 'report.html', 'breakdown_used.toml', 'run_1_log.txt'):
        assert (tmp_path / name).exists(), name
    assert used['statistics']['chart_name'] == tmp_path.name
    # The state at step 0, not the most common one (here the shifted state lasts longer and is TOY_A1)
    assert used['flowcharts']['start_state'] == loadBlock(str(tmp_path), 'runSequences')['run_1.mol']['0'] == 'TOY_A2'
    assert loadConfig(str(tmp_path / 'breakdown_used.toml')) == used

    logTime = os.path.getmtime(tmp_path / 'run_1_log.txt')
    runPipeline(str(tmp_path), config)
    assert os.path.getmtime(tmp_path / 'run_1_log.txt') == logTime

    config['analysis']['break_distance'] = 3.0
    runPipeline(str(tmp_path), config)
    settings = json.loads(open(tmp_path / 'run_1_log.txt').read().splitlines()[1].split(': ', 1)[1])
    assert settings['break_distance'] == 3.0


def test_run_pipeline_without_mol_files(tmp_path):
    with pytest.raises(ConfigError, match='No files matching'):
        runPipeline(str(tmp_path), defaultConfig())
