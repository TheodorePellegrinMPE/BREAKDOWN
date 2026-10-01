import json
import os
import pytest
from breakdown.auto_detect_isomerization import AutoDetectIsomerization
from breakdown.generate_charts_stats import GenerateChartsStats
from breakdown.elemental_stability_chart import ElementalStabilityChart
from breakdown.event_survival_chart import (FIRST_EVENTS_HEADING, generateEventComparison,
                                                loadFirstEvents)
from breakdown.run_summary import (getSurvivalCurve, kaplanMeier, loadSection, loadSummary,
                                       readSectionFromFullLog)
from test_pipeline import makeRun


def test_kaplan_meier_with_censoring():
    curve = kaplanMeier([1, 2, 3, 4], [True, True, False, True])
    assert curve['times'] == [0.0, 1, 2, 4]
    assert curve['survival'] == pytest.approx([1, 0.75, 0.5, 0.0])
    assert curve['median'] == 2
    assert all(l <= s <= u for l, s, u in zip(curve['lower'], curve['survival'], curve['upper']))


def test_no_events_means_no_median():
    curve = kaplanMeier([5, 5, 5], [False, False, False])
    assert curve['survival'] == [1.0] and curve['median'] is None


def test_survival_curve_uses_last_step_for_runs_without_event():
    firstEvents = {'a': {'fragmentation': 100, 'lastStep': 1000},
                   'b': {'fragmentation': None, 'lastStep': 1000}}
    curve = getSurvivalCurve(firstEvents, 'fragmentation', timePerStepFs=5.0)
    assert curve['times'] == [0.0, 0.5]  # 100 steps * 5 fs = 0.5 ps
    assert curve['survival'] == pytest.approx([1.0, 0.5])
    assert (curve['numEvents'], curve['numRuns']) == (1, 2)


@pytest.fixture
def analysedFolder(tmp_path):
    for i, leaves in enumerate([10, None, 20, 25]):
        AutoDetectIsomerization(str(tmp_path), makeRun(tmp_path, f'run_{i}.mol', leaves)).autoAnalyseFile()
    GenerateChartsStats(str(tmp_path), 'toy', 'TOY', 2, ['C', 'H', 'S'], 'thiol', False, False).autoAnalyseInFolder()
    return tmp_path


def test_summary_is_written_and_matches_full_log(analysedFolder):
    summary = loadSummary(str(analysedFolder))
    assert summary['numRuns'] == 4 and summary['timePerStepFs'] == 5.0
    assert (analysedFolder / 'plots' / 'toy_first_event_survival.pdf').exists()
    for heading, block in summary['sections'].items():
        assert readSectionFromFullLog(str(analysedFolder), heading) == block, heading


def test_first_events_content(analysedFolder):
    firstEvents, timePerStep = loadFirstEvents(str(analysedFolder))
    assert timePerStep == 5.0
    assert firstEvents['run_0.mol']['fragmentation'] == 10
    assert firstEvents['run_1.mol']['anyChange'] is None
    assert firstEvents['run_3.mol']['lastStep'] == 29


def test_stale_or_missing_summary_falls_back_to_full_log(analysedFolder):
    folder = str(analysedFolder)
    fromSummary = loadSection(folder, FIRST_EVENTS_HEADING)
    summaryPath = analysedFolder / 'run_summary.json'
    # A summary older than full_log.txt is ignored
    old = os.path.getmtime(analysedFolder / 'full_log.txt') - 100
    os.utime(summaryPath, (old, old))
    assert loadSummary(folder) is None
    assert loadSection(folder, FIRST_EVENTS_HEADING) == fromSummary
    os.remove(summaryPath)
    assert loadSection(folder, FIRST_EVENTS_HEADING) == fromSummary
    assert loadSection(folder, 'No such heading.') is None


def test_elemental_series_same_from_summary_and_log(analysedFolder):
    folder = str(analysedFolder)
    chart = ElementalStabilityChart('x', str(analysedFolder), heteroatom='S')
    fast = chart.loadDissociationData(folder)
    os.remove(analysedFolder / 'run_summary.json')
    slow = chart.loadDissociationData(folder)
    assert fast == slow
    assert fast['cIntact'][0][1] == 4 and fast['cIntact'][-1][1] == 4  # only an H leaves


def test_comparison_needs_first_event_data(analysedFolder, tmp_path_factory):
    out = str(tmp_path_factory.mktemp('out'))
    path = generateEventComparison([str(analysedFolder)] * 2, ['a', 'b'], 'fragmentation', 't', out)
    assert os.path.exists(path)
    with pytest.raises(ValueError):
        generateEventComparison([out], ['x'], 'fragmentation', 't', out)
