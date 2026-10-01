import os
import matplotlib
matplotlib.use('Agg')
import pytest
from breakdown.auto_detect_utils import classifyRun
from breakdown.dissociation_time_chart import DissociationTimeChart
from breakdown.fragments_bar_chart import FragmensBarChart
from breakdown.fragments_pie_chart import FragmentsPieChart
from breakdown.functional_group_comparison import FunctionalGroupComparisonChart
from breakdown.isomer_frag_pie_chart import IsomerFragPieChart
from breakdown.oxygen_fate_pie_chart import OxygenFatePieChart
from breakdown.plot_style import rangeFormulaLabel


def log(*lines):
    return '\n'.join(lines)


def test_recombined_run_is_not_counted_as_fragmented():
    recombined = log('At step 0:               CC', 'At step 1: Fragmentation C.C', 'At step 2: Recombination CC',
                     'At step 3:               CC', 'Final molecules present: C2')
    fragmented = log('At step 0:               CC', 'At step 1: Fragmentation C.C', 'At step 2:               C.C',
                     'Final molecules present: C, C')
    isomerized = log('At step 0:               CC', 'At step 1: Isomerization C=C', 'Final molecules present: C2')
    assert classifyRun(recombined) == (False, True, False)
    assert classifyRun(fragmented) == (False, True, True)
    assert classifyRun(isomerized) == (True, False, False)
    both = log('At step 0: CC', 'At step 1: Isomerization C=C', 'At step 2: Fragmentation C.C')
    assert classifyRun(both) == (True, True, True)


def test_dissociation_legend_format():
    assert rangeFormulaLabel(6, [6, 7], {'O': [0, 1]}) == '$\\mathrm{C_{6}H_{6–7}O_{0–1}}$'
    assert rangeFormulaLabel(6, [6], {'S': [1]}) == '$\\mathrm{C_{6}H_{6}S_{1}}$'
    assert rangeFormulaLabel(10, [8], {'O': [0], 'Cl': [2]}) == '$\\mathrm{C_{10}H_{8}Cl_{2}}$'  # never seen: left out


def test_charts_write_pdf_and_png(tmp_path):
    folder = str(tmp_path)
    isoFrag = IsomerFragPieChart('a', folder)
    isoFrag.addIsomerizedAndFragmented()
    isoFrag.addNotIsomerizedOrFragmented()
    pie = FragmentsPieChart('a', folder)
    pie.addFragment('H')
    pie.addFragment('No fragmentation')
    bar = FragmensBarChart('a', folder)
    bar.addFragments(['C6H6', 'H'])
    dissociation = DissociationTimeChart('a', folder, 6, heteroatoms=['O'])
    dissociation.addRun([6, 6, 4, 4])
    dissociation.addRun([6, 6, 6])
    dissociation.addToHydrogenArray(6, 6)
    dissociation.addToHydrogenArray(6, 5)
    dissociation.addToOxygenArray(6, 1)
    for chart, method in [(isoFrag, 'createIsomerFragPieChart'), (pie, 'createFragmentsPieChart'),
                          (bar, 'createFragmentsBarChart'), (dissociation, 'createDissociationTimeChart')]:
        path = getattr(chart, method)()
        assert path.endswith('.pdf') and os.path.exists(path) and os.path.exists(path[:-3] + 'png')
    assert dissociation.getDissociationArray()[3].tolist()[4] == 50.0
    assert dissociation.makeLegendString(6) == '$\\mathrm{C_{6}H_{5–6}O_{1}}$'


def test_empty_charts_do_nothing(tmp_path):
    folder = str(tmp_path)
    assert IsomerFragPieChart('x', folder).createIsomerFragPieChart() is None
    assert FragmentsPieChart('x', folder).createFragmentsPieChart() is None
    assert FragmensBarChart('x', folder).createFragmentsBarChart() is None
    assert FunctionalGroupComparisonChart('x', folder).createComparisonPlot() is None
    assert OxygenFatePieChart('x', folder).createPiePlots() is None
    assert not os.listdir(tmp_path)


def test_more_than_eight_series_and_one_molecule_do_not_crash(tmp_path):
    chart = FunctionalGroupComparisonChart('many', str(tmp_path))
    for i in range(10):
        chart.processSeries(f'M{i}', 'ketone', {'0': 10, '1': 9, '2': 4})
    assert os.path.exists(chart.createComparisonPlot())
    pies = OxygenFatePieChart('one', str(tmp_path))
    pies.processFateData('A', 'ketone', {'ketone': 3, 'hydroxyl': 2})
    assert os.path.exists(pies.createPiePlots())


def test_nested_pie_labels_follow_the_wedge_order(tmp_path, monkeypatch):
    chart = IsomerFragPieChart('x', str(tmp_path))
    chart.numIsomerizedAndFragmented, chart.numIsomerizedNotFragmented = 1, 2
    chart.numNotIsomerizedOrFragmented, chart.numNotIsomerizedButFragmented = 3, 4
    calls = []
    original = matplotlib.axes.Axes.pie

    def recordingPie(self, values, *args, **kwargs):
        calls.append((list(values), list(kwargs['labels'])))
        return original(self, values, *args, **kwargs)
    monkeypatch.setattr(matplotlib.axes.Axes, 'pie', recordingPie)
    chart.createIsomerFragPieChart()
    values, labels = calls[1]  # the inner ring
    assert values == labels == [1, 2, 3, 4]  # every wedge carries its own number
