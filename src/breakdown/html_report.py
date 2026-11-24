"""report.html: one page per folder with the states, transitions, final fragments, settings and plots.

The page is a single file (drawings are inline SVG, plots are embedded PNGs) so it can be copied or mailed on its
own. It only reads the statistics blocks (run_summary.json, or full_log.txt for older results) and the PNGs in
plots/, so it is quick to make again.
"""
import base64
import glob
import html
import os
from rdkit.Chem.Draw import rdMolDraw2D
from breakdown.molecule_utils import moleculeFromSmiles
from breakdown.run_summary import loadBlock, loadSummary, resolveTimePerStep
from breakdown.state_keys import loadStateKeys

REPORT_FILE = 'report.html'
TOP_LINKS = 3
TOP_TRANSITIONS = 40
DRAWING_SIZE = (220, 160)

STYLE = """
:root { color-scheme: light dark; --fg:#1d1d1f; --muted:#6b6b70; --bg:#fafafa; --card:#fff; --line:#e3e3e6;
        --accent:#2a6fdb; }
@media (prefers-color-scheme: dark) { :root { --fg:#ececef; --muted:#a0a0a8; --bg:#18181a; --card:#222225;
        --line:#38383d; --accent:#7aa7ff; } }
body { font: 14px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif; color: var(--fg); background: var(--bg);
       margin: 0; padding: 24px 32px 64px; }
h1 { font-size: 22px; margin: 0 0 4px; } h2 { font-size: 17px; margin: 32px 0 8px; }
.muted { color: var(--muted); } a { color: var(--accent); }
.facts { display: flex; flex-wrap: wrap; gap: 8px 28px; margin: 12px 0; }
.facts div b { display: block; font-size: 18px; font-variant-numeric: tabular-nums; }
table { border-collapse: collapse; background: var(--card); border: 1px solid var(--line); }
th, td { padding: 6px 10px; border-bottom: 1px solid var(--line); text-align: left; vertical-align: middle; }
th { position: sticky; top: 0; background: var(--card); cursor: pointer; user-select: none; white-space: nowrap; }
th[data-dir="asc"]::after { content: " \\25B2"; } th[data-dir="desc"]::after { content: " \\25BC"; }
td.num { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
td.links { font-size: 12px; } code { font-size: 12px; }
.drawing { background: #fff; border-radius: 6px; display: inline-block; line-height: 0; }
.drawing svg { width: 165px; height: 120px; }
input[type=search] { font: inherit; padding: 6px 10px; width: 320px; margin-bottom: 8px; border: 1px solid var(--line);
       border-radius: 6px; background: var(--card); color: var(--fg); }
.plots { display: grid; grid-template-columns: repeat(auto-fill, minmax(420px, 1fr)); gap: 16px; }
.plots figure { margin: 0; background: #fff; border: 1px solid var(--line); border-radius: 6px; padding: 8px; }
.plots img { width: 100%; } .plots figcaption { color: #555; font-size: 12px; }
pre { background: var(--card); border: 1px solid var(--line); padding: 10px; overflow-x: auto; font-size: 12px; }
"""

SCRIPT = """
document.querySelectorAll('table.sortable').forEach(table => {
  table.querySelectorAll('th').forEach((th, column) => th.addEventListener('click', () => {
    const dir = th.dataset.dir === 'desc' ? 'asc' : 'desc';
    table.querySelectorAll('th').forEach(other => delete other.dataset.dir);
    th.dataset.dir = dir;
    const body = table.tBodies[0];
    const value = row => { const cell = row.cells[column]; return cell.dataset.sort ?? cell.textContent.trim(); };
    const rows = [...body.rows].sort((a, b) => {
      const x = value(a), y = value(b), nx = parseFloat(x), ny = parseFloat(y);
      const order = (!isNaN(nx) && !isNaN(ny)) ? nx - ny : x.localeCompare(y);
      return dir === 'asc' ? order : -order;
    });
    rows.forEach(row => body.appendChild(row));
  }));
});
const filter = document.getElementById('stateFilter');
if (filter) filter.addEventListener('input', () => {
  const words = filter.value.toLowerCase().split(/\\s+/).filter(Boolean);
  document.querySelectorAll('#states tbody tr').forEach(row => {
    const text = row.textContent.toLowerCase();
    row.hidden = !words.every(word => text.includes(word));
  });
});
"""


def drawSvg(smiles):
    """SVG drawing of a state, an empty string if it cannot be drawn."""
    try:
        mol = moleculeFromSmiles(smiles)
        drawer = rdMolDraw2D.MolDraw2DSVG(*DRAWING_SIZE)
        drawer.drawOptions().addStereoAnnotation = False
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg[svg.find('<svg'):]
    except Exception:
        return ''


def formatLinks(counts, net, limit=TOP_LINKS):
    """Text of the most common linked states, like 'A2 (12, net 8)'."""
    items = sorted(((name, count) for name, count in counts.items() if name), key=lambda item: -item[1])
    parts = [f'{html.escape(name)} ({count}, net {net.get(name, 0)})' for name, count in items[:limit]]
    if len(items) > limit:
        parts.append(f'+{len(items) - limit} more')
    return ', '.join(parts)


def numberCell(value, text=None):
    return f'<td class="num" data-sort="{value}">{html.escape(text if text is not None else str(value))}</td>'


def buildStatesTable(names, keys, runs, steps, transitions, net, numRuns, timePerStepFs, maxStates):
    """Rows of the states table, the most common states by number of runs."""
    children, netChildren = {}, {}
    for target, sources in transitions.items():
        for source, count in sources.items():
            children.setdefault(source, {})[target] = count
    for target, sources in net.items():
        for source, count in sources.items():
            netChildren.setdefault(source, {})[target] = count
    smilesOf = {name: smiles for smiles, name in names.items()}
    ordered = sorted(smilesOf, key=lambda name: (-runs.get(name, 0), -steps.get(name, 0), name))[:maxStates]
    psPerStep = timePerStepFs / 1000
    rows = []
    for name in ordered:
        info = keys.get(name, {})
        nRuns = runs.get(name, 0)
        totalPs = steps.get(name, 0) * psPerStep
        meanPs = totalPs / nRuns if nRuns else 0
        rows.append(
            '<tr>'
            f'<td><span class="drawing">{drawSvg(smilesOf[name])}</span></td>'
            f'<td><b>{html.escape(name)}</b></td>'
            f'<td>{html.escape(info.get("formula", ""))}</td>'
            + numberCell(nRuns, f'{nRuns} ({100 * nRuns / numRuns:.0f}%)' if numRuns else str(nRuns))
            + numberCell(round(totalPs, 3), f'{totalPs:.2f}')
            + numberCell(round(meanPs, 3), f'{meanPs:.2f}')
            + f'<td class="links">{formatLinks(transitions.get(name, {}), net.get(name, {}))}</td>'
            f'<td class="links">{formatLinks(children.get(name, {}), netChildren.get(name, {}))}</td>'
            f'<td><code>{html.escape(info.get("key", ""))}</code></td>'
            '</tr>')
    return rows, len(smilesOf)


def buildTransitionRows(transitions, limit=TOP_TRANSITIONS):
    """Rows of the most common transitions with their reverse and net counts."""
    pairs = [(source, target, count) for target, sources in transitions.items()
             for source, count in sources.items() if source]
    pairs.sort(key=lambda pair: -pair[2])
    rows = []
    for source, target, count in pairs[:limit]:
        reverse = transitions.get(source, {}).get(target, 0)
        rows.append(f'<tr><td>{html.escape(source)}</td><td>{html.escape(target)}</td>'
                    + numberCell(count) + numberCell(reverse) + numberCell(count - reverse) + '</tr>')
    return rows, len(pairs)


def embedPlots(folder):
    """Figures of every PNG in folder/plots."""
    figures = []
    for path in sorted(glob.glob(os.path.join(folder, 'plots', '*.png'))):
        with open(path, 'rb') as file:
            data = base64.b64encode(file.read()).decode()
        name = os.path.basename(path)
        figures.append(f'<figure><img alt="{html.escape(name)}" src="data:image/png;base64,{data}">'
                       f'<figcaption>{html.escape(name)}</figcaption></figure>')
    return figures


def linkFiles(folder, pattern):
    """Links, relative to the folder, to the files matching a pattern."""
    paths = sorted(glob.glob(os.path.join(folder, pattern)))
    return [f'<a href="{html.escape(os.path.relpath(path, folder))}">{html.escape(os.path.relpath(path, folder))}</a>'
            for path in paths]


def writeReport(folder, maxStates=100, fileName=REPORT_FILE):
    """Write the report of a folder.

    Args:
        folder (str): Folder with full_log.txt (and plots/).
        maxStates (int): Number of states in the states table, the most common by number of runs.
        fileName (str): Name of the file in the folder.

    Returns:
        str: Path of the report.

    Raises:
        ValueError: Raised if the folder has no statistics.
    """
    names = loadBlock(folder, 'isomerNames')
    if names is None:
        raise ValueError(f'{folder} has no full_log.txt with statistics, run makechartsandstats first.')
    runs = loadBlock(folder, 'isomerRuns') or {}
    steps = loadBlock(folder, 'isomerStats') or {}
    transitions = loadBlock(folder, 'isomerizationStats') or {}
    net = loadBlock(folder, 'netIsomerizationStats') or {}
    endFragments = loadBlock(folder, 'endFragmentStats') or {}
    firstEvents = loadBlock(folder, 'firstEvents') or {}
    settings = loadBlock(folder, 'settings') or {}
    keys = loadStateKeys(folder) or {}
    summary = loadSummary(folder) or {}
    numRuns = summary.get('numRuns') or len(firstEvents) or max(runs.values(), default=0)
    timePerStepFs = resolveTimePerStep([folder])
    molecule = summary.get('molecule') or settings.get('statistics', {}).get('molecule_name') or ''
    title = f'BREAKDOWN report: {molecule or os.path.basename(os.path.abspath(folder))}'

    stateRows, numStates = buildStatesTable(names, keys, runs, steps, transitions, net, numRuns, timePerStepFs,
                                            maxStates)
    transitionRows, numTransitions = buildTransitionRows(transitions)
    events = {event: sum(1 for record in firstEvents.values() if record.get(event) is not None)
              for event in ('anyChange', 'isomerization', 'fragmentation')}
    fragmentRows = [f'<tr><td>{html.escape(formula)}</td>' + numberCell(count)
                    + numberCell(round(100 * count / numRuns, 1) if numRuns else 0,
                                 f'{100 * count / numRuns:.0f}%' if numRuns else '') + '</tr>'
                    for formula, count in sorted(endFragments.items(), key=lambda item: -item[1])]

    analysis = settings.get('analysis') or []
    program = settings.get('program') or {}
    settingsText = '\n'.join(
        [f'Analysis settings: {setting}' for setting in analysis]
        + ([f'WARNING: the runs were analysed with {len(analysis)} different settings.'] if len(analysis) > 1 else [])
        + ([f'{settings["runsWithoutRecordedSettings"]} runs have no recorded settings (older version).']
           if settings.get('runsWithoutRecordedSettings') else [])
        + [f'Statistics options: {settings.get("statistics", {})}',
           f'Program: {program}', f'Made: {settings.get("created", "unknown")}'])

    def fact(value, label):
        return f'<div><b>{html.escape(str(value))}</b><span class="muted">{html.escape(label)}</span></div>'

    flowLinks = linkFiles(folder, 'transitionGraphs/*.pdf')
    frameLinks = linkFiles(folder, 'frames/*.xyz')
    parts = [
        '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">',
        f'<meta name="viewport" content="width=device-width, initial-scale=1"><title>{html.escape(title)}</title>',
        f'<style>{STYLE}</style></head><body>',
        f'<h1>{html.escape(title)}</h1>',
        f'<div class="muted">{html.escape(os.path.abspath(folder))}</div>',
        '<div class="facts">',
        fact(numRuns, 'runs'),
        fact(numStates, 'states'),
        fact(events['anyChange'], 'runs that changed'),
        fact(events['isomerization'], 'runs that isomerized'),
        fact(events['fragmentation'], 'runs that fragmented'),
        fact(f'{timePerStepFs:g} fs', 'per step'),
        '</div>',
        f'<h2>States</h2><p class="muted">The {len(stateRows)} most common of {numStates} states by number of runs. '
        'Time is summed over all runs; "per run" divides it by the runs that visited the state. "Formed from" and '
        '"Turns into" give total transitions and net transitions (forward minus back). The key is the same for the same '
        'connectivity in every folder (<code>breakdown matchstates</code>). Click a heading to sort.</p>',
        '<input type="search" id="stateFilter" placeholder="Filter, e.g. a name, formula or key">',
        '<table class="sortable" id="states"><thead><tr><th>Structure</th><th>Name</th><th>Formula</th><th>Runs</th>'
        '<th>Time (ps)</th><th>Per run (ps)</th><th>Formed from</th><th>Turns into</th><th>Key</th></tr></thead>'
        '<tbody>' + ''.join(stateRows) + '</tbody></table>',
        f'<h2>Most common transitions</h2><p class="muted">{len(transitionRows)} of {numTransitions} transitions. '
        'Net close to zero with large counts means flickering back and forth, not a route.</p>',
        '<table class="sortable"><thead><tr><th>From</th><th>To</th><th>Count</th><th>Back</th><th>Net</th></tr></thead>'
        '<tbody>' + ''.join(transitionRows) + '</tbody></table>',
        '<h2>Fragments at the end of the runs</h2>',
        '<table class="sortable"><thead><tr><th>Fragment</th><th>Runs</th><th>Share</th></tr></thead><tbody>'
        + ''.join(fragmentRows) + '</tbody></table>',
    ]
    if flowLinks:
        parts.append('<h2>Flow charts</h2><p>' + ' · '.join(flowLinks) + '</p>')
    if frameLinks:
        parts.append('<h2>Extracted MD frames</h2><p>' + ' · '.join(frameLinks) + '</p>')
    plots = embedPlots(folder)
    if plots:
        parts.append('<h2>Plots</h2><div class="plots">' + ''.join(plots) + '</div>')
    parts += ['<h2>Settings</h2>', f'<pre>{html.escape(settingsText)}</pre>',
              f'<script>{SCRIPT}</script></body></html>']
    path = os.path.join(folder, fileName)
    with open(path, 'w', encoding='utf-8') as file:
        file.write('\n'.join(parts))
    return path
