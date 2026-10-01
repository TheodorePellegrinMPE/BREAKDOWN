---
name: extending_the_pipeline
description: >
  Checklist of the files to change when adding a functional group, a statistic to full_log.txt, a plot, a CLI command
  or a new event type to breakdown. Use when implementing any new analysis feature.
---

# Extending breakdown

## Add a functional group
1. `isomer_analyser.py`: write `identifyX(mol)` (connectivity only, hydrogens are explicit atoms, use `atom.GetNeighbors()`),
   add `(function, "Label")` to the list in `analyseIsomer`. The label becomes `Label count: N.` in the log.
   Heteroatom variants take an `element` argument (see `partial(identifyKetone, element='S')`).
2. `isomer_dicts.py`: add the lowercase name to `OXYGEN_GROUPS` or `SULFUR_GROUPS` so `endingGroups` counts it.
   Patterns use `\b` so `ether` does not match `thioether`, keep names that are not suffixes of each other.
3. Add the molecule to `MOLECULES` in `tests/test_diverse_molecules.py`.
4. Update the field list in `.agent/skills/full_log_interpreter/SKILL.md`.

## Add a statistic to full_log.txt
1. Compute it in `IsomerDicts` from `self.runs` (parsed once: per run a list of `(step, event, smiles, analysis)` and `final`).
   Do not re-scan `allLogString`.
2. Call `createX()` in `createIsomerDicts`, add `(heading, 'attr')` to the `blocks` list in
   `GenerateChartsStats.addStatsToString`. That writes it to `full_log.txt` and to `run_summary.json` (key `sections`).
   `getFormattedStats` renames SMILES keys to isomer names. Read blocks in consumers with `run_summary.loadSection`.
   Values must be JSON serialisable, the heading line must end with a period and the JSON must follow on the next lines.
3. Consumers find blocks by heading text (`loadSection`), so never rename an
   existing heading. Document the new block in the full log skill.

## Add a plot or command
- Plot class in its own module. Draw inside `with textStyle():` (fonts, file settings, never `rcParams` at import) and finish with
  `savePdfAndPng(folder, name)` (writes pdf + png, closes the figure). Return the pdf path, or `None` when there is no data.
  Do not redesign the existing charts (nested pie, pies, bar, line plots): the author chose those on purpose, only fix bugs and text sizes.
  Check that labels are given in the same order as the values (the nested pie had two swapped labels). `event_survival_chart.py` shows the
  fuller `plotStyle()` look for new charts.
  Take `timePerStepFs` (default 5.0 fs: 0.5 fs time step, output every 10 steps) instead of hard-coding time or the number of runs.
- Command in `main.py`: thin wrapper. Comparison commands take semicolon separated `--folders/--names`.
  Options are click, names auto-lowercased (`makeFlowCharts` -> `makeflowcharts`).
- Streaming beats loading: `full_log.txt` can be 100 MB, iterate `for line in file` and stop at
  `Automatically generated molecule names.` when only run data is needed (see `elemental_stability_chart.py`).

## Add an event type
Add the label (13 characters, same width as `Isomerization`) to `EVENT_LABELS` in `auto_detect_utils.py`, produce it in
`checkIsomerization`, and check every place that tests for labels: `isomer_dicts.parseRuns` (uses `EVENT_LABELS`),
`generate_charts_stats.addToIsomerizationCharts`, `auto_detect_utils.classifyRun`, the full log skill.

## Add an element
`VALENCES` in `molecule_utils.py`, `NOMINAL_MASSES` in `auto_detect_utils.py` (falls back to RDKit isotope masses), formation
distance overrides in `bond_adder.getPairThresholds` only if the default 1.63 Å is clearly wrong for that pair.

## Finally
Run `pytest` and the regression recipe in `testing_and_regression`.
