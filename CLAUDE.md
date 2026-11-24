# BREAKDOWN – guide for agents

BREAKDOWN = Bond Rearrangement, Events And Kinetics: Dissociation Output from Whole-trajectory aNalysis (formerly autoastrochem).

Poetry project (`src/breakdown`, CLI `breakdown`, click 8 lowercases command names:
`analysedemon`, `makechartsandstats`, ...). It turns deMon MD trajectories into SMILES logs, statistics and plots.
Read `README.md` for the workflow, and the skills in `.agent/skills/` before changing anything:

- `testing_and_regression` – environment, example data, how to check nothing changed. **Read first.**
- `smiles_and_bond_orders` – why the SMILES strings look the way they do, and the RDKit implicit hydrogen trap.
- `extending_the_pipeline` – which files to touch to add a functional group, statistic, plot or command.
- `full_log_interpreter` – format of `full_log.txt`. These files are huge (millions of lines), never read them directly.

## Pipeline

```
.mol trajectory --analysedemon--> X_log.txt (one per run) --makechartsandstats--> full_log.txt + plots/
full_log.txt --> makeflowcharts / comparefunctionalgroups / compareelementalstability / plotoxygenfates
full_log.txt --> makereport (report.html) / matchstates (across folders) / extractframes (+ .mol -> frames/*.xyz)
breakdown run FOLDER = analysedemonfolder + makechartsandstats + flow charts + makereport, options from breakdown.toml
```

| Module | Role |
|---|---|
| `main.py` | click commands only, no logic |
| `run_log.py` | the run log format (`_log.txt`, and the runs at the start of `full_log.txt`): `formatRunLog`/`formatStepLine` write it, `parseRuns` reads it into `Run`s (`steps`, `final`, `settings`, `isomerized`, `fragmentedAtEnd`), `iterSteps` streams it |
| `auto_detect_isomerization.py` | reads a `.mol` file into NumPy arrays, writes the `_log.txt` |
| `bond_adder.py` | vectorised bond finding with formation/break hysteresis, one bond per hydrogen |
| `create_smiles_from_dataframe.py` | bonds -> SMILES per step (cached per bond pattern), event labels, min lifetime |
| `isomer_analyser.py` | functional groups of a state (ring sizes, hydroxyl, ketone, thiol, ...) |
| `molecule_utils.py` | `moleculeFromSmiles`: SMILES -> sanitized RDKit mol without extra hydrogens |
| `generate_charts_stats.py`, `isomer_dicts.py` | combine logs into `full_log.txt` (read runs -> statistics -> write -> charts), all statistics dictionaries |
| `isomer_frag_pie_chart.py`, `fragments_pie_chart.py`, `fragments_bar_chart.py`, `dissociation_time_chart.py`, `functional_group_comparison.py`, `elemental_stability_chart.py`, `oxygen_fate_pie_chart.py` | the author's chart designs (nested pie, pies, line plots). Do not change what kind of chart they draw, only fix bugs and adjust text sizes/formats. The per-folder charts take one `Run` at a time (`addRun(run)`) and do their own counting |
| `plot_style.py`, `event_survival_chart.py` | `textStyle()` (fonts) and `savePdfAndPng` used by all charts so every figure is a pdf + png; the full `plotStyle()` is only used by the newer survival chart |
| `molecule_transition_graph.py`, `transition_best_path.py` | flow charts (graphviz, needs `dot`) |
| `run_summary.py`, `event_survival_chart.py` | the statistics blocks (`BLOCKS`: key -> heading), `formatBlocks` writes them, `loadBlock(folder, key)` reads one from `run_summary.json` or streams `full_log.txt`; Kaplan-Meier time-to-first-event statistics and plots |
| `auto_detect_utils.py` | shared helpers, cached atom counting |
| `state_keys.py` | connectivity keys (`stateKey`: hash of the bond graph, ignores bond orders/atom order), the `stateKeys` block, `matchStates` across folders (`loadStateKeys` computes keys from `isomerNames` for older results) |
| `frame_extractor.py` | `extractframes`: MD frames of a state or around a transition from the `.mol` file (step = frame index), bonds re-found with the run's settings and checked against the log |
| `html_report.py` | `makereport`: single-file `report.html` from the blocks (via `loadBlock`) and `plots/*.png` |
| `pipeline.py` | `breakdown run`: TOML config (`CONFIG_OPTIONS` is the single list of options and defaults, unknown keys are errors), `analyseFiles` (the parallel analysis also used by `analysedemonfolder`), keeps logs whose recorded settings match |

`cluster_helpers/` holds the author's scripts for running many simulations on a cluster (settings in `cluster_config.py`, which is git-ignored; copy
`cluster_config.example.py`). They are plain scripts with `--dry-run`, see README section 11. Test them with stub `ssh`/`rsync` commands in `PATH`,
never against the real cluster. Every run log has a `Settings:` line, see `run_settings.py`.

## Rules that are easy to break

- Log format is an interface: `At step N: <13 char event> <smiles>`, the analysis line directly after an event line,
  `Final molecules present: ...`. Only `run_log.py` writes or parses it, everything else works on `Run` objects.
  Never parse lines elsewhere. Old logs must stay readable, update the full log skill doc when the format changes.
- Never read SMILES with `Chem.MolFromSmiles` when you need a molecule: use `moleculeFromSmiles`.
- Do not put per-step Python loops over pandas rows back into the bond finding, it was 100x slower.
- Consumers of `full_log.txt` blocks must use `run_summary.loadBlock(folder, key)`, never read the whole file (up to 100 MB)
  or write heading text. `run_summary.BLOCKS` is the single place that defines blocks and their order, it feeds both
  `full_log.txt` and `run_summary.json`. Never change a heading's text: older files are found by it.
- Defaults reproduce the reference logs: formation 1.63 Å (C-S 1.9 Å, pairs with elements other than H/C/N/O use 1.1 x covalent radii), break 2.5 Å. Do not add element specific options,
  users raise `--break_distance` when a molecule needs it.
- Naming style is camelCase for functions/variables, docstrings in Google style with Args/Returns, as in the existing code.
- Files in the working tree may have CRLF line endings (`.gitattributes` normalises). Keep edits minimal so diffs stay small.

## Quick commands

```bash
conda activate aac            # has rdkit, networkx, graphviz, matplotlib, pytest
pytest                        # from this folder, uses src/ via pyproject (pythonpath)
python -m breakdown.main --help    # with PYTHONPATH=src, see testing skill
```

The `aac` environment still has the old `autoastrochem` package installed in editable mode from a **different checkout**
(`~/Documents/autoastrochem`) and does not contain `breakdown`, so use `PYTHONPATH=<this folder>/src` when running the CLI,
or `pip install -e .` from this folder. The project folder on disk is still called `autoastrochem`, the Python package is `breakdown`.
