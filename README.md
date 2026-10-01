# BREAKDOWN

**B**ond **R**earrangement, **E**vents **A**nd **K**inetics: **D**issociation **O**utput from **W**hole-trajectory a**N**alysis
(the N is a cheat, we are not sorry).

A tool for analysing molecular dynamics (MD) simulations of molecules that get hit with energy,
for example a PAH cation in space that absorbs a UV photon. During such a simulation the molecule
can shuffle its atoms around (**isomerization**), lose pieces (**fragmentation**) or take pieces back
(**recombination**). Reading that from thousands of coordinates by hand is impossible, this program does it for you:

1. It reads the atom coordinates of every simulation step and works out which atoms are bonded.
2. It writes the molecule at every step as a **SMILES string** (a line of text that describes the structure).
3. It collects statistics over many simulation runs: which structures appear, how often one turns into another,
   what the molecule breaks into, and how long it takes.
4. It draws plots and flow charts of all that.

This README is written for someone who has never used the program. Work through it top to bottom the first time.

Contents: [Background](#1-background-what-you-need-to-know) · [Install](#2-installation) ·
[Input files](#3-what-you-need-input-files) · [Quick start](#4-quick-start-tutorial) ·
[Understanding the output](#5-understanding-the-output) · [Command reference](#6-command-reference) ·
[Settings that matter](#7-settings-that-matter) · [Troubleshooting](#8-troubleshooting) ·
[Glossary](#9-glossary) · [For developers](#10-for-developers)

---

## 1. Background: what you need to know

**The simulation.** The molecule is simulated with the deMon-Nano program. The atoms move for a few tens of
picoseconds (ps). The program writes the position of every atom every few steps into a trajectory file
(`deMon.mol`). One such simulation is called a **run**. Because the start conditions are random, you always
do many runs of the same molecule (for example 100) and look at what happens on average.

**Steps and time.** The analysis works on the *output steps* in the trajectory file, not on the internal MD steps.
With the usual settings (0.5 fs time step, output every 10 steps) one output step is 5 fs, so 10 000 steps are 50 ps.
You can check this in the simulation output: `TIME STEP (FEMTOSECONDS)` and `STEPS BETWEEN OUTPUT`.
All plots use time in ps, the conversion is the option `--time_per_step_fs` (default 5).

**Bonds.** The simulation only gives positions, not bonds. The program decides that two atoms are bonded when they
come closer than a *formation distance* (1.63 Å, and 1.9 Å for carbon–sulfur; longer for heavier atoms such as Cl or Br)
and that they stay bonded until they are further apart than a *break distance* (2.5 Å by default). Two different distances
stop a bond from flickering on and off when atoms vibrate around one value. A hydrogen can only be bonded to one atom.

**SMILES.** A SMILES string is text that describes a molecule, for example `[H]C([H])([H])O[H]` is methanol
(every hydrogen is written out as `[H]`). Two steps have the same SMILES exactly when the same atoms are bonded, so a
change of the SMILES string means the molecule has changed. Fragments are separated by a dot:
`[H]C1C(=S)C([H])C([H])C([H])([H])C1[H].[H]CC[H]` is a ring molecule plus a separate C2H2 piece.
Only the connectivity is reliable: bond orders (single/double) are guessed, see [section 7](#7-settings-that-matter).

**Events.** When the SMILES differs from the step before, that step is labelled:

| Label | Meaning |
|---|---|
| `Isomerization` | same number of pieces, but the atoms are bonded differently (e.g. a hydrogen moved) |
| `Fragmentation` | one more piece than before |
| `Recombination` | one piece fewer than before (fragments joined again) |

---

## 2. Installation

You need Python 3.10 or newer, and [Graphviz](https://graphviz.org/download/) if you want the flow charts
(the `dot` command must work in your terminal). The Python packages are numpy, matplotlib, rdkit, networkx, graphviz and click,
they are installed automatically.

### Option A: conda (recommended if you use Anaconda/Miniconda)

```bash
conda create -n aac python=3.12
conda activate aac
conda install -c conda-forge graphviz        # only needed for flow charts
git clone https://github.com/TheodorePellegrinMPE/BREAKDOWN.git breakdown
cd breakdown      # the folder with pyproject.toml in it
pip install -e .
pip install pytest                            # only needed to run the tests
```

### Option B: Poetry

```bash
cd breakdown
poetry install
poetry run breakdown --help
```
With Poetry, put `poetry run` in front of every `breakdown` command below.

### Check that it works

```bash
breakdown --help              # lists the commands
pytest                        # should print "passed", takes a few seconds
```

`pip install -e .` means "editable": the program runs directly from this folder, so when you (or someone else) change the code the
changes are used immediately. If you have another copy of the project installed in the same environment, only one of them is used:
`python -c "import breakdown; print(breakdown.__file__)"` shows which.

Type command names in lowercase (`analysedemon`, not `analyseDemon`).

---

## 3. What you need: input files

Put all runs of **one molecule at one energy** in **one folder**, each run as a `.mol` file:

```
my_molecule_15ev/
    deMon_1.mol
    deMon_2.mol
    ...
    deMon_100.mol
```
(File names can be anything ending in `.mol`. A file with `keep` in its name is skipped.)

A `.mol` file is the trajectory written by deMon. Every step is a block: a line with the number of atoms, a line
starting with `BOMD INFORMATION`, and one line per atom with its element and x, y, z in Å (the rest of the line, velocities, is ignored):

```
13
 BOMD INFORMATION     -14.65821624      0.18377746      0.00000000 ...
  C     -3.401529    2.370002   -0.000483     -0.073450 ...
  C     -3.498345    0.981338   -0.010944      0.032666 ...
  ...
  S     -2.009500    4.749031   -0.003921 ...
```
The atoms must be in the same order in every step (deMon does this). The number of atoms and the elements are read from the file.
The `.out` and `.inp` files of a simulation are not needed for the analysis, only `.out` files are used by the
optional `makeenergiestable` command.

The folder `example_out_and_log_files/` contains one run for each of nine molecules (`*_traj7.mol`), the input file of
the simulation (`*.inp`), its output (`*.out`) and the log this program produces for it (`*_traj7_log.txt`). The `.mol`, `.out` and `full_log.txt` files are
very large (tens of MB or more), so a normal text editor may be slow, use `grep` or `less` to look for specific things.

---

## 4. Quick start tutorial

Always work in a copy, the program writes its results next to the input files.

```bash
mkdir ~/try_it
cp example_out_and_log_files/fluorene_25ev_traj7.mol ~/try_it/
```
In the commands below `~/try_it` is the **path of the folder** that holds the `.mol` files. Use your own path for your own data, for example
`~/simulations/phenol_5ev` or `/home/anna/thesis/fluorene/run1`. (`.` means "the folder I am in now", so after `cd ~/try_it` you could write `.` instead.)

### Step 1: from trajectory to a log of SMILES strings

```bash
breakdown analysedemonfolder ~/try_it
```
This reads every `.mol` file in the folder (in parallel) and writes a `..._log.txt` next to each. It takes a second or two per file.
To do a single file: `breakdown analysedemon ~/try_it fluorene_25ev_traj7.mol` (folder first, then the file name).

Open `~/try_it/fluorene_25ev_traj7_log.txt` in a text editor (a log is a few MB; the numbers below are from this file). It looks like this:

```
In file: deMon_7.mol
At step 0:               [H]C1C([H])C([H])C2C(C1[H])C1C([H])C([H])C([H])C([H])C1C2([H])[H]
At step 1:               [H]C1C([H])C([H])C2C(C1[H])C1C([H])C([H])C([H])C([H])C1C2([H])[H]
...
At step 79: Isomerization [H]C1C([H])C([H])C2C(C1[H])C([H])C1([H])C([H])C([H])C([H])C([H])C21
Number of 5-carbon rings: 1. Number of 6-carbon rings: 2. Tertiary carbon count: 1.
...
Final molecules present: C13H9, H
```
Each line is one step, followed by an event label if something changed. Directly after an event the program adds one line that lists the structure
of the new state (ring sizes and functional groups such as `Hydroxyl group count`, `Ketone count`, `Thiol group count`). The last line lists the
molecules present at the end as chemical formulas (the biggest first).

### Step 2: statistics and plots over all runs

```bash
breakdown makechartsandstats ~/try_it fluorene_25ev FLU
```
The three arguments are: the path of the folder with the `_log.txt` files (`~/try_it`), a name for the plot titles and files, and a short name for the molecule that is used to name its structures
(`FLU_A1`, `FLU_A2`, ...). The number of carbons and the elements are read from the logs.
Useful extras: `--draw` makes a picture of every structure (`drawings/`), `--coordinates` makes an `.xyz` file of every structure
(`coordinates/`), `-f ketone` follows one functional group of the starting molecule through the runs.

This creates in `~/try_it`: `full_log.txt` (all runs plus the statistics), `run_summary.json` (a small copy of the statistics that the
comparison commands read), and a `plots/` folder. With a single run the statistics only illustrate the format, in real use the folder holds many runs (about 100) of one molecule.

### Step 3 (optional): compare molecules, draw flow charts

See the [command reference](#6-command-reference).

---

## 5. Understanding the output

### Files

| File | What it is |
|---|---|
| `<run>_log.txt` | The step-by-step SMILES of one run, from `analysedemon`. |
| `full_log.txt` | All run logs joined, followed by the statistics blocks (below). It can be 100 MB or more (millions of lines), so a normal text editor may be slow. To look for something specific use `grep` (e.g. `grep -n "Fragmentation" full_log.txt`), a pager such as `less` (press `G` to jump to the statistics at the end), or a lightweight editor that handles large files. |
| `run_summary.json` | The statistics blocks and a few per-step counts in a small file, so that commands that compare molecules do not have to read `full_log.txt`. Delete it if you edit `full_log.txt` by hand, it is then ignored automatically. |
| `plots/` | Figures, see below. |
| `drawings/`, `coordinates/` | Pictures and `.xyz` files of the structures (only with `--draw` / `--coordinates`). |
| `transitionGraphs/` | Flow charts from `makeflowcharts`. |

### The statistics blocks at the end of `full_log.txt`

Each block starts with a heading line, followed by a JSON dictionary. Structures get names like `FLU_A3`: the letter is
how many carbons were lost (`A` = none, `B` = one, ...), the number ranks the structures by how long they were seen.
The block `Automatically generated molecule names.` translates names back to SMILES.

| Heading | Content |
|---|---|
| Total number of steps at which an isomer is present | Steps spent in each structure, summed over all runs. |
| Total number of runs in which an isomer appeared | In how many runs each structure was seen. |
| Number of runs at which each chemical fragment was present at the end | Final fragments (formulas) and in how many runs. |
| Total / Net number of transitions to an isomer from an isomer | `{"to": {"from": count}}`. "Net" subtracts the reverse transitions. |
| Number of transitions only including fragmentation | Transitions where the number of pieces changes. |
| Chemical names of fragments / Molecular masses | Formulas and masses (whole numbers) of the structures. |
| Runs containing original functional group ... at step | Only with `-f`: for each step, in how many runs the group is still there. |
| Runs with each oxygen/sulfur in each final group and/or fragment | Only with `-f`: where the heteroatom ended up. `other` means none of the tracked groups. |
| First step of each event type in each run | For every run, the first step of an isomerization, fragmentation, recombination and of any change, and the last step. `null` means it never happened. |

### The plots in `plots/`

Every figure is saved twice: `.pdf` (vector, for manuscripts) and `.png` (300 dpi, for slides).

| File | How to read it |
|---|---|
| `<name>_frags_isos_pie` | Nested pie: outer ring = runs with / without an isomerization, inner ring = runs that were / were not in more than one piece *at the end*. A run that lost something for a moment and took it back is not counted as fragmented. |
| `<name>_fragmentations_pie` | What the runs ended as: no fragmentation, or the smaller piece that came off (`H`, `C2H2`, classes such as `CmHnOk`, or multiple fragments). |
| `<name>_fragmentations_bar` | Every molecule present at the end of a run and in how many runs, the main molecule included. |
| `<name>_dissociation` | Percentage of runs (y) in which the largest fragment has a given formula, against time. The legend shows the formula as C<sub>m</sub>H<sub>n</sub>O<sub>k</sub> with ranges, e.g. C<sub>6</sub>H<sub>4-6</sub>O<sub>0-1</sub>: 6 carbons, 4 to 6 hydrogens and 0 to 1 oxygens were seen in the largest fragment of that size. |
| `<name>_first_event_survival` | **Survival curves.** y = percentage of runs in which nothing has happened yet, x = time. One curve for "any change", one for isomerization and one for fragmentation, with the median time (when half of the runs had it, n.r. if fewer than half did) and the number of runs with the event. |

Read the survival plot like this: a curve that falls quickly means the molecule is unstable at this energy, a curve that stays at 100 % means
the event never happened. To claim that molecule A is less stable than molecule B, their 95 % bands (drawn when a plot has two curves or fewer) should not overlap.
The method (Kaplan–Meier) counts runs where the event did not happen as "still waiting until the end of the run", not as failed runs.

### Log in a nutshell

```
step 0 ────── same structure ─────────────── step 147: Isomerization ──── ... ─── step 5652: Fragmentation ─── ... ─── end
```
A useful way to work is: find interesting events in `full_log.txt`'s statistics, look up the structure name in `drawings/`, and go to the
step in the run log to see how it happened.

### Flow charts (`makeflowcharts`)

A flow chart shows the route the molecule takes from its starting structure to another structure, as boxes (structures) joined by arrows (transitions).
It is built from the statistics blocks in `full_log.txt`, summed over **all** runs:

1. Every structure that was seen is a possible box. Every transition from structure A to structure B that happened in any run, and how often, is a possible arrow.
2. For each structure you want to reach (the *end states*) the program looks for the best path from the start structure. "Best" means the path whose weakest arrow is as strong as
   possible: the path is as strong as its rarest step, so it prefers routes that many runs actually took over rare shortcuts. The strength of an arrow is its *net* count, transitions A to B minus transitions B to A,
   so a pair of structures that just flips back and forth is not a strong route. A path that has to go against the net flow somewhere gets the worst possible score and is only used when there is no other path. Ties are broken by the average strength.
   The search first looks up to 5 steps deep and goes deeper (up to 20) only if it finds no path (`-mb`, `-mp`).
3. The end states are the 10 structures the molecule spent the most steps in (`-n`), up to 5 structures where the molecule has broken into several carbon containing pieces and that were seen
   for more than 30 steps in total (`-ns`), and any structure you name with `-s`. One chart is made for each end state (`transitions_graph_<state>.pdf`) and one with all the paths together (`full_transitions_graph.pdf`).
   `--direct` instead draws only the start structure and the structures it turned into directly (`direct_transitions_<start>.pdf`).

What the numbers mean:

- **Number in a box** (below the structure name): the total number of output steps the molecule was in that structure, added up over all runs. Multiply by the time per step
  (5 fs by default) to get a time. It is not a number of runs, one run can contribute thousands of steps.
- **Number on an arrow**: how many times that transition happened, added up over all runs. In the path charts this is the total count of the transition in that direction
  (when the reverse transition also happened it gets its own arrow with its own count). In the `--direct` chart it is the *net* count (A to B minus B to A).
- A transition that many runs took is not the same as a likely one: a rare structure can still be reached through many arrows if the runs spent a long time near it.

If you made the statistics with `--draw`, each box contains a picture of the structure. `makeflowcharts` also writes `visited_states.txt` next to `full_log.txt`, the names of all
structures in the charts, which `optimizeisomer` uses.

---

## 6. Command reference

Get the options of any command with `breakdown <command> --help`.

`FOLDER` is always the **path of a folder**, for example `~/try_it` or `/home/anna/thesis/phenol_5ev`, or `.` for the folder you are currently in.
Examples: `breakdown makeflowcharts ~/try_it FLU_A1`, `breakdown makechartsandstats ~/simulations/phenol_5ev phenol_5ev PHE --draw`.

| Command | What it does |
|---|---|
| `analysedemon FOLDER FILE.mol` | One `.mol` file to `_log.txt`. |
| `analysedemonfolder FOLDER` | All `.mol` files in a folder, `-j 8` for 8 at a time, `--skip_existing` to skip finished ones, `-p "run*.mol"` to select files. |
| `makechartsandstats FOLDER NAME MOLNAME` | Joins the logs in a folder, writes `full_log.txt`, `run_summary.json` and the plots. Options: `--draw`, `--coordinates`, `-f GROUP` (`ketone`, `hydroxyl`, `thiol`, ..., add `_2` to require at least two, e.g. `hydroxyl_2`), `-s` number of carbons, `-e C,H,O` elements, `-t` time per step in fs, `-lc 14,12` carbon counts to show in the dissociation legend. |
| `makeflowcharts FOLDER START_STATE` | Flow charts (PDF) of how the molecule gets from the start structure (e.g. `FLU_A1`) to the most common structures, from `full_log.txt`. `--direct` only draws the structures reached directly from the start. Needs Graphviz. |
| `comparefunctionalgroups -d DIR -f "mol1;mol2" -n "A;B" -g "ketone;hydroxyl"` | Plots how fast the tracked functional group disappears in each molecule (`T50`: time when half of the runs lost it). The folders must have been analysed with `-f`. |
| `compareelementalstability -d DIR -f "mol1;mol2" -n "A;B"` | Compares how long the carbon skeleton and the oxygen (`--heteroatom S` for sulfur) stay together. |
| `comparefirstevents -d DIR -f "mol1;mol2" -n "A;B" -e fragmentation` | Survival curves of several molecules in one plot. `-e` is `anyChange`, `isomerization`, `fragmentation` or `recombination`. |
| `plotoxygenfates -d DIR -f "mol1;mol2" -n "A;B" -g "ketone;hydroxyl"` | Pie charts of where the oxygen went in runs that lost the starting group. |
| `makeenergiestable FOLDER` | Table of the DFTB electronic energies in the `.out` files of *optimizations* (not of MD runs). **Do not use these energies for results, see the warning below.** |
| `optimizeisomer FOLDER NAME` | Optimizes the geometry of a structure with deMon-Nano (DFTB), only for structures in `visited_states.txt` (needs the deMon executable `deMon.static-website.x` in the folder, and the structure as `NAME.xyz`). The input it writes assumes charge +1 and multiplicity 2. **See the warning below.** |

> **Warning: DFTB is not accurate enough for energies and optimized geometries.** The simulations (and `optimizeisomer` and `makeenergiestable`, which use the same
> method) are done with DFTB, a fast approximate method that is fine for following where atoms go during hundreds of thousands of steps. It is **not** good enough for
> reliable relative energies of isomers, reaction energies, barriers or final optimized geometries. Use these two commands only to get a first geometry, and to
> look at trends. For numbers you want to report, redo the geometry optimization and energies with DFT in an external program such as Gaussian or ORCA, starting from the
> `.xyz` files that `--coordinates` writes (those are only rough force-field geometries, so they have to be optimized anyway).

For the `compare…` commands `-d` is a base directory, `-f` lists folders inside it separated by `;`, and the result goes to `DIR/comparisons/`.

---

## 7. Settings that matter

| Option | Default | Meaning |
|---|---|---|
| `--formation_distance`, `-fd` | 1.63 Å | Two atoms closer than this become bonded. Carbon–sulfur uses 1.9 Å, pairs with heavier atoms (Cl, Br, S–S, ...) use 1.1 × the sum of their covalent radii when that is larger. |
| `--break_distance`, `-bd` | 2.5 Å | A bond breaks when the atoms are further apart than this. **This is the most important setting.** Too small and vibrating bonds look like fragmentation, too large and real dissociation is noticed late. The default reproduces the logs in `example_out_and_log_files/`. |
| `--min_lifetime`, `-ml` | 1 | A new structure must last at least this many steps to count. `-ml 5` removes structures that flicker for a few steps. Use it only when you understand why you need it and report it. |
| `--time_per_step_fs`, `-t` | 5 | Time between output steps in fs. |

How to choose `--break_distance`: run one trajectory with a few values (2.5, 3.0, 3.5), and compare the number of events (`grep -c Isomerization`).
If the count changes a lot, the result depends on your choice and you should say so in your thesis.

**About bond orders.** The program does not know if a bond is single or double. In the SMILES, only a terminal O or S atom is drawn as double bonded (`C=O`),
all other bonds are single. So `[H]C1C([H])...` is a description of *which atoms are bonded to which*. To turn it into a real molecule (drawing, 3D coordinates) the program
assigns double bonds so that the valences fit, and marks what is left as radicals. **Never read these SMILES strings with RDKit's `Chem.MolFromSmiles`**, it adds
hydrogens to every atom with a missing bond. Use `breakdown.molecule_utils.moleculeFromSmiles`. Two structures that differ only in bond order count as the same structure.

---

## 8. Troubleshooting

| Message / problem | Cause and fix |
|---|---|
| `invalid file: ... (expected an existing .mol file)` | The file does not exist or does not end in `.mol`. |
| `... has 13 atoms per step, but the molecule size is 12` | The `-s` you gave is not the number of atoms in the file. Leave `-s` out, it is read from the file. |
| `... ends with an incomplete step` | The simulation stopped in the middle of writing a step. The last, incomplete step is ignored, the rest is analysed. |
| `No _log.txt files found` | Run `analysedemon` / `analysedemonfolder` first. |
| `... is incomplete and is ignored` | A log has no `Final molecules present` line (crashed analysis). Delete it and re-run `analysedemon` for that file. |
| `... has no first event data, run makechartsandstats again` | The folder was analysed with an older version. Run `makechartsandstats` again. |
| `... is not a state name` | The start state given to `makeflowcharts` is not in `Automatically generated molecule names.` of `full_log.txt`. |
| Flow charts fail with `dot not found` | Install Graphviz and make sure `dot -V` works. |
| Results are different from a colleague's | Check `--break_distance`, `--formation_distance` and `--min_lifetime`, and that both use the same version of the program. |
| `breakdown: command not found` | The environment is not activated (`conda activate aac`), or with Poetry use `poetry run breakdown`. |
| Plots look empty or have only one line | Too few runs, or nothing happened in the simulation (at low energy most runs never change). |

---

## 9. Glossary

- **MD (molecular dynamics)**: simulating how atoms move by calculating forces at every time step.
- **Trajectory**: the positions of all atoms at every output step. **Run**: one simulation.
- **SMILES**: text notation for a molecule. **Canonical**: written in a unique order, so equal molecules give equal text.
- **Isomer**: same atoms, different bonding. **Fragmentation**: the molecule splits into pieces.
- **Radical**: an atom with an unpaired electron (a bond "missing"). Very common in these simulations.
- **Å (ångström)**: 0.1 nm. Bonds between C, H, O are about 1–1.5 Å. **fs / ps**: femto- / picosecond, 10⁻¹⁵ s / 10⁻¹² s.
- **T50**: time at which half of the runs have lost something (a functional group, or intactness).
- **Kaplan–Meier / survival curve**: estimate of the fraction of runs without an event yet, handling runs that ended without one.

---

## 10. For developers

```bash
pip install -e . && pip install pytest
pytest                                  # unit tests plus synthetic molecules (about 50 tests, a few seconds)
```

- `CLAUDE.md` and `.agent/skills/` explain the structure of the code, how to test it against the example data, how the SMILES/bond-order design works
  and which files to touch when adding a functional group, statistic, plot or command. They are written for AI coding agents, but are a good guide for people too.
- `tools/compare_logs.py new_log.txt reference_log.txt` compares two run logs by event and by bonded graph.
- Source layout: `src/breakdown/` (one module per job, `main.py` holds only the command line), tests in `tests/`.
