"""Make the deMon-Nano MD input and the Slurm array script of every molecule and energy, in LOCAL_BASE/<molecule>/.

The start geometry is read from <molecule>.out (the optimized structure). Nothing is sent to the cluster here,
queue_md_jobs.py does that. Also writes md_settings_<energy>ev.json, the record of the settings the runs are made with.
"""
import json
import os
from datetime import datetime, timezone

import cluster_config as cfg
from cluster_common import (inpName, localDir, makeParser, runDirName, selectedJobs, timePerStepFs, totalPs,
                            energyText, nodeLine)

KB = 1.380649e-23  # J/K
EV_TO_J = 1.6021766e-19
GEOMETRY_MARKER = 'OPTIMIZED STRUCTURE IN ANGSTROM (INPUT ORDER)'


def readOptimizedGeometry(outPath):
    """Heavy atoms first, then hydrogens, as formatted geometry lines, from the optimization output."""
    with open(outPath) as file:
        lines = file.readlines()
    starts = [i for i, line in enumerate(lines) if GEOMETRY_MARKER in line]
    if not starts:
        return None
    heavy, hydrogens = [], []
    for line in lines[starts[-1] + 4:]:
        if line.strip() == '':
            break
        parts = line.split()
        if len(parts) >= 6:
            formatted = f'{parts[1]:<10} {parts[3]:>12} {parts[4]:>12} {parts[5]:>12}'
            (hydrogens if parts[1].upper() == 'H' else heavy).append(formatted)
    return heavy + hydrogens


def initialTemperature(numAtoms, energyEv):
    """Temperature at which 3N degrees of freedom hold the kinetic energy: E = N_dof * kB * T / 2."""
    return 2 * energyEv * EV_TO_J / (3 * numAtoms * KB)


def mdInput(mol, geometry, temperature):
    return f"""# {mol} md test
MDYNAMICS RANDOM={temperature:.0f}
TIMESTEP {cfg.MD_TIMESTEP_FS}
MDSTEPS MAX={cfg.MD_STEPS}
DFTB SCC FERMI={cfg.MD_FERMI_K} DISP=2
PARAM PTYPE=BIO
{cfg.REMOTE_BASE}/basis/
MULTIPLICITY {cfg.MULTIPLICITY}
CHARGE {cfg.CHARGE}
GEOMETRY
""" + '\n'.join(geometry)


def slurmScript(mol, abbr, energy, inpFile):
    """One array task is one trajectory. A finished run leaves a marker file and is skipped when resubmitted."""
    ev = energyText(energy)
    runDir = f'{cfg.REMOTE_BASE}/{mol}/{runDirName(energy)}'
    scratch = f'${{TMPDIR:-{cfg.REMOTE_BASE}/tmp}}'
    return f"""#!/bin/bash
#SBATCH --job-name={abbr}_{ev}ev
{nodeLine()}#SBATCH --nodes 1
#SBATCH --time {cfg.MD_TIME}
#SBATCH --mem {cfg.MD_MEMORY}
#SBATCH --array=1-{cfg.N_RUNS}%{cfg.MAX_PARALLEL_RUNS}
#SBATCH --output={cfg.REMOTE_BASE}/error_files/deMonNano_%A_%a.out
#SBATCH --error={cfg.REMOTE_BASE}/error_files/deMonNano_%A_%a.err

set -euo pipefail
module purge

RUN_ID=$SLURM_ARRAY_TASK_ID
RUN_DIR="{runDir}"
mkdir -p "$RUN_DIR" "{cfg.REMOTE_BASE}/error_files" "{scratch}"

if [ -f "$RUN_DIR/.done_$RUN_ID" ]; then
    echo "Run $RUN_ID is already finished, skipping"
    exit 0
fi

echo "### Starting run $RUN_ID on $(hostname) ###"
WORK_DIR=$(mktemp -d "{scratch}/deMon.XXXXXX")
trap 'rm -rf "$WORK_DIR"' EXIT

cp "{cfg.REMOTE_BASE}/{mol}/{inpFile}" "$WORK_DIR/deMon.inp"
cp "{cfg.REMOTE_BASE}/deMon.static-website.x" "$WORK_DIR/"
cd "$WORK_DIR"

# deMon-Nano is not known to return 0 on success, so the result is judged by its output files
./deMon.static-website.x deMon.inp || echo "WARNING: deMon exited with status $?" >&2
if [ ! -s deMon.out ] || [ ! -s deMon.mol ]; then
    echo "ERROR: run $RUN_ID produced no deMon.out / deMon.mol" >&2
    exit 1
fi

mv deMon.out "$RUN_DIR/deMon_$RUN_ID.out"
mv deMon.mol "$RUN_DIR/deMon_$RUN_ID.mol"
[ -f deMon.keep.mol ] && mv deMon.keep.mol "$RUN_DIR/deMon.keep_$RUN_ID.mol"
[ -f deMon.qmd ] && cp deMon.qmd "$RUN_DIR/deMon.qmd_$RUN_ID"
touch "$RUN_DIR/.done_$RUN_ID"
echo "### Finished run $RUN_ID ###"
"""


def main():
    args = makeParser('Make MD inputs and Slurm scripts.').parse_args()
    for mol, abbr, energy in selectedJobs(args):
        outPath = os.path.join(localDir(mol), f'{mol}.out')
        if not os.path.exists(outPath):
            print(f'{mol}: {outPath} not found, skipped (optimize the molecule first, see submit_opt_jobs.py)')
            continue
        geometry = readOptimizedGeometry(outPath)
        if geometry is None:
            print(f'{mol}: no "{GEOMETRY_MARKER}" in {outPath}, was the optimization finished?')
            continue
        temperature = initialTemperature(len(geometry), energy)
        print(f'{mol}, {energyText(energy)} eV: {len(geometry)} atoms, initial temperature {temperature:.0f} K')

        inpFile = inpName(mol, energy)
        slurmFile = f'{abbr}_{energyText(energy)}ev.slurm'
        settings = {
            'molecule': mol, 'abbreviation': abbr, 'kinetic_energy_ev': energy, 'atoms': len(geometry),
            'initial_temperature_k': round(temperature, 1), 'timestep_fs': cfg.MD_TIMESTEP_FS, 'md_steps': cfg.MD_STEPS,
            'output_interval': cfg.MD_OUTPUT_INTERVAL, 'time_per_output_step_fs': timePerStepFs(),
            'simulated_ps': totalPs(), 'charge': cfg.CHARGE, 'multiplicity': cfg.MULTIPLICITY,
            'fermi_k': cfg.MD_FERMI_K, 'runs': cfg.N_RUNS, 'input': inpFile, 'slurm': slurmFile,
            'created': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        }
        if args.dry_run:
            print(f'  [dry run] would write {inpFile}, {slurmFile} and the settings file in {localDir(mol)}')
            continue
        with open(os.path.join(localDir(mol), inpFile), 'w') as file:
            file.write(mdInput(mol, geometry, temperature))
        with open(os.path.join(localDir(mol), slurmFile), 'w') as file:
            file.write(slurmScript(mol, abbr, energy, inpFile))
        with open(os.path.join(localDir(mol), f'md_settings_{energyText(energy)}ev.json'), 'w') as file:
            json.dump(settings, file, indent=2)

    print(f'\nDone. Each Slurm script is an array of {cfg.N_RUNS} single trajectories, at most '
          f'{cfg.MAX_PARALLEL_RUNS} at the same time.')


if __name__ == '__main__':
    main()
