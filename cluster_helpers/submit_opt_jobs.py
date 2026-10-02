"""Optimize the geometry of each molecule with deMon-Nano (DFTB) on this machine, starting from LOCAL_BASE/<molecule>/<molecule>.xyz.

The result, <molecule>.out, is the start geometry of the MD runs. DFTB geometries are fine as a starting point of an MD
simulation, but do not use them or their energies as results: optimize and compute energies with DFT in Gaussian or ORCA.
"""
import os
import shutil
import subprocess
import sys

import cluster_config as cfg
from cluster_common import localDir, makeParser, selectedJobs

SUCCESS_STRING = 'OPTIMIZED STRUCTURE IN ANGSTROM (INPUT ORDER)'


def createOptInput(mol, molDir):
    """Read {mol}.xyz and write {mol}_opt.inp with the deMon settings."""
    xyzPath = os.path.join(molDir, f'{mol}.xyz')
    inpPath = os.path.join(molDir, f'{mol}_opt.inp')
    if not os.path.exists(xyzPath):
        raise FileNotFoundError(f'Source XYZ file not found: {xyzPath}')
    with open(xyzPath) as file:
        geometry = ''.join(file.readlines()[2:])  # the first two lines are the atom count and a comment
    with open(inpPath, 'w') as file:
        file.write(f"""# {mol} optimisation
OPTIMIZATION
DFTB SCC FERMI={cfg.OPT_FERMI_K} DISP=2
PARAM PTYPE=BIO
{cfg.LOCAL_BASIS}
MULTIPLICITY {cfg.MULTIPLICITY}
CHARGE {cfg.CHARGE}
GEOMETRY
{geometry}
""")
    return inpPath


def optimize(mol, timeoutHours):
    molDir = localDir(mol)
    if not os.path.exists(molDir):
        print(f'Directory not found: {molDir}')
        return False
    print(f'--- Processing: {mol} ---')
    inpSource = createOptInput(mol, molDir)
    inpDest = os.path.join(molDir, 'deMon.inp')
    execDest = os.path.join(molDir, 'deMon.static-website.x')
    outFile = os.path.join(molDir, 'deMon.out')
    molFile = os.path.join(molDir, 'deMon.mol')

    # An output of an earlier run must not be mistaken for the result of this one
    for stale in (outFile, molFile):
        if os.path.exists(stale):
            os.remove(stale)
    shutil.copy2(inpSource, inpDest)
    shutil.copy2(cfg.LOCAL_EXECUTABLE, execDest)

    print(f'  Running deMon for {mol}...')
    try:
        result = subprocess.run(['./deMon.static-website.x', 'deMon.inp'], cwd=molDir, timeout=timeoutHours * 3600)
    except subprocess.TimeoutExpired:
        print(f'  Error: deMon did not finish within {timeoutHours} hours, stopped')
        return False
    if result.returncode != 0:
        print(f'  Warning: deMon exited with status {result.returncode}')

    if not os.path.exists(outFile):
        print(f'  Error: deMon.out was not generated for {mol}.')
        return False
    with open(outFile) as file:
        converged = SUCCESS_STRING in file.read()
    if not converged:
        print(f'  Warning: the optimization of {mol} did not finish, see {outFile}')
        return False

    print(f'  Success! Renaming and cleaning up {mol}...')
    os.rename(outFile, os.path.join(molDir, f'{mol}.out'))
    if os.path.exists(molFile):
        os.rename(molFile, os.path.join(molDir, f'{mol}.mol'))
    for temporary in (inpDest, execDest, os.path.join(molDir, 'deMon.keep.mol')):
        if os.path.exists(temporary):
            os.remove(temporary)
    return True


def main():
    parser = makeParser('Optimize geometries locally.', energies=False)
    parser.add_argument('--timeout-hours', type=float, default=6, help='stop a deMon run after this long')
    args = parser.parse_args()
    jobs = selectedJobs(args)
    molecules = list(dict.fromkeys(mol for mol, _, _ in jobs))
    if args.dry_run:
        for mol in molecules:
            print(f'  [dry run] would optimize {mol} from {localDir(mol)}/{mol}.xyz')
        return
    results = {mol: optimize(mol, args.timeout_hours) for mol in molecules}
    failed = [mol for mol, ok in results.items() if not ok]
    print('\nPipeline finished.' + (f' Failed: {", ".join(failed)}' if failed else ''))
    if failed:
        sys.exit(1)


if __name__ == '__main__':
    main()
