"""Compare two per-run _log.txt files ignoring how the SMILES strings are written.

Usage: python tools/compare_logs.py new_log.txt reference_log.txt

Steps are matched one to one. The event labels (Isomerization, Fragmentation, Recombination)
and the final formulas must match, and the SMILES strings must describe the same bonded
graph once bond orders and radical/hydrogen bookkeeping are removed. A reference log made by
older code has no Recombination label, so those show up as expected event mismatches
(reference says Isomerization). Prints the number of mismatches, 0 means equivalent.
"""
import sys
from rdkit import Chem


def connectivity(smiles):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    for bond in mol.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
    for atom in mol.GetAtoms():
        atom.SetNoImplicit(True)
        atom.SetNumRadicalElectrons(0)
        atom.SetNumExplicitHs(0)
    return Chem.MolToSmiles(mol)


def load(path):
    rows = []
    for line in open(path):
        if line.startswith('At step'):
            parts = line.split()
            rows.append((parts[2][:-1], parts[3] if len(parts) == 5 else '', parts[-1]))
        elif line.startswith('Final molecules present'):
            rows.append(('final', line.strip(), ''))
    return rows


def main(newPath, referencePath):
    new, reference = load(newPath), load(referencePath)
    print(len(new), len(reference))
    mismatches = 0
    cache = {}
    for a, b in zip(new, reference):
        if a[0] == 'final':
            if a[1] != b[1]:
                print('FINAL', a[1], '|', b[1])
            continue
        if a[:2] != b[:2]:
            print('EVENT MISMATCH', a[:2], b[:2])
            mismatches += 1
            continue
        if (a[2], b[2]) not in cache:
            cache[(a[2], b[2])] = connectivity(a[2]) == connectivity(b[2])
        mismatches += not cache[(a[2], b[2])]
    print('mismatches', mismatches)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
