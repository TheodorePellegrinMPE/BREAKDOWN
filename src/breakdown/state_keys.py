"""Connectivity keys: names for molecule states that are the same in every folder.

The automatic state names (`PHE_A3`) number the states by how common they are in one folder, so the same structure
has different names at different energies. The key of a state only depends on which atoms are bonded to which, not
on the bond orders, the atom order in the trajectory or the folder, so it can be used to match states between folders
(and later DFT structures to MD states).
"""
import hashlib
from rdkit import Chem
from breakdown.auto_detect_utils import countAtoms, processSmileToFragment
from breakdown.run_summary import loadBlock

KEY_LENGTH = 12


def canonicalConnectivity(smiles):
    """Canonical SMILES of the bond graph only: every bond single, no hydrogens added, so bond orders and the
    hydrogens RDKit writes in brackets like [SH] do not matter.

    Args:
        smiles (str): SMILES string made by breakdown, may contain several fragments.

    Returns:
        str: Canonical SMILES of the connectivity.

    Raises:
        ValueError: Raised if the SMILES string cannot be read.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        raise ValueError(f'Invalid SMILES string: {smiles}')
    mol = Chem.RWMol(mol)
    for atom in mol.GetAtoms():
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)
        atom.SetFormalCharge(0)
        atom.SetNumRadicalElectrons(0)
        atom.SetIsAromatic(False)
    for bond in mol.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    mol.UpdatePropertyCache(strict=False)
    return Chem.MolToSmiles(mol, canonical=True)


def stateKey(smiles):
    """Short key of a state's connectivity.

    Args:
        smiles (str): SMILES string made by breakdown.

    Returns:
        str: First KEY_LENGTH hex digits of the SHA-1 of the canonical connectivity.
    """
    return hashlib.sha1(canonicalConnectivity(smiles).encode()).hexdigest()[:KEY_LENGTH]


def largestFragment(smiles):
    """The fragment with the most heavy atoms, then the most atoms, then the smallest key.

    Args:
        smiles (str): SMILES string, fragments separated by dots.

    Returns:
        str: SMILES string of the largest fragment.
    """
    def size(fragment):
        counts = countAtoms(fragment)
        heavy = sum(count for element, count in counts.items() if element != 'H')
        return -heavy, -sum(counts.values()), stateKey(fragment)
    return min(smiles.split('.'), key=size)


def makeStateKeys(isomerNames):
    """Keys of every state of a folder.

    Args:
        isomerNames (dict): SMILES string to state name, the isomerNames block.

    Returns:
        dict: State name to {key, largestFragmentKey, formula, largestFragmentFormula}. The formula lists the
        fragments joined by ' + ', largest first.
    """
    keys = {}
    for smiles, name in isomerNames.items():
        largest = largestFragment(smiles)
        fragments = sorted(smiles.split('.'), key=lambda f: sum(countAtoms(f).values()), reverse=True)
        keys[name] = {
            'key': stateKey(smiles),
            'largestFragmentKey': stateKey(largest),
            'formula': ' + '.join(processSmileToFragment(fragment) for fragment in fragments),
            'largestFragmentFormula': processSmileToFragment(largest),
        }
    return keys


def loadStateKeys(folder):
    """Keys of the states of a folder, from its stateKeys block, or computed from its names for older results.

    Args:
        folder (str): Folder with full_log.txt.

    Returns:
        dict or None: The result of makeStateKeys, None if the folder has no statistics.
    """
    keys = loadBlock(folder, 'stateKeys')
    if keys is not None:
        return keys
    isomerNames = loadBlock(folder, 'isomerNames')
    return makeStateKeys(isomerNames) if isomerNames is not None else None


def matchStates(folders, labels, byLargestFragment=False):
    """Line up the states of several folders by connectivity key.

    Args:
        folders (list): Folders with full_log.txt.
        labels (list): Short label of each folder, used as column names.
        byLargestFragment (bool): Match by the key of the largest fragment instead of the whole state, so
            for example C12H8 + H and C12H8 + H2 count as the same.

    Returns:
        list: One dictionary per key, with key, formula, and per label a list of (name, runs, steps) of the states
        with that key, sorted by the total number of runs, then steps, most common first.

    Raises:
        ValueError: Raised if a folder has no statistics.
    """
    keyField = 'largestFragmentKey' if byLargestFragment else 'key'
    formulaField = 'largestFragmentFormula' if byLargestFragment else 'formula'
    rows = {}
    for folder, label in zip(folders, labels):
        keys = loadStateKeys(folder)
        if keys is None:
            raise ValueError(f'{folder} has no full_log.txt with statistics, run makechartsandstats first.')
        runs = loadBlock(folder, 'isomerRuns') or {}
        steps = loadBlock(folder, 'isomerStats') or {}
        for name, info in keys.items():
            row = rows.setdefault(info[keyField], {'key': info[keyField], 'formula': info[formulaField],
                                                   'states': {l: [] for l in labels}})
            row['states'][label].append((name, runs.get(name, 0), steps.get(name, 0)))
    for row in rows.values():
        for states in row['states'].values():
            states.sort(key=lambda state: -state[1])
    def total(row, index):
        return sum(state[index] for states in row['states'].values() for state in states)
    return sorted(rows.values(), key=lambda row: (-total(row, 1), -total(row, 2), row['key']))
