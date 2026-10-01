import networkx as nx
from rdkit import Chem

# Number of bonds an atom can have. The lowest listed value that is at least the
# number of bonds it has is taken as its valence, any bonding missing to reach it
# is either made up by a double bond to a neighbour with the same problem, or it
# is a radical electron.
VALENCES = {'H': (1,), 'C': (4,), 'N': (3, 5), 'O': (2,), 'S': (2, 4, 6), 'P': (3, 5),
            'F': (1,), 'Cl': (1,), 'Br': (1,), 'I': (1,)}
DIVALENT_ATOMS = ('O', 'S')
BOND_TYPES = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}


def moleculeFromSmiles(smiles):
    """Make a sanitized RDKit molecule from a SMILES string made by breakdown.

    The SMILES strings from the simulation only reliably give which atoms are
    bonded, and all hydrogens are explicit atoms. RDKit would add hydrogens to atoms
    with too few bonds when the string is read, which is wrong, these atoms are radicals.
    So here no implicit hydrogens are allowed. To get reasonable bond orders, an
    oxygen or sulfur atom with more than one neighbour has only single bonds, then
    neighbouring atoms that both have unused valence are joined by double (or triple)
    bonds, and whatever valence is left is a radical electron.

    Args:
        smiles (str): SMILES string of a single molecule or of several fragments.

    Returns:
        RDKit Mol: Sanitized molecule with explicit hydrogens only.

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
    makeDivalentAtomBondsSingle(mol)
    joinAtomsWithFreeValence(mol)
    for atom in mol.GetAtoms():
        atom.SetNumRadicalElectrons(getFreeValence(atom))
    mol = mol.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def getBondOrderSum(atom):
    """Sum of the orders of all bonds of an atom."""
    return int(round(sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())))


def getFreeValence(atom):
    """Number of bonds an atom could still make without being over its valence.

    Args:
        atom (RDKit Atom): Atom in a molecule.

    Returns:
        int: Free valence, 0 if the atom is not a known element or is full.
    """
    used = getBondOrderSum(atom)
    valences = VALENCES.get(atom.GetSymbol())
    if valences is None:
        return 0
    valence = min((v for v in valences if v >= used), default=used)
    return valence - used


def makeDivalentAtomBondsSingle(mol):
    """Set all bonds of oxygen and sulfur atoms that have several neighbours to single.

    Args:
        mol (RDKit RWMol): Molecule that will be edited.
    """
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in DIVALENT_ATOMS and atom.GetDegree() > 1:
            for bond in atom.GetBonds():
                bond.SetBondType(Chem.BondType.SINGLE)


def joinAtomsWithFreeValence(mol):
    """Raise the order of bonds between neighbouring atoms that both have free valence.

    The bonds are chosen as a maximum matching, so as many atoms as possible use up
    their free valence. Repeated so that triple bonds can form.

    Args:
        mol (RDKit RWMol): Molecule that will be edited.
    """
    for _ in range(2):
        graph = nx.Graph()
        for bond in mol.GetBonds():
            atomA, atomB = bond.GetBeginAtom(), bond.GetEndAtom()
            if getFreeValence(atomA) > 0 and getFreeValence(atomB) > 0:
                graph.add_edge(atomA.GetIdx(), atomB.GetIdx())
        matching = nx.max_weight_matching(graph, maxcardinality=True)
        if not matching:
            return
        for idxA, idxB in matching:
            bond = mol.GetBondBetweenAtoms(idxA, idxB)
            newOrder = min(int(bond.GetBondTypeAsDouble()) + 1, 3)
            bond.SetBondType(BOND_TYPES[newOrder])
