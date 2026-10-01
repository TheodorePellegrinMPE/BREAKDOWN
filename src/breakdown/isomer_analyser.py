from functools import partial
from rdkit import Chem

def analyseIsomer(smilesString):
    """Creates a RDKit molecule from a SMILES string and then looks for
    different atomic groups within it, returning a readable string
    reporting the counts of present groups.

    Args:
        smilesString (str): A string in the SMILES format representing a molecule.

    Returns:
        str: A string giving information about the chemical groups present.
    """
    mol = Chem.MolFromSmiles(smilesString, sanitize=False)
    if not mol:
        return ""
        
    isomerString = ''
    isomerString += identifyAllRings(mol)
    
    # List of functions and their display names
    functions = [
        (identifyHydroxylGroups, "Hydroxyl group"),
        (identifyMethyleneGroups, "Methylene group"),
        (identifyTertiaryCarbon, "Tertiary carbon"),
        (identifyKetone, "Ketone"),
        (identifyKetene, "Ketene"),
        (identifyEther, "Ether"),
        (identifyAldehyde, "Aldehyde"),
        (identifyCarboxyl, "Carboxyl"),
        (identifyEpoxide, "Epoxide"),
        (partial(identifyHydroxylGroups, element='S'), "Thiol group"),
        (partial(identifyKetone, element='S'), "Thioketone"),
        (partial(identifyAldehyde, element='S'), "Thioaldehyde"),
        (partial(identifyEther, element='S'), "Thioether"),
        (partial(identifyEpoxide, element='S'), "Thiirane"),
        (identifyDisulfide, "Disulfide"),
    ]
    
    for func, label in functions:
        count = func(mol)
        if count > 0:
            isomerString += f"{label} count: {count}. "
            
    return isomerString

def identifyAllRings(mol):
    """Gets a dictionary for ring lengths and reporting counts."""
    ringLengths = findRingLengths(mol)
    ringsString = ''
    for ringKey in sorted(ringLengths.keys()):
        ringsString += f'Number of {ringKey}-carbon rings: {ringLengths[ringKey]}. '
    return ringsString

def findRingLengths(mol):
    """Generates a dictionary of ring lengths and their counts."""
    ringInfo = Chem.GetSymmSSSR(mol)
    ringLengths = {}
    for ring in ringInfo:
        ringLength = len(ring)
        ringLengths[ringLength] = ringLengths.get(ringLength, 0) + 1
    return ringLengths

def identifyHydroxylGroups(mol, element='O'):
    """Counts oxygens attached to one H and one C, excluding carboxylic acids via connectivity."""
    count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == element:
            neighbors = atom.GetNeighbors()

            if len(neighbors) != 2:
                continue

            neighbor_symbols = [n.GetSymbol() for n in neighbors]

            if 'H' in neighbor_symbols and 'C' in neighbor_symbols:
                carbon = [n for n in neighbors if n.GetSymbol() == 'C'][0]
                is_carboxyl = False
                
                for c_neighbor in carbon.GetNeighbors():
                    if (c_neighbor.GetSymbol() == element and 
                        c_neighbor.GetIdx() != atom.GetIdx()):
                        if len(c_neighbor.GetNeighbors()) == 1:
                            is_carboxyl = True
                            break
                
                if not is_carboxyl:
                    count += 1
    return count

def identifyMethyleneGroups(mol):
    """Counts carbons attached to exactly 2 hydrogens."""
    count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            h_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'H']
            if len(h_neighbors) == 2:
                count += 1
    return count

def identifyTertiaryCarbon(mol):
    """Counts carbons attached to 3 carbons and 1 hydrogen."""
    count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            c_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'C']
            h_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'H']
            if len(c_neighbors) >= 3 and len(h_neighbors) >= 1:
                count += 1
    return count

def identifyKetone(mol, element='O'):
    """Counts carbons bonded to 2 carbons and 1 terminal oxygen."""
    count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            neighbors = atom.GetNeighbors()
            c_neighbors = [n for n in neighbors if n.GetSymbol() == 'C']
            o_neighbors = [n for n in neighbors if n.GetSymbol() == element]
            if len(c_neighbors) == 2 and len(o_neighbors) == 1:
                if len(o_neighbors[0].GetNeighbors()) == 1:
                    count += 1
    return count

def identifyKetene(mol):
    """Counts ketene functional groups."""
    count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            c_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'C']
            if len(c_neighbors) == 3:
                for cn in c_neighbors:
                    cn_neighbors = cn.GetNeighbors()
                    if len(cn_neighbors) == 2:
                        symbols = [n.GetSymbol() for n in cn_neighbors]
                        if 'C' in symbols and 'O' in symbols:
                            oxy = [n for n in cn_neighbors if n.GetSymbol() == 'O'][0]
                            if len(oxy.GetNeighbors()) == 1:
                                count += 1
    return count

def identifyEther(mol, element='O'):
    """Counts oxygens bonded to 2 non-adjacent carbons."""
    count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == element:
            c_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'C']
            if len(c_neighbors) == 2 and len(atom.GetNeighbors()) == 2:
                if not mol.GetBondBetweenAtoms(c_neighbors[0].GetIdx(), c_neighbors[1].GetIdx()):
                    count += 1
    return count

def identifyAldehyde(mol, element='O'):
    """Counts carbons bonded to 1 carbon, 1 terminal oxygen, and 1 hydrogen."""
    count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            c_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'C']
            o_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == element]
            h_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'H']
            if len(c_neighbors) == 1 and len(o_neighbors) == 1 and len(h_neighbors) == 1:
                if len(o_neighbors[0].GetNeighbors()) == 1:
                    count += 1
    return count

def identifyCarboxyl(mol):
    """Counts carboxyl groups."""
    count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            c_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'C']
            o_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'O']
            if len(c_neighbors) == 1 and len(o_neighbors) == 2:
                has_term = any(len(o.GetNeighbors()) == 1 for o in o_neighbors)
                has_oh = any(any(n.GetSymbol() == 'H' for n in o.GetNeighbors()) for o in o_neighbors)
                if has_term and has_oh:
                    count += 1
    return count

def identifyEpoxide(mol, element='O'):
    """Counts oxygens in a 3-membered ring with 2 carbons."""
    count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == element:
            c_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'C']
            if len(c_neighbors) == 2 and len(atom.GetNeighbors()) == 2:
                if mol.GetBondBetweenAtoms(c_neighbors[0].GetIdx(), c_neighbors[1].GetIdx()):
                    count += 1
    return count

def identifyDisulfide(mol):
    """Counts bonds between 2 sulfur atoms."""
    count = 0
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetSymbol() == 'S' and bond.GetEndAtom().GetSymbol() == 'S':
            count += 1
    return count
