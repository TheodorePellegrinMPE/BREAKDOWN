import numpy as np
from rdkit import Chem
from breakdown.isomer_analyser import analyseIsomer

# Atoms that are drawn double bonded to a carbon when they have no other bond.
TERMINAL_DOUBLE_BOND_ATOMS = ('O', 'S')
BLANK_EVENT = ' ' * 13
MANY_FRAGMENTS_LIMIT = 20


def analyseTrajectory(atomTypes, firstAtoms, secondAtoms, bonded, minLifetime=1,
                      manyFragmentsLimit=MANY_FRAGMENTS_LIMIT):
    """Analyse a molecule at every step to construct a string with SMILES
    representations and information about all isomerizations and fragmentations.

    RDKit is only asked for a SMILES string when the bonds differ from a
    bond pattern that has already been seen, most steps do not change the bonds.

    Args:
        atomTypes (np.ndarray): Element symbol of each atom in the molecule.
        firstAtoms (np.ndarray): Index of the first atom of each pair of atoms.
        secondAtoms (np.ndarray): Index of the second atom of each pair of atoms.
        bonded (np.ndarray): Boolean array (steps, pairs), True where a pair is bonded.
        minLifetime (int): Number of steps a new molecule state has to last for
            it to be reported, shorter lived states are treated as noise and
            replaced by the state before them.
        manyFragmentsLimit (int): Analysis stops after this many steps
            with many carbon containing fragments.

    Returns:
        str: Result of analyzing every step into one readable string which will be
        sent back to and saved by auto_detect_isomerizations.py.
    """
    smilesPerStep = getSmilesPerStep(atomTypes, firstAtoms, secondAtoms, bonded, minLifetime)
    analysisCache = {}
    reportLines = []
    prevSmiles = ''
    endCount = 0
    for step, smiles in enumerate(smilesPerStep):
        if checkManyFragments(smiles):
            endCount += 1
        eventString = checkIsomerization(smiles, prevSmiles)
        reportString = reportSmiles(step, eventString, smiles)
        if eventString.strip():
            if smiles not in analysisCache:
                analysisCache[smiles] = analyseIsomer(smiles)
            reportString += '\n' + analysisCache[smiles]
        reportLines.append(reportString)
        prevSmiles = smiles
        if endCount > manyFragmentsLimit:
            break
    return '\n'.join(reportLines) + '\n'


def getSmilesPerStep(atomTypes, firstAtoms, secondAtoms, bonded, minLifetime=1):
    """Make a SMILES string for every step of a trajectory.

    Args:
        atomTypes (np.ndarray): Element symbol of each atom in the molecule.
        firstAtoms (np.ndarray): Index of the first atom of each pair of atoms.
        secondAtoms (np.ndarray): Index of the second atom of each pair of atoms.
        bonded (np.ndarray): Boolean array (steps, pairs), True where a pair is bonded.
        minLifetime (int): Shortest number of steps a state must last to be kept.

    Returns:
        list: SMILES string for every step.
    """
    packed = np.packbits(bonded, axis=1)
    stateIds = np.empty(len(bonded), dtype=int)
    stateBytes = {}
    stateSmiles = []
    for step in range(len(bonded)):
        key = packed[step].tobytes()
        if key not in stateBytes:
            stateBytes[key] = len(stateSmiles)
            pairs = np.flatnonzero(bonded[step])
            stateSmiles.append(makeSmiles(atomTypes, firstAtoms[pairs], secondAtoms[pairs]))
        stateIds[step] = stateBytes[key]
    if minLifetime > 1:
        stateIds = removeShortLivedStates(stateIds, minLifetime)
    return [stateSmiles[stateId] for stateId in stateIds]


def removeShortLivedStates(stateIds, minLifetime):
    """Replace states that last for less than minLifetime steps with the state before.

    Args:
        stateIds (np.ndarray): Id of the molecule state at every step.
        minLifetime (int): Shortest number of steps a state must last to be kept.

    Returns:
        np.ndarray: Ids with short lived states removed.
    """
    cleaned = stateIds.copy()
    changes = np.flatnonzero(cleaned[1:] != cleaned[:-1]) + 1
    starts = np.concatenate([[0], changes])
    ends = np.concatenate([changes, [len(cleaned)]])
    for start, end in zip(starts[1:], ends[1:]):
        # Cleaning earlier runs may have made the run before this one identical to it.
        if end - start < minLifetime and end != len(cleaned):
            cleaned[start:end] = cleaned[start - 1]
    return cleaned


def checkManyFragments(smiles):
    """Check if a molecule has many carbon containing fragments.
    If more than 3 carbon containing fragments are present, the molecule is considered to have many fragments.

    Args:
        smiles (str): SMILES string representing the molecule.

    Returns:
        bool: True if the molecule has many fragments, False otherwise.
    """
    fragments = smiles.split('.')
    carbonFragmentCount = sum('C' in s for s in fragments)
    manyFragments = carbonFragmentCount > 3
    return manyFragments


def makeSmiles(atomTypes, firstAtoms, secondAtoms):
    """Generate a SMILES string representing the connectivity of a molecule.

    The hydrogens are all explicit atoms in the simulation. The bond orders are not
    known, so the string only encodes which atoms are bonded. RDKit writes brackets
    like [SH] around atoms with an unusual number of bonds, the H there is not real.
    See moleculeFromSmiles in molecule_utils.py for turning the string into a
    molecule with sensible valences and no extra hydrogens.

    Args:
        atomTypes (np.ndarray): Element symbol of each atom in the molecule.
        firstAtoms (np.ndarray): Index of the first atom of each bond.
        secondAtoms (np.ndarray): Index of the second atom of each bond.

    Returns:
        str: SMILES string representing the molecule.
    """
    mol = Chem.RWMol()
    for atomType in atomTypes:
        rdAtom = Chem.Atom(str(atomType))
        if atomType == 'H':
            # Otherwise a lone hydrogen atom is written as [HH]
            rdAtom.SetNoImplicit(True)
        mol.AddAtom(rdAtom)
    degrees = np.bincount(np.concatenate([firstAtoms, secondAtoms]), minlength=len(atomTypes))
    for a, b in zip(firstAtoms, secondAtoms):
        bondOrder = getBondOrder(atomTypes[a], atomTypes[b], degrees[a], degrees[b])
        mol.AddBond(int(a), int(b), Chem.BondType.DOUBLE if bondOrder == 2 else Chem.BondType.SINGLE)
    return Chem.MolToSmiles(mol.GetMol())


def getBondOrder(typeA, typeB, degreeA, degreeB):
    """Get the bond order that is used to make the SMILES string. The real order
    can't be known because the structure is not known and the distance is
    unreliable as bonds get stretched, so only bonds that are certainly double are
    drawn as such: a terminal oxygen or sulfur atom bonded to something that is not
    hydrogen. All other bonds are single.

    Args:
        typeA (str): Element symbol of the first atom.
        typeB (str): Element symbol of the second atom.
        degreeA (int): Number of bonds of the first atom.
        degreeB (int): Number of bonds of the second atom.

    Returns:
        int: 2 for a terminal O or S atom attached to a heavy atom, otherwise 1.
    """
    if typeA == 'H' or typeB == 'H':
        return 1
    if typeA in TERMINAL_DOUBLE_BOND_ATOMS and degreeA == 1:
        return 2
    if typeB in TERMINAL_DOUBLE_BOND_ATOMS and degreeB == 1:
        return 2
    return 1


def checkIsomerization(smiles, prevSmiles):
    """Determine if the molecule's SMILES has changed, indicating isomerization,
    fragmentation or recombination of fragments.

    Args:
        smiles (str): SMILES string for the molecule at the current step.
        prevSmiles (str): SMILES string for the molecule at the previous step.

    Returns:
        str: String indicating what happened since the previous step, blank if nothing did.
    """
    if prevSmiles == '' or prevSmiles == smiles:
        return BLANK_EVENT
    numFragments = smiles.count('.')
    prevNumFragments = prevSmiles.count('.')
    if numFragments > prevNumFragments:
        return 'Fragmentation'
    if numFragments < prevNumFragments:
        return 'Recombination'
    return 'Isomerization'


def reportSmiles(step, eventString, smiles):
    """Generate a string containing relevant information about the simulation step,
    the current molecule's SMILES, and if there was isomerization or fragmentation.

    Args:
        step (int): The number of steps since the start of the simulation.
        eventString (str): String, isomerization, fragmentation, recombination or
            blank, specifying if the molecule changed since the last step.
        smiles (str): SMILES string representing the molecule at this step.

    Returns:
        str: Readable string showing information of the molecule at the current step.
    """
    reportString = f'At step {step}: {eventString} {smiles}'
    return reportString
