import numpy as np
from rdkit import Chem

# Pairs involving a C-S bond form at a longer distance than the default.
CARBON_SULFUR_FORMATION_DISTANCE = 1.9
# Pairs of elements that are not both in this set form bonds up to this factor times
# the sum of their covalent radii, when that is longer than the formation distance.
# Bonds like C-Cl (1.8 A), C-Br (1.9 A) or S-S (2.05 A) are longer than the default formation distance.
STANDARD_ELEMENTS = ('H', 'C', 'N', 'O')
COVALENT_RADIUS_FACTOR = 1.1
# Number of time steps whose distances are held in memory at once.
CHUNK_SIZE = 4096


class BondAdder:
    def __init__(self, moleculeSize, formationDistance, breakDistance):
        """Initialize BondAdder object.

        Args:
            moleculeSize (int): Number of atoms in starting molecule.
            formationDistance (float): Distance in angstroms where a new bond will form.
            breakDistance (float): Distance in angstroms where a bond will dissociate.
        """
        self.moleculeSize = moleculeSize
        self.formationDistance = formationDistance
        self.breakDistance = breakDistance

    def getPairThresholds(self, atomTypes):
        """Get the formation and break distances of every pair of atoms.

        Args:
            atomTypes (np.ndarray): Element symbol of each atom in the molecule.

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]: Tuple containing:
            - firstAtoms (np.ndarray): Index of the first atom of each pair.
            - secondAtoms (np.ndarray): Index of the second atom of each pair, always larger.
            - formation (np.ndarray): Distance under which each pair forms a bond.
            - breaking (np.ndarray): Distance over which each bonded pair breaks.
        """
        firstAtoms, secondAtoms = np.triu_indices(len(atomTypes), k=1)
        typesA = atomTypes[firstAtoms]
        typesB = atomTypes[secondAtoms]
        formation = np.full(len(firstAtoms), self.formationDistance, dtype=float)
        breaking = np.full(len(firstAtoms), self.breakDistance, dtype=float)

        isCSBond = ((typesA == 'C') & (typesB == 'S')) | ((typesA == 'S') & (typesB == 'C'))
        periodicTable = Chem.GetPeriodicTable()
        for elementA in set(atomTypes):
            for elementB in set(atomTypes):
                if elementA in STANDARD_ELEMENTS and elementB in STANDARD_ELEMENTS:
                    continue
                radiiSum = periodicTable.GetRcovalent(str(elementA)) + periodicTable.GetRcovalent(str(elementB))
                pairs = (typesA == elementA) & (typesB == elementB)
                formation[pairs] = np.maximum(formation[pairs], COVALENT_RADIUS_FACTOR * radiiSum)
        formation[isCSBond] = CARBON_SULFUR_FORMATION_DISTANCE
        breaking = np.maximum(breaking, formation)
        return firstAtoms, secondAtoms, formation, breaking

    def findBonds(self, atomTypes, coords):
        """Find which pairs of atoms are bonded at every time step.

        Two atoms become bonded when closer than the formation distance and stay
        bonded until they are further apart than the break distance. Hydrogen can
        only be bonded to one atom, the closest currently bonded non-hydrogen
        atom, or if there is none, the closest bonded hydrogen of lower index
        (so an H2 pair is bonded once). This is done for all time steps at once,
        only the bond persistence needs a loop over time steps.

        Args:
            atomTypes (np.ndarray): Element symbol of each atom in the molecule.
            coords (np.ndarray): Coordinates with shape (steps, atoms, 3).

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray]: Tuple containing:
            - firstAtoms (np.ndarray): Index of the first atom of each pair.
            - secondAtoms (np.ndarray): Index of the second atom of each pair, always larger.
            - bonded (np.ndarray): Boolean array with shape (steps, pairs), True if
              the pair is bonded at that step.
        """
        nSteps, nAtoms, _ = coords.shape
        firstAtoms, secondAtoms, formation, breaking = self.getPairThresholds(atomTypes)
        bonded = np.empty((nSteps, len(firstAtoms)), dtype=bool)
        pairIndex = np.full((nAtoms, nAtoms), -1, dtype=int)
        pairIndex[firstAtoms, secondAtoms] = np.arange(len(firstAtoms))
        pairIndex[secondAtoms, firstAtoms] = pairIndex[firstAtoms, secondAtoms]
        previous = np.zeros(len(firstAtoms), dtype=bool)

        for start in range(0, nSteps, CHUNK_SIZE):
            chunk = coords[start:start + CHUNK_SIZE]
            diff = chunk[:, firstAtoms, :] - chunk[:, secondAtoms, :]
            distances = np.sqrt((diff ** 2).sum(axis=2))
            positive = distances > 0
            canForm = (distances <= formation) & positive
            canStay = (distances <= breaking) & positive
            bound = np.empty_like(canForm)
            for t in range(len(chunk)):
                previous = canForm[t] | (previous & canStay[t])
                bound[t] = previous
            bonded[start:start + len(chunk)] = self.resolveHydrogenBonds(
                atomTypes, bound, distances, pairIndex)
        return firstAtoms, secondAtoms, bonded

    def resolveHydrogenBonds(self, atomTypes, bound, distances, pairIndex):
        """Limit every hydrogen atom to at most one bond.

        Args:
            atomTypes (np.ndarray): Element symbol of each atom in the molecule.
            bound (np.ndarray): Boolean array (steps, pairs) of pairs close enough to be bonded.
            distances (np.ndarray): Distance of each pair at each step, shape (steps, pairs).
            pairIndex (np.ndarray): Matrix mapping two atom indices to the index of their pair.

        Returns:
            np.ndarray: Copy of bound where each hydrogen only has its closest bond.
        """
        isH = atomTypes == 'H'
        hydrogens = np.flatnonzero(isH)
        heavyAtoms = np.flatnonzero(~isH)
        result = bound.copy()
        for h in hydrogens:
            heavyPairs = pairIndex[h, heavyAtoms]
            lowerHydrogenPairs = pairIndex[h, hydrogens[hydrogens < h]]
            owned = np.concatenate([heavyPairs, lowerHydrogenPairs])
            result[:, owned] = False
            chosenHeavy, hasHeavy = self.closestBound(bound, distances, heavyPairs)
            # A lower hydrogen that is already bonded to a heavy atom does not take a second bond
            lower = hydrogens[hydrogens < h]
            allowedLower = bound[:, lowerHydrogenPairs].copy()
            for column, otherHydrogen in enumerate(lower):
                allowedLower[:, column] &= ~result[:, pairIndex[otherHydrogen, heavyAtoms]].any(axis=1)
            chosenLower, hasLower = self.closestBound(bound, distances, lowerHydrogenPairs, allowedLower)
            rows = np.arange(len(bound))
            useHeavy = hasHeavy
            useLower = ~hasHeavy & hasLower
            result[rows[useHeavy], chosenHeavy[useHeavy]] = True
            result[rows[useLower], chosenLower[useLower]] = True
        return result

    def closestBound(self, bound, distances, pairs, allowed=None):
        """Find the closest pair with a bond at every step, among some pairs.

        Args:
            bound (np.ndarray): Boolean array (steps, pairs) of pairs close enough to be bonded.
            distances (np.ndarray): Distance of each pair at each step, shape (steps, pairs).
            pairs (np.ndarray): Indices of the pairs to consider.
            allowed (np.ndarray or None): Boolean array (steps, len(pairs)) of the pairs that may be chosen at
                each step, the bonded ones if None.

        Returns:
            Tuple[np.ndarray, np.ndarray]: The index of the closest bonded pair
            at each step, and whether there was any bonded pair at that step.
        """
        if len(pairs) == 0:
            empty = np.zeros(len(bound), dtype=int)
            return empty, np.zeros(len(bound), dtype=bool)
        masked = np.where(bound[:, pairs] if allowed is None else allowed, distances[:, pairs], np.inf)
        closest = masked.argmin(axis=1)
        hasBond = np.isfinite(masked[np.arange(len(bound)), closest])
        return pairs[closest], hasBond
