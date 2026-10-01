import numpy as np
from rdkit import Chem
from breakdown.bond_adder import BondAdder
from breakdown.create_smiles_from_dataframe import (analyseTrajectory, getSmilesPerStep,
                                                        checkIsomerization)
from breakdown.molecule_utils import moleculeFromSmiles


def bondPairs(atomTypes, coords, **kwargs):
    adder = BondAdder(len(atomTypes), 1.63, 2.5, **kwargs)
    first, second, bonded = adder.findBonds(atomTypes, coords)
    return first, second, bonded


def bondedSet(first, second, bondedRow):
    return {(int(a), int(b)) for a, b in zip(first[bondedRow], second[bondedRow])}


def test_bonds_of_toy_molecule(atomTypes, frames):
    first, second, bonded = bondPairs(atomTypes, frames(3))
    assert bondedSet(first, second, bonded[0]) == {(0, 1), (0, 2), (0, 3), (1, 4)}


def test_bond_persists_until_break_distance(frames):
    atomTypes = np.array(['C', 'C'])
    distances = [1.5, 2.0, 2.4, 2.6, 2.0, 1.5]
    coords = np.array([[[0, 0, 0], [d, 0, 0]] for d in distances], dtype=float)
    _, _, bonded = bondPairs(atomTypes, coords)
    # Stays bonded while stretched to 2.4, breaks at 2.6, does not reform until 1.63
    assert bonded[:, 0].tolist() == [True, True, True, False, False, True]


def test_carbon_sulfur_forms_at_longer_distance():
    atomTypes = np.array(['C', 'S', 'C', 'C'])
    coords = np.array([[[0, 0, 0], [1.8, 0, 0], [0, 5, 0], [1.8, 5, 0]]], dtype=float)
    _, _, bonded = bondPairs(atomTypes, coords)
    assert bonded[0].tolist().count(True) == 1  # C-S at 1.8 bonds, C-C at 1.8 does not


def test_hydrogen_bonds_only_closest_atom(atomTypes, frames):
    coords = frames(1).copy()
    coords[0, 3] = [0.7, 0.5, 0.0]  # H between both carbons, closer to C0
    first, second, bonded = bondPairs(atomTypes, coords)
    hydrogenBonds = [pair for pair in bondedSet(first, second, bonded[0]) if 3 in pair]
    assert len(hydrogenBonds) == 1


def test_atom_order_does_not_change_graph(atomTypes, frames):
    order = [3, 0, 4, 1, 2]  # hydrogens before their heavy atoms
    coords = frames(2)
    smiles = getSmilesPerStep(atomTypes, *bondPairs(atomTypes, coords))
    shuffledTypes = atomTypes[order]
    shuffled = getSmilesPerStep(shuffledTypes, *bondPairs(shuffledTypes, coords[:, order]))
    assert smiles == shuffled


def test_h2_fragment_and_recombination(atomTypes, frames):
    coords = frames(4).copy()
    coords[1:3, 3] = [10.0, 0.0, 0.0]  # H flies away at steps 1 and 2
    coords[1:3, 4] = [10.5, 0.0, 0.0]  # and forms H2 with the other H
    first, second, bonded = bondPairs(atomTypes, coords)
    log = analyseTrajectory(atomTypes, first, second, bonded)
    events = [line.split()[3] for line in log.splitlines()
              if line.startswith('At step') and len(line.split()) == 5]
    assert events == ['Fragmentation', 'Recombination']


def test_min_lifetime_removes_flicker(atomTypes, frames):
    coords = frames(10).copy()
    coords[4, 3] = [10.0, 0, 0]  # one step blip of an H leaving
    args = bondPairs(atomTypes, coords)
    noisy = analyseTrajectory(atomTypes, *args)
    smoothed = analyseTrajectory(atomTypes, *args, minLifetime=3)
    assert 'Fragmentation' in noisy
    assert 'Fragmentation' not in smoothed


def test_no_implicit_hydrogens_added(atomTypes, frames):
    smiles = getSmilesPerStep(atomTypes, *bondPairs(atomTypes, frames(1)))[0]
    mol = moleculeFromSmiles(smiles)
    assert mol.GetNumAtoms() == len(atomTypes)
    assert sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()) > 0


def test_molecule_from_smiles_gives_double_bonds():
    # Only connectivity is known for these, sp2 carbons should be joined by double bonds, not left as radicals
    smiles = '[H]C1C([H])C([H])C([H])C([H])C1[H]'
    mol = moleculeFromSmiles(smiles)
    assert sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()) == 0
    assert mol.GetNumAtoms() == 12
    # A ring with an odd number of sp2 carbons has to keep one radical, and gets no extra hydrogen
    radical = moleculeFromSmiles('[H]C1C([H])C([H])C([H])C1[H]')
    assert radical.GetNumAtoms() == 10
    assert sum(atom.GetNumRadicalElectrons() for atom in radical.GetAtoms()) == 1
    # The H inside [SH] is an RDKit artifact
    assert moleculeFromSmiles('[H]C1[SH]C1').GetNumAtoms() == 4


def test_check_isomerization_labels():
    assert checkIsomerization('C', '').strip() == ''
    assert checkIsomerization('C.C', 'C') == 'Fragmentation'
    assert checkIsomerization('C', 'C.C') == 'Recombination'
    assert checkIsomerization('CC', 'C=C') == 'Isomerization'


def test_hydrogen_never_has_two_bonds():
    # H1 is bonded to the carbon, H2 sits close to H1 but further from the carbon: it must not bond to H1
    atomTypes = np.array(['C', 'H', 'H'])
    coords = np.array([[[0, 0, 0], [1.1, 0, 0], [1.5, 0.9, 0]]], dtype=float)
    first, second, bonded = bondPairs(atomTypes, coords)
    pairs = bondedSet(first, second, bonded[0])
    degree = {i: sum(i in pair for pair in pairs) for i in range(3)}
    assert degree[1] <= 1 and degree[2] <= 1
    assert (0, 1) in pairs
