"""Runs the whole .mol -> log pipeline on small molecules built with RDKit, with thermal noise and
scripted events, so that many kinds of atoms and functional groups are covered."""
import os
import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from breakdown.auto_detect_isomerization import AutoDetectIsomerization
from breakdown.isomer_analyser import analyseIsomer
from conftest import writeMol


def embed(smiles):
    """Get atom symbols and 3D coordinates of a molecule with explicit hydrogens."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol, [a.GetSymbol() for a in mol.GetAtoms()], mol.GetConformer().GetPositions()


def graphKey(mol):
    """Canonical SMILES of only the connectivity, hydrogens explicit."""
    mol = Chem.RWMol(Chem.MolFromSmiles(Chem.MolToSmiles(mol), sanitize=False))
    for bond in mol.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
    for atom in mol.GetAtoms():
        atom.SetNoImplicit(True)
        atom.SetNumRadicalElectrons(0)
        atom.SetNumExplicitHs(0)
        atom.SetIsAromatic(False)
    return Chem.MolToSmiles(mol)


def trajectory(positions, nSteps=30, noise=0.03, seed=0):
    rng = np.random.default_rng(seed)
    return positions + rng.normal(0, noise, (nSteps,) + positions.shape)


def analyse(tmp_path, symbols, coords, name='m.mol', **kwargs):
    writeMol(os.path.join(tmp_path, name), symbols, coords)
    analyser = AutoDetectIsomerization(str(tmp_path), name, **kwargs)
    return open(analyser.autoAnalyseFile()).read().splitlines(), analyser


def firstSmiles(lines):
    return next(line.split()[-1] for line in lines if line.startswith('At step 0'))


MOLECULES = {
    'ethanol': ('CCO', 'Hydroxyl group count: 1', 'C2H6O'),
    'acetone': ('CC(=O)C', 'Ketone count: 1', 'C3H6O'),
    'acetaldehyde': ('CC=O', 'Aldehyde count: 1', 'C2H4O'),
    'acetic acid': ('CC(=O)O', 'Carboxyl count: 1', 'C2H4O2'),
    'diethyl ether': ('CCOCC', 'Ether count: 1', 'C4H10O'),
    'ethylene oxide': ('C1CO1', 'Epoxide count: 1', 'C2H4O'),
    'ethanethiol': ('CCS', 'Thiol group count: 1', 'C2H6S'),
    'thioacetone': ('CC(=S)C', 'Thioketone count: 1', 'C3H6S'),
    'thioacetaldehyde': ('CC=S', 'Thioaldehyde count: 1', 'C2H4S'),
    'dimethyl sulfide': ('CSC', 'Thioether count: 1', 'C2H6S'),
    'thiirane': ('C1CS1', 'Thiirane count: 1', 'C2H4S'),
    'pyridine': ('c1ccncc1', 'Number of 6-carbon rings: 1', 'C5H5N'),
    'benzene': ('c1ccccc1', 'Number of 6-carbon rings: 1', 'C6H6'),
    'methylamine': ('CN', None, 'CH5N'),
    'fluoroethane': ('CCF', None, 'C2H5F'),
    'trifluoroethane': ('FC(F)(F)C', None, 'C2H3F3'),
    'chlorobenzene': ('Clc1ccccc1', 'Number of 6-carbon rings: 1', 'C6H5Cl'),
    'dichloroethane': ('ClCCCl', None, 'C2H4Cl2'),
    'bromoethane': ('CCBr', None, 'C2H5Br'),
    'iodomethane': ('CI', None, 'CH3I'),
    'acetyl chloride': ('CC(Cl)=O', None, 'C2H3ClO'),
    'dimethyl disulfide': ('CSSC', 'Disulfide count: 1', 'C2H6S2'),
    'chloro thioether': ('ClCSC', 'Thioether count: 1', 'C2H5ClS'),
}


@pytest.mark.parametrize('name', MOLECULES)
def test_static_molecule(tmp_path, name):
    smiles, group, formula = MOLECULES[name]
    mol, symbols, positions = embed(smiles)
    lines, analyser = analyse(tmp_path, symbols, trajectory(positions))
    assert not any('Isomerization' in l or 'Fragmentation' in l for l in lines[:-1]), 'noise made an event'
    assert lines[-1] == f'Final molecules present: {formula}'
    assert graphKey(Chem.MolFromSmiles(firstSmiles(lines), sanitize=False)) == graphKey(mol)
    if group:
        # The analysis is only written after events, so call it on the step 0 string
        assert group in analyseIsomer(firstSmiles(lines))
    assert analyser.moleculeElements == sorted(set(symbols), key=lambda e: ({'C': 0, 'H': 1}.get(e, 2), e))


def test_hydrogen_shift_keto_enol(tmp_path):
    """Acetaldehyde to vinyl alcohol: an H moves from carbon to the oxygen."""
    mol, symbols, positions = embed('CC=O')
    carbon = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'C' and a.GetDegree() == 3 and
              any(n.GetSymbol() == 'O' for n in a.GetNeighbors())][0]
    oxygen = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'O'][0]
    hydrogen = [n.GetIdx() for n in mol.GetAtomWithIdx(carbon).GetNeighbors() if n.GetSymbol() == 'H'][0]
    coords = trajectory(positions, 40)
    outward = positions[oxygen] - positions[carbon]
    coords[20:, hydrogen] = positions[oxygen] + 0.96 * outward / np.linalg.norm(outward)
    lines, _ = analyse(tmp_path, symbols, coords)
    events = [l for l in lines if l.startswith('At step') and 'Isomerization' in l]
    assert len(events) == 1 and events[0].startswith('At step 20:')
    assert 'Hydroxyl group count: 1' in lines[lines.index(events[0]) + 1]
    assert lines[-1] == 'Final molecules present: C2H4O'


def test_fragmentation_and_final_formulas(tmp_path):
    """Acetic acid loses its OH group, then the pieces are named by formula."""
    mol, symbols, positions = embed('CC(=O)O')
    oh = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'O' and a.GetTotalNumHs(includeNeighbors=True) == 1]
    hydroxylH = [n.GetIdx() for n in mol.GetAtomWithIdx(oh[0]).GetNeighbors() if n.GetSymbol() == 'H']
    moved = oh + hydroxylH
    coords = trajectory(positions, 40)
    coords[15:, moved] += np.array([12.0, 0, 0])
    lines, _ = analyse(tmp_path, symbols, coords)
    assert any(l.startswith('At step 15: Fragmentation') for l in lines)
    assert lines[-1] == 'Final molecules present: C2H3O, HO'


def test_nitrogen_and_atom_order(tmp_path):
    """Hydrogens first, heteroatoms scattered through the file, gives the same graph."""
    mol, symbols, positions = embed('CC(=O)N')
    order = list(np.random.default_rng(3).permutation(len(symbols)))
    shuffledSymbols = [symbols[i] for i in order]
    lines, analyser = analyse(tmp_path, shuffledSymbols, trajectory(positions[order]))
    assert graphKey(Chem.MolFromSmiles(firstSmiles(lines), sanitize=False)) == graphKey(mol)
    assert lines[-1] == 'Final molecules present: C2H5NO'


@pytest.mark.parametrize('smiles, halogen, formula', [('CCCl', 'Cl', 'C2H5, Cl'), ('CCBr', 'Br', 'C2H5, Br'),
                                                       ('CCF', 'F', 'C2H5, F')])
def test_halogen_leaves(tmp_path, smiles, halogen, formula):
    mol, symbols, positions = embed(smiles)
    index = symbols.index(halogen)
    coords = trajectory(positions, 30)
    coords[12:, index] += np.array([15.0, 0, 0])
    lines, _ = analyse(tmp_path, symbols, coords)
    assert any(l.startswith('At step 12: Fragmentation') for l in lines)
    assert lines[-1] == f'Final molecules present: {formula}'
