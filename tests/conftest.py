import numpy as np
import pytest

# A 6 atom chain-like toy molecule: C-C bonded pair, an H on each C, and an S
# bonded to the first carbon.
ATOMS = ['C', 'C', 'S', 'H', 'H']
POSITIONS = np.array([[0.0, 0.0, 0.0],
                      [1.5, 0.0, 0.0],
                      [-1.8, 0.0, 0.0],
                      [0.0, 1.1, 0.0],
                      [1.5, -1.1, 0.0]])


@pytest.fixture
def atomTypes():
    return np.array(ATOMS)


@pytest.fixture
def frames():
    """Function that repeats the toy molecule for a number of steps."""
    def make(nSteps):
        return np.tile(POSITIONS, (nSteps, 1, 1))
    return make


def writeMol(path, atoms, coords):
    """Write coordinates of shape (steps, atoms, 3) in the deMon .mol layout."""
    with open(path, 'w') as file:
        for step in coords:
            file.write(f'{len(atoms)}\n BOMD INFORMATION     -14.658  0.18\n')
            for atom, (x, y, z) in zip(atoms, step):
                file.write(f'  {atom:2}{x:14.6f}{y:12.6f}{z:12.6f}   0.0  0.0  0.0  0.0\n')
