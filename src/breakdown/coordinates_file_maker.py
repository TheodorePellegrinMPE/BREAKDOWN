from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
from breakdown.auto_detect_utils import ensureFolderExists
from breakdown.molecule_utils import moleculeFromSmiles

class CoordinatesFileMaker:
    """
    A class for generating 3D coordinates for molecules with explicit hydrogens.
    """
    def __init__(self, folder):
        """
        Initialize the CoordinatesFileMaker with a folder.

        Args:
            folder (str): Folder to save files to.
        """
        self.mol = None
        self.forceFieldName = ''
        self.folder = folder
        ensureFolderExists(self.folder)
        RDLogger.DisableLog('rdApp.warning')

    def makeXYZForSmiles(self, smiles, isomerName):
        """Make an xyz file for a molecule, if it fails a message is printed.

        Args:
            smiles (str): SMILES string representing the molecule.
            isomerName (str): Name of molecule, for saving.

        """
        try:
            self.mol = moleculeFromSmiles(smiles)
            energy = self.generate3DCoordinates()
            self.saveXYZFile(isomerName, energy)
        except (RuntimeError, ValueError) as error:
            print(f'Failed to generate 3D coordinates for {isomerName}: {error}')

    def generate3DCoordinates(self):
        """
        Generate 3D coordinates for the molecule and optimise them with the MMFF
        force field, or UFF for molecules MMFF has no parameters for, such as radicals.

        Returns:
            float: Energy of the molecule in kcal/mol in the force field that was used.

        Raises:
            ValueError: Raised if no 3D coordinates could be embedded.
        """
        if AllChem.EmbedMolecule(self.mol, randomSeed=42) != 0:
            if AllChem.EmbedMolecule(self.mol, randomSeed=42, useRandomCoords=True) != 0:
                raise ValueError('3D embedding failed')
        if AllChem.MMFFHasAllMoleculeParams(self.mol):
            self.forceFieldName = 'MMFF94'
            forceField = AllChem.MMFFGetMoleculeForceField(
                self.mol, AllChem.MMFFGetMoleculeProperties(self.mol))
        else:
            self.forceFieldName = 'UFF'
            forceField = AllChem.UFFGetMoleculeForceField(self.mol)
        forceField.Initialize()
        forceField.Minimize(maxIts=2000)
        return forceField.CalcEnergy()

    def saveXYZFile(self, isomerName, energy):
        """
        Save atomic coordinates to an XYZ file.

        Args:
            isomerName (str): Name of molecule, for saving.
            energy (float): Force field energy, written in the comment line.
        """
        fileName = f'{self.folder}/{isomerName}.xyz'
        atomicInfo = self.getAtomicInformation()
        with open(fileName, 'w') as file:
            file.write(f"{len(atomicInfo)}\n")
            file.write(f"\tEnergy ({self.forceFieldName}, kcal/mol): {energy}\n")
            for atom in atomicInfo:
                file.write(f"{atom['element']:2} {atom['x']: .6f} {atom['y']: .6f} {atom['z']: .6f}\n")

    def getAtomicInformation(self):
        """
        Get atomic information (element, x, y, z) for the molecule, with hydrogens last.

        Returns:
            list of dicts: Atomic information for each atom.
        """
        positions = self.mol.GetConformer().GetPositions()
        symbols = [atom.GetSymbol() for atom in self.mol.GetAtoms()]

        atomic_info = []
        for atom_idx, (pos, symbol) in enumerate(zip(positions, symbols)):
            x, y, z = pos
            atomic_info.append({
                'atomIndex': atom_idx + 1,
                'element': symbol,
                'x': x,
                'y': y,
                'z': z
            })
        atomic_info.sort(key=lambda x: (x['element'] == 'H', x['element']))

        return atomic_info
