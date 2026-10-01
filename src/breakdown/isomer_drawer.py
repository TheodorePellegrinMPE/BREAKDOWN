from rdkit import Chem
from rdkit.Chem import Draw
from breakdown.auto_detect_utils import ensureFolderExists
from breakdown.molecule_utils import moleculeFromSmiles

class IsomerDrawer:
    def __init__(self, saveFolder):
        """Initialize IsomerDrawer.

        Args:
            saveFolder (str): Path to folder where images will be saved.
        """
        self.saveFolder = saveFolder
        ensureFolderExists(self.saveFolder)

    def drawAndSaveSmiles(self, smiles, isomerName):
        """Draw a molecule from a SMILES string and save it as an image.

        Args:
            smiles (str): SMILES string representing the molecule.
            isomerName (str): Name for the isomer.

        Raises:
            ValueError: Raised if the SMILES string is invalid.
        """
        mol = moleculeFromSmiles(smiles)
        img = Draw.MolToImage(mol, kekulize=True, wedgeBonds=True)
        img.save(f'{self.saveFolder}/{isomerName}.png')
