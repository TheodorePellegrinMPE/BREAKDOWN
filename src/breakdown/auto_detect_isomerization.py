import re
import os
import warnings
import numpy as np
from breakdown.bond_adder import BondAdder
from breakdown.create_smiles_from_dataframe import analyseTrajectory
from breakdown.auto_detect_utils import countLetters, processSmileToFragment

NUMBER_OR_WORD_REGEX = re.compile(r'-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?|[A-Za-z]+')


class AutoDetectIsomerization:
    def __init__(self, folder, fileName, moleculeSize=None, moleculeElements=None,
                 formationDistance=1.63, breakDistance=2.5, minLifetime=1):
        """Initialize AutoDetectIsomerization object.

        Args:
            folder (str): Path to folder with .mol file.
            fileName (str): Name of .mol file.
            moleculeSize (int or None): Number of atoms in starting molecule, if None
                it is read from the file.
            moleculeElements (list or None): Elements present in the starting molecule,
                if None they are read from the file.
            formationDistance (float): Distance in angstroms where a new bond will form.
            breakDistance (float): Distance in angstroms where a bond will dissociate.
            minLifetime (int): Number of steps a new molecule state has to last to be reported.
        """
        self.moleculeSize = moleculeSize
        self.moleculeElements = list(moleculeElements) if moleculeElements else None
        self.folder = folder
        self.fileName = fileName
        self.formationDistance = formationDistance
        self.breakDistance = breakDistance
        self.minLifetime = minLifetime

    def autoAnalyseFile(self):
        """Analyse data in a file, make a log file of SMILES strings at each step.
        Returns nothing if it is not a valid mol file.

        Returns:
            str or None: Path of the log file that was made, None if the file was skipped.
        """
        file = os.path.join(self.folder, self.fileName)
        if (not os.path.isfile(file) or not self.fileName.endswith('.mol') or
                'keep' in self.fileName):
            print(f'invalid file: {file} (expected an existing .mol file)')
            return None
        saveName = f'{os.path.splitext(self.fileName)[0]}_log.txt'
        saveFile = os.path.join(self.folder, saveName)
        fileString = self.processFile(file)
        fragmentsString = self.getFragments(fileString)
        stringToSave = f'In file: {self.fileName} \n' + fileString + fragmentsString
        self.saveStringAsFile(stringToSave, saveFile)
        return saveFile

    def processFile(self, file):
        """Analyse data in a file by finding the bonds between atoms at each time
        step, create a molecule from these and analyse how this molecule changes
        over the steps.

        Args:
            file (str): Full path URL of file to be read.

        Returns:
            str: String output of analyseTrajectory, contains SMILES
            strings for the molecule at each step and info about isomerization
            and/or fragmentation.
        """
        atomTypes, coords = self.readTrajectory(file)
        bondAdder = BondAdder(len(atomTypes), self.formationDistance, self.breakDistance)
        firstAtoms, secondAtoms, bonded = bondAdder.findBonds(atomTypes, coords)
        return analyseTrajectory(atomTypes, firstAtoms, secondAtoms, bonded, self.minLifetime)

    def readTrajectory(self, file):
        """Read a file into atom types and coordinates, and fill in the molecule
        size and elements if they were not given.

        Args:
            file (str): Full path URL of file to be read.

        Returns:
            Tuple[np.ndarray, np.ndarray]: Tuple containing:
            - atomTypes (np.ndarray): Element symbol of each atom in the molecule.
            - coords (np.ndarray): Coordinates with shape (steps, atoms, 3).

        Raises:
            ValueError: Raised if the file has no atoms or if the atoms do not
                match the molecule size.
        """
        headerSize, symbols, xyz = self.readFile(file)
        if self.moleculeSize is None:
            self.moleculeSize = headerSize
        elif headerSize is not None and headerSize != self.moleculeSize:
            raise ValueError(f'{file} has {headerSize} atoms per step, but the molecule size is '
                             f'{self.moleculeSize}.')
        if not self.moleculeSize or len(symbols) < self.moleculeSize:
            raise ValueError(f'No complete step of {self.moleculeSize} atoms found in {file}.')
        nSteps, extra = divmod(len(symbols), self.moleculeSize)
        if extra:
            warnings.warn(f'{file} ends with an incomplete step of {extra} atoms, it is ignored.')
        symbols = symbols[:nSteps * self.moleculeSize].reshape(nSteps, self.moleculeSize)
        coords = xyz[:nSteps * self.moleculeSize].reshape(nSteps, self.moleculeSize, 3)
        if (symbols != symbols[0]).any():
            raise ValueError(f'The order of the atoms changes between steps in {file}, '
                             'check that the molecule size is correct.')
        atomTypes = symbols[0]
        self.fillMoleculeElements(atomTypes)
        return atomTypes, coords

    def fillMoleculeElements(self, atomTypes):
        """Set the elements in the molecule from the atoms if they were not given,
        C and H first followed by the rest alphabetically. If they were given, warn if
        elements that are present are missing from them.

        Args:
            atomTypes (np.ndarray): Element symbol of each atom in the molecule.
        """
        present = sorted(set(atomTypes), key=lambda e: ({'C': 0, 'H': 1}.get(e, 2), e))
        if self.moleculeElements is None:
            self.moleculeElements = present
        else:
            missing = [e for e in present if e not in self.moleculeElements]
            if missing:
                warnings.warn(f'Elements {missing} are in {self.fileName} but not in the '
                              'list of elements, they will be missing from formulas.')

    def readFile(self, file):
        """Read a file, filter out BOMD INFORMATION and molecule size lines.
        The atom lines are then split into element and coordinates.

        Args:
            file (str): Full path URL of file to be read.

        Returns:
            Tuple[int or None, np.ndarray, np.ndarray]: Tuple containing the number of
            atoms given by the first line of the file, an array of the element of every
            atom line, and an array of their coordinates with shape (atom lines, 3).
        """
        symbols = []
        xyz = []
        headerSize = None
        with open(file) as f:
            for lineNumber, line in enumerate(f):
                if lineNumber == 0 and line.strip().isdigit():
                    headerSize = int(line)
                elif line.startswith('  ') and line[2:3].isupper():
                    symbol, coordinates = self.processLine(line)
                    symbols.append(symbol)
                    xyz.append(coordinates)
        return headerSize, np.array(symbols, dtype=str), np.array(xyz, dtype=float).reshape(-1, 3)

    def processLine(self, line):
        """Given a string line of an atom, get the element and the coordinates.
        Lines are split on spaces, or if that fails on dashes as well, for when negative
        4 digit coordinate numbers cause the columns to not be separated by a space.

        Args:
            line (str): String of a line from .mol file.

        Returns:
            Tuple[str, list]: Element symbol and the x, y and z coordinates.
        """
        parts = line.split()
        try:
            return parts[0], [float(parts[1]), float(parts[2]), float(parts[3])]
        except (ValueError, IndexError):
            parts = NUMBER_OR_WORD_REGEX.findall(line.strip())
            return parts[0], [float(parts[1]), float(parts[2]), float(parts[3])]

    def saveStringAsFile(self, stringToSave, saveFile):
        """Save the string stringToSave at the URL saveFile.

        Args:
            stringToSave (str): String to be saved in a file.
            saveFile (str): Full path URL for the file to be saved at.
        """
        with open(saveFile, 'w') as file:
            file.write(stringToSave)

    def getFragments(self, fileString):
        """Make a string listing all fragments present at the last step of the
        run in chemical formula format.

        Args:
            fileString (str): String of SMILES and isomerization details for a run.

        Returns:
            str: String readably listing the fragments for the save file.
        """
        lastStepLine = self.getLastStepLine(fileString)
        lastSmiles = lastStepLine.split()[-1]
        fragmentSmiles = lastSmiles.split('.')
        fragmentSmiles = sorted(fragmentSmiles,
                                key=countLetters, reverse=True)
        fragments = []
        for smile in fragmentSmiles:
            fragment = processSmileToFragment(smile, self.moleculeElements)
            fragments.append(fragment)
        fragmentsString = 'Final molecules present: ' + ', '.join(fragments)
        return fragmentsString

    def getLastStepLine(self, fileString):
        """Gets the last line of the string to start with At step.

        This is done so that getFragments works if the last line is details about
        isomerization.

        Args:
            fileString (str): String of SMILES and isomerization details for a run.

        Returns:
            str: Last line of fileString to start with 'At step'.

        Raises:
            ValueError: Raised if there is no such line.
        """
        lines = fileString.split('\n')
        for line in reversed(lines):
            if line.startswith('At step'):
                return line
        raise ValueError(f'No steps were analysed in {self.fileName}.')
