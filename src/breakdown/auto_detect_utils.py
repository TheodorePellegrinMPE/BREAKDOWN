import numpy as np
import os
import json
from collections import Counter
from functools import lru_cache
from rdkit import Chem

EVENT_LABELS = ('Isomerization', 'Fragmentation', 'Recombination')
NOMINAL_MASSES = {'H': 1, 'C': 12, 'N': 14, 'O': 16, 'S': 32}


def isEventLine(line):
    """Check if a line of a log file reports a change of the molecule.

    Args:
        line (str): Line from a log file.

    Returns:
        bool: True if the line is an At step line with isomerization, fragmentation or recombination.
    """
    return line.startswith('At step') and any(label in line for label in EVENT_LABELS)


@lru_cache(maxsize=None)
def countAtoms(smile):
    """Counts the atoms of each element in a SMILES string. All hydrogens in the
    strings are atoms of their own, an H inside brackets, like [SH], is added by RDKit
    to atoms with an unusual number of bonds and is not counted.

    Args:
        smile (str): SMILES string representation of a molecule.

    Returns:
        dict: Element symbols as keys and atom counts as values.
    """
    mol = Chem.MolFromSmiles(smile, sanitize=False)
    if mol is None:
        raise ValueError(f'Invalid SMILES string: {smile}')
    counts = Counter()
    for atom in mol.GetAtoms():
        counts[atom.GetSymbol()] += 1
    return dict(counts)


def orderElements(elements):
    """Sort element symbols with C first, then H, then the rest alphabetically.

    Args:
        elements (iterable): Element symbols.

    Returns:
        list: Sorted list of the distinct symbols.
    """
    return sorted(set(elements), key=lambda e: ({'C': 0, 'H': 1}.get(e, 2), e))


def processSmileToFragment(smile, moleculeElements=None):
    """
    Converts a SMILES string to a chemical formula string.

    Args:
        smile (str): SMILES string representation of a molecule.
        moleculeElements (list, optional): List of chemical elements, sets which
            elements are in the formula and their order. All elements present
            in the SMILES are used if None.

    Returns:
        str: Chemical formula string for matching smile.
    """
    counts = countAtoms(smile)
    if moleculeElements is None:
        moleculeElements = orderElements(counts)
    fragment = ''
    for element in moleculeElements:
        numElement = counts.get(element, 0)
        if numElement == 1:
            fragment += element
        elif numElement > 1:
            fragment += f'{element}{numElement}'
    return fragment


def calculateNominalMass(smile):
    """Calculate the mass of a SMILES string as a whole number, from the most
    common isotope of each atom.

    Args:
        smile (str): SMILES string representation of a molecule.

    Returns:
        int: Nominal mass.
    """
    table = Chem.GetPeriodicTable()
    mass = 0
    for element, count in countAtoms(smile).items():
        mass += count * NOMINAL_MASSES.get(element, round(table.GetMostCommonIsotopeMass(element)))
    return mass


def countLetters(string, specificLetter=None):
    """
    Counts the total number of English alphabet letters in a string.
    If a specific letter is provided, counts only that letter.

    Args:
        string (str): The string to be counted.
        specific_letter (str, optional): The specific letter to count.
            If None, counts all English alphabet letters.

    Returns:
        int: The total number of letters in the string or the count of the specific letter.
    """
    if specificLetter is not None:
        specific_count = sum(1 for char in string if char.lower() == specificLetter.lower() and char.isalpha())
        return specific_count
    else:
        letter_count = sum(1 for char in string if char.isalpha())
        return letter_count

def formatIntegerList(numbers):
    """
    Formats a list of numbers as a string with complete ranges.

    Args:
        numbers (list): List of integers.

    Returns:
        str: String listing the numbers.
    """
    if not numbers:
        return ''
    numbers = sorted(numbers)
    result = [str(numbers[0])]
    InRange = False
    for i in range(1, len(numbers)):
        if numbers[i] == numbers[i - 1] + 1:
            InRange = True
            continue
        if InRange:
            result[-1] += f'-{numbers[i - 1]}'
            InRange = False
        result.append(str(numbers[i]))
    if InRange:
        result[-1] += f'-{numbers[-1]}'
    integerString = ', '.join(result)
    return integerString

def makeArrayOfEmptyLists(length):
    """
    Builds a numpy array of a given length, where each element is an empty list.

    Args:
        length (int): Integer size of the desired array.

    Returns:
        numpy.ndarray: Numpy array of empty lists.
    """
    array = np.empty(length, dtype=object)
    for i in range(len(array)):
        array[i] = []
    return array

def ensureFolderExists(folderPath):
    """
    Creates a folder at folderPath if it does not exist.

    Args:
        folderPath (str): Full path to the desired folder.
    """
    if not os.path.exists(folderPath):
        os.makedirs(folderPath)

def extractToDict(inputString, startSubstring, endSubstring=None):
    """
    Extracts a substring from the input string that is located between the specified start
    and end substrings (or until the end of the input string if no end substring is provided),
    and converts the extracted substring into a dictionary assuming a key-value pair format.

    Parameters:
    - inputString (str): The input string from which to extract the substring.
    - startSubstring (str): The starting substring to search for.
    - endSubstring (str or None): The ending substring to search for. If None, extraction will
                                   continue until the end of the input string.

    Returns:
    - dict or None: A dictionary containing key-value pairs extracted from the substring
                   between the start and end substrings, or until the end of the input string
                   if no end substring is provided. Returns None if no match is found.
    """
    if startSubstring not in inputString:
        raise ValueError(f'The heading "{startSubstring}" is not in the log, was it made by an older version?')
    extractedString = inputString.split(startSubstring)[1]
    if endSubstring:
        extractedString = extractedString.split(endSubstring)[0]
    extractedString = extractedString.lstrip('\n').rstrip('\n')
    resultDict = json.loads(extractedString)
    return resultDict

def getFragmentNames(smiles, name):
    """Get names for all fragments of a smiles string.

    Args:
        smiles (str): smiles string of a state.
        name (str): name of the molecule state.

    Returns:
        list: list of fragment smiles.
        list: list of fragment names.
    """
    fragmentSmiles = smiles.split('.')
    fragmentSmiles = sorted(fragmentSmiles, key=countLetters, reverse=True)
    multipleLargeFragments = False
    fragments = []
    isomerNames = []
    if len(fragmentSmiles) > 1:
        multipleLargeFragments = countLetters(fragmentSmiles[1], 'C') > 2
    for index, fragment in enumerate(fragmentSmiles):
        isomerName = fragment
        if multipleLargeFragments:
            if countLetters(fragment, 'C') > 2:
                isomerName = f'{name}_{index}'
        else:
            if index == 0:
                isomerName = name
        fragments.append(fragment)
        isomerNames.append(isomerName)
    return fragments, isomerNames


def classifyRun(logString):
    """Say if a run had an isomerization, and if it ended fragmented.

    A run that loses a hydrogen for a moment and takes it back (fragmentation followed by recombination)
    is not fragmented, only a run whose last step is in more than one piece is.

    Args:
        logString (str): Log of one run.

    Returns:
        Tuple[bool, bool, bool]: Whether the run had any isomerization, whether it ever fragmented,
        and whether it is fragmented at the end.
    """
    isomerized = False
    everFragmented = False
    lastSmiles = ''
    for line in logString.split('\n'):
        if line.startswith('At step'):
            parts = line.split()
            lastSmiles = parts[-1]
            if len(parts) == 5:
                isomerized = isomerized or parts[3] == 'Isomerization'
                everFragmented = everFragmented or parts[3] == 'Fragmentation'
    return isomerized, everFragmented, '.' in lastSmiles
