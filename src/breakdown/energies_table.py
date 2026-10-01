import os
import csv

class EnergiesTable:
    def __init__(self, folder):
        """Initialize the EnergiesTable with the folder.
        """
        self.folder = folder
        self.energyData = []

    def processOutFiles(self):
        """Process all .out files in the folder and extract DFTB electronic energy.
        """
        outFiles = [file for file in os.listdir(self.folder) if file.endswith('.out')]
        for outFile in outFiles:
            energyValue = self.extractEnergyFromFile(outFile)
            if energyValue is not None:
                self.energyData.append({'File': outFile, 'Energy [Hartree]': energyValue})

    def extractEnergyFromFile(self, fileName):
        """Extract DFTB electronic energy from a given .out file.
        """
        outFilepath = os.path.join(self.folder, fileName)
        energyKeyword = 'DFTB electronic energy    [Hartree] ='
        with open(outFilepath, 'r') as file:
            for line in file:
                if energyKeyword in line:
                    energyValue = float(line.split('=')[-1].strip())
                    return energyValue
        return None

    def saveToCsv(self, outputCsv='energies_table.csv'):
        """Save the extracted energy data to a CSV file in alphabetical order.
        """
        if not self.energyData:
            print('No "DFTB electronic energy" lines found in the .out files, they are only written by single point and optimization runs, not MD runs.')
            return
        csvPath = os.path.join(self.folder, outputCsv)
        sortedData = sorted(self.energyData, key=lambda x: x['File'].lower())
        with open(csvPath, 'w', newline='') as csvfile:
            fieldnames = ['File', 'Energy [Hartree]']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for data in sortedData:
                writer.writerow(data)
        print(f'Energy data saved to {csvPath}')