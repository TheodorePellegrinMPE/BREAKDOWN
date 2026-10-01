import subprocess

class IsomerOptimizer:
    def __init__(self, folder, isomerName):
        """Initialize the IsomerOptimizer with the folder and set file paths.
        """
        self.executable = f'{folder}/deMon.static-website.x'
        self.coordinatesFile = f'{folder}/{isomerName}.xyz'
        self.inputFile = f'{folder}/deMon.inp'
        self.outputFile = f'{folder}/deMon.out'
        self.isomerName = isomerName
        self.visitedStates = []
        self.maxRepetition = 20
        with open(f'{folder}/visited_states.txt', 'r') as file:
            self.visitedStates = [line.strip() for line in file.readlines()]

    def optimizeIsomerGeometry(self):
        """Optimize isomer geometry using the deMon program.
        """
        if self.isomerName in self.visitedStates:
            self.createInputFile()
            attemptsCount = 0
            while attemptsCount < self.maxRepetition:
                self.runDemon()
                if self.checkConvergence():
                    print("Optimization converged!")
                    break
                else:
                    print('Optimization did not converge. Retrying...')
                    optimizedCoordinates = self.readCoordinates()
                    self.writeCoordinates(optimizedCoordinates)
                    attemptsCount = attemptsCount + 1
            if attemptsCount >= self.maxRepetition:
                print('Optimization failed!')
    
    def createInputFile(self):
        """Create the initial input file using coordinates from the XYZ file.
        """
        initialCoordinates = self.readCoordinatesFromXYZ()
        self.writeInitialInput(initialCoordinates)
    
    def readCoordinatesFromXYZ(self):
        """Read coordinates from the XYZ file starting from the third line.
        """
        with open(self.coordinatesFile, 'r') as file:
            lines = file.readlines()[2:]
        coordinates = [line.split()[0:4] for line in lines]
        return coordinates

    def writeInitialInput(self, coordinates):
        """Write the initial input file with provided header and coordinates.
        """
        with open(self.inputFile, 'w') as file:
            file.write('# Phen isomer optimisation\n')
            file.write('OPTIMIZATION\n')
            file.write('DFTB SCC FERMI=1500 DISP=2\n')
            file.write('PARAM PTYPE=BIO\n')
            file.write('~/basis\n')
            file.write('MULTIPLICITY 2\n')
            file.write('CHARGE 1\n')
            file.write('GEOMETRY\n')
            for coord in coordinates:
                file.write(f'{coord[0]:<2} {coord[1]:>15} {coord[2]:>15} {coord[3]:>15}\n')

    def readCoordinates(self):
        """Read optimized coordinates from the output file.
        """
        with open(self.outputFile, 'r') as file:
            lines = file.readlines()
        startIndex = lines.index(' OPTIMIZED STRUCTURE IN ANGSTROM (INPUT ORDER)\n')
        startIndex += 4
        coordinates = []
        for line in lines[startIndex:]:
            if not line.strip():
                break
            coordinates.append(line.split()[3:6])
        return coordinates

    def writeCoordinates(self, coordinates):
        """Write optimized coordinates to the input file.
        """
        with open(self.inputFile, 'r') as file:
            lines = file.readlines()
        startWrite = lines.index('GEOMETRY\n') + 1
        for i, coord in enumerate(coordinates):
            elementName = lines[startWrite + i].split()[0]
            lines[startWrite + i] = f'{elementName:>2} {coord[0]:>5} {coord[1]:>5} {coord[2]:>5}\n'
        with open(self.inputFile, 'w') as file:
            file.writelines(lines)

    def runDemon(self):
        """Run the deMon program.
        """
        subprocess.Popen([self.executable]).wait()

    def checkConvergence(self):
        """Check if the optimization has converged by examining the output file.
        """
        with open(self.outputFile, 'r') as file:
            return 'not converged' not in file.read()
