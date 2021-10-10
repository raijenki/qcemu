# Imports
import math
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import random 
import sys

# Main class which contain all gates and write/read functions
class QuantumSystem:
    def __init__(self, noqubits=1):
        self.noqubits = noqubits
        self.nstates = 2**noqubits ## Number of possible states is equal to 2^(number of qubits)

        # Print some useful message
        print("Initiation of allocation")
        print("Number of qubits = ", self.noqubits)
        print("Number of possible states = ", self.nstates)

        # Set matrices as the size of number of states
        self.qubit = np.zeros((self.nstates), dtype=np.complex)
        self.ampliqubit = np.zeros((self.nstates), dtype=np.float)
        self.phase = np.zeros((self.nstates), dtype=np.float)

        # Basics matrix
        self.identitymatrix = np.array([(1, 0), (0, 1)], dtype=np.float)
        self.cNotmatrix = np.array([(1, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 0), (1, 0, 0, 0)], dtype=np.float)
    
    # Hadamard gate, equals probabilities 
    def had(self):
        self.hadamardmatrix = (1/math.sqrt(2))*np.array([(1, 1), (1, -1)], dtype=np.float)
        self.allcircuit = self.hadamardmatrix

        for i in range(self.noqubits -1):
            self.allcircuit = np.kron(self.hadamardmatrix, self.allcircuit)
        self.qubit = self.allcircuit.dot(self.qubit)
    
    # Not gate, transforms |0> into |1> and vice-versa
    def Not(self):
        self.notmatrix = np.array([(0, 1), (1,0)], dtype=np.float)
        self.allcircuit = self.notmatrix

        for i in range(self.noqubits - 1):
            self.allcircuit = np.kron(self.notmatrix, self.allcircuit)
        self.qubit = self.allcircuit.dot(self.qubit)

    # Read the system (although it uses the random module, running one time might not be random enough)
    def read(self):
        self.possibleoutcome = np.arange(self.nstates)
        self.probqubit = np.square(np.absolute(self.qubit))
        x = random.choices(self.possibleoutcome, weights=self.probqubit)
        print("Read Quantum States: ", x)

    # Read, but reads multiple times
    def readmultiple(self, shots=1000):
        self.possibleoutcome = np.arange(self.nstates)
        self.probqubit = np.square(np.absolute(self.qubit))
        self.measurementstate = np.zeros((self.nstates), dtype=np.int)

        for i in range(shots):
            x = random.choices(self.possibleoutcome, weights=self.probqubit)
            self.measurementstate[x] = self.measurementstate[x] + 1

        # Plotting shots results in matplotlib
        plt.grid(b=True)
        plt.bar(self.possibleoutcome, (self.measurementstate/shots))
        plt.xlabel("Quantum States")
        plt.ylabel("Probability")
        plt.xticks(self.possibleoutcome, rotation='65')
        plt.show()
    
    def write(self, initstate):
        self.initstate = initstate
        print("Written state = ", self.initstate)

        if(self.initstate > self.nstates-1):
            print("Initial state can't be represented in the system")
            sys.exit(1)
        self.qubit[self.initstate] = 1

    def writebin(self, binarystring):
        self.initstate = int(binarystring, 2)
        if (self.initstate > self.nstates-1):
            print("Initial state can't be represented in the system")
            sys.exit(1)
        self.qubit[self.initstate] = 1.0

    def viz2(self):
        self.probqubit = np.absolute(self.qubit)
        self.phasequbit = np.angle(self.qubit)

        rows = int(math.ceil(self.nstates / 8.0))
        cols = min(self.nstates, 8)
        fig, axs = plt.subplots(rows, cols)
        
        for col in range(cols):
            circleExt = matplotlib.patches.Circle((0.5, 0.5), 0.5, color='gray', alpha=0.1)
            circleInt = matplotlib.patches.Circle((0.5, 0.5), self.probqubit[col]/2, color='b', alpha=0.3)
            axs[col].add_patch(circleExt)
            axs[col].add_patch(circleInt)
            axs[col].set_aspect('equal')
            statenumber = "|" + str(col) +  ">"
            axs[col].set_title(statenumber)
            xl = [0.5, 0.5 + 0.5*self.probqubit[col]*math.cos(self.phasequbit[col] + np.pi/2)]
            yl = [0.5, 0.5 + 0.5*self.probqubit[col]*math.sin(self.phasequbit[col] + np.pi/2)]
            axs[col].plot(xl,yl,'r')
            axs[col].axis('off')
        plt.show()