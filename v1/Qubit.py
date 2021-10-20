import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import random

## THE QUBIT CLASS
## This file is meant for unitary gates
## For complex circuits, refer to QuantumCircuit.py

class Qubit():
    def __init__(self):
        self.noqubits = 1
        self.qubit = np.zeros((2), dtype=np.complex)
        self.ampliqubit = np.zeros((2), dtype=np.float)
        self.phase = np.zeros((2), dtype=np.float)

    # Hadamard gate, equals probabilities 
    def had(self):
        self.hadamardmatrix = (1/math.sqrt(2))*np.array([(1, 1), (1, -1)], dtype=np.complex)
        self.allcircuit = self.hadamardmatrix
        for i in range(self.noqubits - 1):
            self.allcircuit = np.kron(self.hadamardmatrix, self.allcircuit)
        self.qubit = self.allcircuit.dot(self.qubit)

        # Not gate, transforms |0> into |1> and vice-versa

    def quarray(self):
        #print(self.qubit)
        return np.array(self.qubit, dtype=np.complex)

    def read(self):
        self.possibleoutcome = np.arange(2)
        self.probqubit = np.square(np.absolute(self.qubit))
        x = random.choices(self.possibleoutcome, weights=self.probqubit)
        print("Read Quantum States: ", x)
        return x

    def readmultiple(self, nshots=1000):
        self.possibleoutcome = np.arange(2)
        self.probqubit = np.square(np.absolute(self.qubit))
        self.measurementstate = np.zeros((2), dtype=np.int)

        for i in range(nshots):
            x = random.choices(self.possibleoutcome, weights=self.probqubit)
            self.measurementstate[x] = self.measurementstate[x] + 1

        # Plotting shots results in matplotlib
        plt.grid(b=True)
        plt.bar(self.possibleoutcome, (self.measurementstate/nshots))
        plt.xlabel("Quantum States")
        plt.ylabel("Probability")
        plt.xticks(self.possibleoutcome, rotation='65')
        plt.show()

    def init(self, state):
        if state == 0:
            self.qubit[0] = 1
            self.qubit[1] = 0
        if state == 1:
            self.qubit[0] = 0
            self.qubit[1] = 1

    
    def identity(self):
        self.imatrix = np.array([(1, 0), (0,1)], dtype=np.float)
        self.allcircuit = self.imatrix
        for i in range(self.noqubits - 1):
            self.allcircuit = np.kron(self.imatrix, self.allcircuit)
        self.qubit = self.allcircuit.dot(self.qubit)

    def rotx(self):
        self.rotxmatrix = np.array([(0,1), (1, 0)], dtype=np.complex)
        self.allcircuit = self.rotxmatrix
        for i in range(self.noqubits - 1):
            self.allcircuit = np.kron(self.rotxmatrix, self.allcircuit)
        self.qubit = self.allcircuit.dot(self.qubit)

    def roty(self):
        self.rotymatrix = np.array([(0, 0 - 1.j), (0 + 1.j, 0)], dtype=np.complex)
        self.allcircuit = self.rotymatrix
        for i in range(self.noqubits - 1):
            self.allcircuit = np.kron(self.rotymatrix, self.allcircuit)
        self.qubit = self.allcircuit.dot(self.qubit)

    def rotz(self):
        self.rotzmatrix = np.array([(1, 0), (0, - 1)], dtype=np.float)
        self.allcircuit = self.rotzmatrix
        for i in range(self.noqubits - 1):
            self.allcircuit = np.kron(self.rotzmatrix, self.allcircuit)
        self.qubit = self.allcircuit.dot(self.qubit)

    def p(self, angle):
        rad = math.radians(angle)
        self.pmatrix = np.array([(1, 0), (0, 0 + math.e**(1.j*rad))], dtype=np.complex)
        self.allcircuit = self.pmatrix
        for i in range(self.noqubits - 1):
            self.allcircuit = np.kron(self.pmatrix, self.allcircuit)
        self.qubit = self.allcircuit.dot(self.qubit)

    def t(self):
        self.tmatrix = np.array([(1, 0), (0, 0 + math.e**(1.j*(math.pi)/4))], dtype=np.complex)
        self.allcircuit = self.tmatrix
        for i in range(self.noqubits - 1):
            self.allcircuit = np.kron(self.tmatrix, self.allcircuit)
        self.qubit = self.allcircuit.dot(self.qubit)

    def tdg(self):
        self.tdgmatrix = np.array([(1, 0), (0, 0 + math.e**(-1.j*(math.pi)/4))], dtype=np.complex)
        self.allcircuit = self.tdgmatrix
        for i in range(self.noqubits - 1):
            self.allcircuit = np.kron(self.tdgmatrix, self.allcircuit)
        self.qubit = self.allcircuit.dot(self.qubit)

    def writeraw(self, pos, val):
        self.qubit[pos] = val 

    def readraw(self, pos):
        return self.qubit[pos]

    # Sanity check
    def test(self):
        print("Hello!")