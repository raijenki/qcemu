## This file contains all the quantum gates as well the IO operations
import copy
from functools import reduce
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import random
import sys
from constants import STATES

class QuantumCircuit():
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


        # Since most of code is redundant, define all matrixes here
        self.identitymatrix = np.array([(1, 0), (0, 1)], dtype=np.float)
        self.X = np.array([[0, 1], [1, 0]])

##  __  __ _    _ _   _______ _____       _____       _______ ______ 
## |  \/  | |  | | | |__   __|_   _|     / ____|   /\|__   __|  ____|
## | \  / | |  | | |    | |    | |______| |  __   /  \  | |  | |__   
## | |\/| | |  | | |    | |    | |______| | |_ | / /\ \ | |  |  __|  
## | |  | | |__| | |____| |   _| |_     | |__| |/ ____ \| |  | |____ 
## |_|  |_|\____/|______|_|  |_____|     \_____/_/    \_\_|  |______|
##
##                   (CCNOT, TOFFOLI, SWAP, ARITHMETICS)


    # This code only works on 2-qubits
    def swap(self):
        self.swapMatrix = np.array([(1, 0, 0, 0), (0, 0, 1, 0), (0, 1, 0, 0), (0, 0, 0, 1)])
        self.qubit = self.swapMatrix.dot(self.qubit)

    # This code (for generating CNOT matrixes) was based on 
    # https://github.com/adamisntdead/QuSimPy/blob/master/QuSim.py
    # My modification was to change the qubit count and transform the gateMatrix into array
    def cx(self, control=0, target=1):
        self.identity = np.eye(2)
        self.notmatrix = np.matrix([[0, 1], [1, 0]], dtype = np.float)
        C = np.mat([
                [float('nan'), 0],
                [0, 1]])
        gateOrder = []
        for i in range(0, self.noqubits):
            if (i == control):
                gateOrder.append(C)
            elif (i == target):
                gateOrder.append(self.notmatrix)
            else:
                gateOrder.append(self.identity)
        
        # Generate the gate and then replace the NaNs to Id gates
        newGate = reduce(np.kron, gateOrder)
        n = newGate.shape[0]
        gateMatrix = np.mat([[newGate[i, j] if not np.isnan(newGate[i, j]) else 1 if i == j else 0 for j in range(n)] for i in range(n)])
        gateMatrix = gateMatrix.tolist()
        self.qubit = np.dot(self.qubit, gateMatrix)

    # def cz(self, control=0, target=1):
    #     self.identity = np.eye(2)
    #     self.czmatrix = np.matrix([[1, 0], [0, -1]], dtype = np.float)
    #     C = np.mat([
    #             [float('nan'), 0],
    #             [0, 1]])
    #     gateOrder = []
    #     for i in range(0, self.noqubits):
    #         if (i == control):
    #             gateOrder.append(C)
    #         elif (i == target):
    #             gateOrder.append(self.czmatrix)
    #         else:
    #             gateOrder.append(self.identity)
        
    #     # Generate the gate and then replace the NaNs to Id gates
    #     newGate = reduce(np.kron, gateOrder)
    #     n = newGate.shape[0]
    #     gateMatrix = np.mat([[newGate[i, j] if not np.isnan(newGate[i, j]) else 1 if i == j else 0 for j in range(n)] for i in range(n)])
    #     print(gateMatrix)
    #     gateMatrix = gateMatrix.tolist()
    #     self.qubit = np.dot(self.qubit, gateMatrix)

    def rk(self, index):
        self.rkMatrix = np.array([(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, math.e**((2*math.pi*1.j)/(2**index)))])       
        self.qubit = self.ccnotMatrix.dot(self.qubit)

    # CCNOT gate (TOFFOLI)
    def ccnot(self):
        self.ccnotMatrix = np.array([(1, 0, 0, 0, 0, 0, 0, 0),
                                     (0, 1, 0, 0, 0, 0, 0, 0),
                                     (0, 0, 1, 0, 0, 0, 0, 0),
                                     (0, 0, 0, 1, 0, 0, 0, 0),
                                     (0, 0, 0, 0, 1, 0, 0, 0),
                                     (0, 0, 0, 0, 0, 1, 0, 0),
                                     (0, 0, 0, 0, 0, 0, 0, 1),
                                     (0, 0, 0, 0, 0, 0, 1, 0)])
        self.qubit = self.ccnotMatrix.dot(self.qubit)
    
    def increment(self):
        copied = copy.deepcopy(self.qubit)
        # Just to make things easier
        self.qubit[0] = copied[self.nstates-1]
        for i in range(1, self.nstates):
            self.qubit[i] = copied[i-1]

    def decrement(self):
        copied = copy.deepcopy(self.qubit)
        # Just to make things easier
        self.qubit[self.nstates-1] = copied[0]
        for i in range(0, self.nstates - 2):
            self.qubit[i] = copied[i+1]

##   _____ _____ _   _  _____ _      ______       _____       _______ ______ 
##  / ____|_   _| \ | |/ ____| |    |  ____|     / ____|   /\|__   __|  ____|
## | (___   | | |  \| | |  __| |    | |__ ______| |  __   /  \  | |  | |__   
##  \___ \  | | | . ` | | |_ | |    |  __|______| | |_ | / /\ \ | |  |  __|  
##  ____) |_| |_| |\  | |__| | |____| |____     | |__| |/ ____ \| |  | |____ 
## |_____/|_____|_| \_|\_____|______|______|     \_____/_/    \_\_|  |______|
##                                                                           
##                (HAD, ROTX/NOT, ROTY, ROTZ, PHASE, T, TDG)


# Define
    def apply(self, qbit, matrix):

        self.applymatrix = matrix
        if(qbit == 0):
            self.allcircuit = self.applymatrix
            for i in range(self.noqubits - 1):
                self.allcircuit = np.kron(self.identitymatrix, self.allcircuit)
        else:
            # Identity tensors
            self.allcircuit = self.identitymatrix
            for i in range(qbit - 1):
                self.allcircuit = np.kron(self.identitymatrix, self.allcircuit)
            # Hadamard tensors
            self.allcircuit = np.kron(self.applymatrix, self.allcircuit)
            # Identity tensors
            for i in range(qbit, self.noqubits - 1):
                self.allcircuit = np.kron(self.identitymatrix, self.allcircuit)
        # Final tensor
        self.qubit = self.allcircuit.dot(self.qubit)
       

# Defines the Hadamard function
    def h(self, qbit):
        self.hadamardmatrix = (1/math.sqrt(2))*np.array([(1, 1), (1, -1)], dtype=np.complex)
        self.apply(qbit, self.hadamardmatrix)

# Defines the ROTX (NOT) function
    def rotx(self, qbit):
        self.rotxmatrix = np.array([(0, 1), (1, 0)], dtype = np.float)
        self.apply(qbit, self.rotxmatrix)

# Defines the ROTY function
    def roty(self, qbit):
        self.rotymatrix = np.array([(0, 0 - 1.j), (0 + 1.j, 0)], dtype=np.complex)
        self.apply(qbit, self.rotymatrix)

# Defines the ROTZ function
    def rotz(self, qbit):
        self.rotzmatrix = np.array([(1, 0), (0, - 1)], dtype=np.float)
        self.apply(qbit, self.rotzmatrix)

# Defines the PHASE function
    def p(self, qbit, angle):
        # Define PHASE matrix
        rad = math.radians(angle)
        self.pmatrix = np.array([(1, 0), (0, 0 + math.e**(1.j*rad))], dtype=np.complex)
        self.apply(qbit, self.pmatrix)

# Defines the T-GATE function
    def t(self, qbit):
        # Define T-GATE matrix
        self.tmatrix = np.array([(1, 0), (0, 0 + math.e**(1.j*(math.pi)/4))], dtype=np.complex)
        self.apply(qbit, self.tmatrix)

# Defines the TDG-GATE function
    def tdg(self, qbit):
        # Define T-GATE matrix
        self.tdgmatrix = np.array([(1, 0), (0, 0 + math.e**(-1.j*(math.pi)/4))], dtype=np.complex)
        self.apply(qbit, self.tdgmatrix)

 ##  _____     ______  
 ## |_   _|   / / __ \ 
 ##   | |    / / |  | |
 ##   | |   / /| |  | | (Read, Write, Visualization)
 ##  _| |_ / / | |__| |
 ## |_____/_/   \____/
 ##
    # Read the system (although it uses the random module, running one time might not be random enough)
    def read(self):
        self.possibleoutcome = np.arange(self.nstates)
        self.probqubit = np.square(np.absolute(self.qubit))
        x = random.choices(self.possibleoutcome, weights=self.probqubit)
        print("Read Quantum States: ", x)

    # Read, but reads multiple times
    def readmultiple(self, shots=1000, save=0):
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
        if (save == 0):
            plt.show()
        else:
            plt.savefig("reads.png")

    
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

    def viz2(self, save=0):
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
        if (save == 0):
            plt.show()
        else:
            plt.savefig("bloch.png")

    def show(self):
        print(self.qubit)