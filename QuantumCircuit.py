import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import random

from numpy.core.fromnumeric import swapaxes
from Qubit import Qubit
import sys

class QuantumCircuit():
    def __init__(self, noqubits=1):
        self.noqubits = noqubits
        self.nstates = 2**noqubits
        # Print some useful message
        print("Initiation of allocation")
        print("Number of qubits = ", self.noqubits)
        print("Number of possible states = ", self.nstates)

        # Create the desired number of qubits objects
        self.objs = [Qubit() for i in range(noqubits)]
        
        # Set matrices as the size of number of states
        self.qubit = np.zeros((self.nstates), dtype=np.complex)
        self.ampliqubit = np.zeros((self.nstates), dtype=np.float)
        self.phase = np.zeros((self.nstates), dtype=np.float)

        # Basics matrixp
        self.identitymatrix = np.array([(1, 0), (0, 1)], dtype=np.float)

    # Get the entire state matrix
    def circuitMatrix(self):
        if(self.noqubits == 1):
            self.qubit = self.objs[0].quarray()
            return 0
        for i in range(0, self.noqubits - 1, 1):
            a = self.objs[i].quarray()
            b = self.objs[i + 1].quarray()
            if(i == 0):
                aux = np.kron(b, a)
            else:
                aux = np.kron(aux, b)
        self.qubit = aux.reshape(4, 1)
        print(self.qubit)
        return 0

    def readmultiple(self, shots=1000):
        self.circuitMatrix()
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
    
    def viz2(self):
        self.circuitMatrix()
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


    ###########################################################
    #### THESE FUNCTIONS ARE RELATED TO TWO-GATES          ####
    ###########################################################


    def cx(self, tgt, control):
        if ((tgt or control) >= self.noqubits):
            print("Invalid value!")
            sys.exit(1)
        else:
            # swap
            tgtbit = self.objs[tgt].quarray()
            ctrlbit = self.objs[control].quarray()
            self.objs[tgt].writeraw(1, ctrlbit[1])
            self.objs[control].writeraw(1, tgtbit[1])


    ###########################################################
    #### THESE FUNCTIONS ARE FOR ACCESSING THE QUBIT CLASS ####
    ###########################################################

    # Call hadamard gate
    def had(self, numb):
        if(numb > self.noqubits - 1):
            print("Cannot read this qubit!")
            sys.exit(1)
        self.objs[numb].had()
        print("HAD gate applied on qubit ", numb)

    # Call not function
    def rotx(self, numb):
        if(numb > self.noqubits - 1):
            print("Cannot read this qubit!")
            sys.exit(1)
        self.objs[numb].rotx()
        print("ROTX (NOT) gate applied on qubit ", numb)

    # Call ROTY function
    def roty(self, numb):
        if(numb > self.noqubits - 1):
            print("Cannot read this qubit!")
            sys.exit(1)
        self.objs[numb].roty()
        print("ROTY gate applied on qubit ", numb)

    # Call ROTZ function
    def rotz(self, numb):
        if(numb > self.noqubits - 1):
            print("Cannot read this qubit!")
            sys.exit(1)
        self.objs[numb].rotz()
        print("ROTZ gate applied on qubit ", numb)

    # Call phase gate
    def p(self, numb, angle):
        if(numb > self.noqubits - 1):
            print("Cannot read this qubit!")
            sys.exit(1)
        self.objs[numb].p(angle)
        print("PHASE gate applied on qubit ", numb)

    # Read a single qubit 
    def read_qubit(self, numb):
        if(numb > self.noqubits - 1):
            print("Cannot read this qubit!")
            sys.exit(1)
        self.objs[numb].read()

    # Read a single qubit, but reads multiple times
    def readmultiple_qubit(self, numb, nshots=1000):
        if(numb > self.noqubits - 1):
            print("Cannot read this qubit!")
            sys.exit(1)
        self.objs[numb].readmultiple(nshots)       
    
    # Write into a single qubit
    def write_qubit(self, numb, initstate):
        self.objs[numb].write(initstate)

    # Write into a single qubit, but in binary
    def writebin_qubit(self, numb, initstate):
        self.objs[numb].writebin(initstate)

    def quarray(self, numb):
        print(self.objs[numb].quarray())

    # Sanity check
    def test(self, numb):
        return self.objs[numb].quarray()