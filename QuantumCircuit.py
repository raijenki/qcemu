## This file contains all the quantum gates as well the IO operations
import copy
from functools import reduce
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import random
import sys
import time

class QuantumCircuit():
    def __init__(self, noqubits=1):

        # This will ensure that our system will function properly, as CNOT is designed for only
        # 2 or 3 qubits
        if(noqubits > 3):
            print("This number of qubits is not allowed!")
            print("Leaving the program...")
            sys.exit(1)

        # Don't even need to put an else
        self.noqubits = noqubits
        self.nstates = 2**noqubits ## Number of possible states is equal to 2^(number of qubits)
        self.figimage = 1

        # Print some useful messages
        print("Initiation of allocation")
        print("Number of qubits = ", self.noqubits)
        print("Number of possible states = ", self.nstates)

        # Set matrices as the size of number of states
        self.qubit = np.zeros((self.nstates), dtype=np.complex)
        self.ampliqubit = np.zeros((self.nstates), dtype=np.float)
        self.phase = np.zeros((self.nstates), dtype=np.float)

##  __  __ _    _ _   _______ _____       _____       _______ ______ 
## |  \/  | |  | | | |__   __|_   _|     / ____|   /\|__   __|  ____|
## | \  / | |  | | |    | |    | |______| |  __   /  \  | |  | |__   
## | |\/| | |  | | |    | |    | |______| | |_ | / /\ \ | |  |  __|  
## | |  | | |__| | |____| |   _| |_     | |__| |/ ____ \| |  | |____ 
## |_|  |_|\____/|______|_|  |_____|     \_____/_/    \_\_|  |______|
##
##                   (CCNOT, TOFFOLI, ARITHMETICS)


# Took this idea from https://github.com/corbett/QuantumComputing
# Computing matrices for cnot gates with >2 qubits is something non-trivial, 
# so the idea was simply to hardcode into the project.
# It was adapted to use in this code (i.e., transform a matrix into np.array and IFs).
# Comments are also mine.
# It also appears that target and control are inverted there, but I won't be changing the variable names. lol
    def cx(self, control=0, target=1):
        H = 1./math.sqrt(2)*np.matrix('1 1; 1 -1') # Hadamard
        CNOT2_01=np.matrix('1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0') 
        CNOT2_10 = np.kron(H,H)*CNOT2_01*np.kron(H,H)  # This uses the hadamard trick to invert function
        eye=np.eye(2,2)

        # For two entangled qubits
        if(self.noqubits==2):
            if(control==0 and target==1):
                self.cnot=CNOT2_01
                #print("CNOT02_01")
            if(control==1 and target==0):
                self.cnot=CNOT2_10
               # print("CNOT02_10")
        # For three entangled qubits
        if(self.noqubits == 3):
            if(control==0 and target==1):
                self.cnot=np.kron(CNOT2_01,eye)
            if(control==1 and target==0):
                self.cnot=np.kron(CNOT2_10,eye)
            if(control==1 and target==2):
                self.cnot=np.kron(eye,CNOT2_01)
            if(control==2 and target==1):
                self.cnot=np.kron(eye,CNOT2_10)
            if(control==0 and target==2):
                self.cnot=np.matrix('1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 1; 0 0 0 0 0 0 1 0')
            if(control==2 and target==0):
                self.cnot=np.matrix('1 0 0 0 0 0 0 0; 0 0 0 0 0 1 0 0; 0 0 1 0 0 0 0 0; 0 0 0 0 0 0 0 1; 0 0 0 0 1 0 0 0; 0 1 0 0 0 0 0 0; 0 0 0 0 0 0 1 0; 0 0 0 1 0 0 0 0')

        self.cnot = np.array(self.cnot.tolist())
        self.qubit = self.cnot.dot(self.qubit)

    # This allows the implementation of R-k gate
    # so we can attempt the QFT later
    def rk(self, r_index):
        r_index = math.e**((2*math.pi*1.j)/(2*math.r_index))
        rk=np.matrix('1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 1 0')


    # CCNOT gate (TOFFOLI)
    # By definition, target is qubit 2, controls are 0 and 1
    def ccx(self, control1=0, control2=1, target=2):
        # Ensure we have three qubits
        if (self.noqubits < 3):
            print("This operation can't be done with less than 3 qubits!")
            sys.exit(1)
        else:
            # Toffoli matrix
            self.toffoli = np.array([(1,0,0,0,0,0,0,0), 
                                (0,1,0,0,0,0,0,0), 
                                (0,0,1,0,0,0,0,0), 
                                (0,0,0,1,0,0,0,0), 
                                (0,0,0,0,1,0,0,0), 
                                (0,0,0,0,0,1,0,0), 
                                (0,0,0,0,0,0,0,1), 
                                (0,0,0,0,0,0,1,0)], dtype=np.float)
            self.qubit = np.dot(self.toffoli, self.qubit)
            # # This produces Toffoli from CNOTs, also not working god-knows-why
            # self.h(target)
            # self.cx(target, control1)
            # #self.cx(control1, target)
            # self.tdg(target)
            # self.cx(target, control2)
            # #self.cx(control2, target)
            # self.t(target)
            # self.cx(target, control1)
            # #self.cx(control1, target)
            # self.tdg(target)
            # self.cx(target, control2)
            # #self.cx(control2, target)
            # self.t(target)
            # self.cx(control1, control2)
            # #self.cx(control2, control1)
            # self.h(target)
            # self.t(control2)
            # self.tdg(control1)
            # self.cx(control1, control2)
            # self.cx(control2, control1)
    
    # Increment function
    # Transport the value of the STATE to STATE+1
    def increment(self):
        copied = copy.deepcopy(self.qubit)
        # Handle exception of last state to the first
        self.qubit[0] = copied[self.nstates-1]
        # Loop
        for i in range(1, self.nstates):
            self.qubit[i] = copied[i-1]

    # Decrement function
    # Transport the value of the STATE to STATE-1
    def decrement(self):
        copied = copy.deepcopy(self.qubit)
        # Handle exception of last state to the first
        self.qubit[self.nstates-1] = copied[0]
        # Loop
        for i in range(0, self.nstates - 1 ):
            self.qubit[i] = copied[i+1]

##   _____ _____ _   _  _____ _      ______       _____       _______ ______ 
##  / ____|_   _| \ | |/ ____| |    |  ____|     / ____|   /\|__   __|  ____|
## | (___   | | |  \| | |  __| |    | |__ ______| |  __   /  \  | |  | |__   
##  \___ \  | | | . ` | | |_ | |    |  __|______| | |_ | / /\ \ | |  |  __|  
##  ____) |_| |_| |\  | |__| | |____| |____     | |__| |/ ____ \| |  | |____ 
## |_____/|_____|_| \_|\_____|______|______|     \_____/_/    \_\_|  |______|
##                                                                           
##                (HAD, ROTX/NOT, ROTY, ROTZ, PHASE, T, TDG)


# This function is the general function to apply the unitary gates
# This avoids repeated code. Also, I noticed that this reusage is
# a common practice among quantum simulator codes.

    def apply(self, qbit, matrix):
        self.applymatrix = matrix
        self.identitymatrix = np.array([(1, 0), (0, 1)], dtype=np.float)
        # start = time.time()
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
        # end = time.time()
        print(end - start)

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

# Defines the PHASE function (depends on the angle)
    def p(self, qbit, angle):
        rad = math.radians(angle)
        self.pmatrix = np.array([(1, 0), (0, 0 + math.e**(rad*1.j))], dtype=np.complex)
        self.apply(qbit, self.pmatrix)

# Defines the T-GATE function
    def t(self, qbit):
        self.tmatrix = np.array([(1, 0), (0, 0 + math.e**(1.j*(math.pi)/4))], dtype=np.complex)
        self.apply(qbit, self.tmatrix)

# Defines the TDG-GATE function
    def tdg(self, qbit):
        self.tdgmatrix = np.array([(1, 0), (0, 0 + math.e**(-1.j*(math.pi)/4))], dtype=np.complex)
        self.apply(qbit, self.tdgmatrix)

 ##  _____     ______  
 ## |_   _|   / / __ \ 
 ##   | |    / / |  | |
 ##   | |   / /| |  | | (Read, Write, Visualization)
 ##  _| |_ / / | |__| |
 ## |_____/_/   \____/
 ##
    # Pretty much copied from Stefano's class, except for the save part
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
        # This will allow for saving multiple times
        else:
            fname = "reads" + str(self.figimage) + ".png"
            plt.savefig(fname)

    # This has been useless for me
    # def write(self, initstate):
    #     self.initstate = initstate
    #     print("Written state = ", self.initstate)

    #     if(self.initstate > self.nstates-1):
    #         print("Initial state can't be represented in the system")
    #         sys.exit(1)
    #     self.qubit[self.initstate] = 1

    # Write state in binary
    def writebin(self, binarystring):
        self.initstate = int(binarystring, 2)
        if (self.initstate > self.nstates-1):
            print("Initial state can't be represented in the system")
            sys.exit(1)
        self.qubit[self.initstate] = 1.0

    # Visualization in Bloch circles (not spheres!)
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
        # This will allow for saving multiple times
        else:
            fname = "bloch" + str(self.figimage) + ".png"
            plt.savefig(fname)

    # Helper function to show state vector and allows copy to variable
    def show(self):
        print(self.qubit)
        return self.qubit

