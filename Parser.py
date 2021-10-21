# This file contains the parser for a domain specific language of QCEMU
# For more information, refer the pdf file in the github repository

from QuantumCircuit import QuantumCircuit

def readFile(arq):
    qc = 0
    f = open(arq, "r")
    for line in f:
        parsed = line.split(" ")
        # This ensures we can comment in the file
        if parsed[0][0] == "#":
            continue
        if len(parsed) == 1:
            execute(parsed[0], qc)
        if len(parsed) == 2:
            qc = execute(parsed[0], qc, parsed[1])
        if len(parsed) == 3:
            qc = execute(parsed[0], qc, parsed[1], parsed[2])

def execute(cmd, qc, val1=0, val2=0):
    val1 = int(val1)
    val2 = int(val2)

    # START THE CIRCUIT
    if(cmd == 'INIT'):
        qc = QuantumCircuit(val1)
    # WRITE CIRCUIT
    if(cmd == 'WRITE'):
        qc.writebin(str(val1))
    # HADAMARD
    if(cmd == 'HAD'):
        qc.h(val1)
    # NOT (ROTX) GATE
    if(cmd == 'NOT'):
        qc.rotx(val1)
    # ROTY GATE
    if(cmd == 'ROTY'):
        qc.roty(val1)
    # ROTZ GATE
    if(cmd == 'ROTZ'):
        qc.rotz(val1)
    # PHASE GATE
    if(cmd == 'PHASE'):
        qc.p(val1, val2)
    # T GATE
    if(cmd == 'T'):
        qc.t(val1)
    # T-DAGGER GATE
    if(cmd == 'TDG'):
        qc.tdg(val1)
    # CX GATE
    if(cmd == 'CX'):
        qc.cx(val1, val2)
    # READ STATES
    if(cmd == 'READ'):
        qc.readmultiple(1000, 1)
    # BLOCH SPHERES
    if(cmd == 'BLOCH'):
        qc.viz2(1)
    return qc             

