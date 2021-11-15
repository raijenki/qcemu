# qcemu
## What is this?
This is a Quantum Computer simulator written in Python 3 for the course `FDD3280 - Quantum Computing to Computer Scientists` at KTH Royal Institute of Technology. Part of the code used was shown in class.

## What does this implements?
* Quantum Instructions: HAD, ROTX (NOT), ROTY, PHASE, CNOT, SWAP and CCNOT (TOFFOLI);
* Quantum Arithmetic operators (INCREMENT and DECREMENT);
* Read and write operations, alongside output graphs;
* A domain-specific language (QCLang) for usage.

## How do I use QCEmu/QCLang?
Please refer to `report.pdf` file for both QCEmu usage and the QCLang grammar.

## Are there bugs?
There are bugs in the CNOT for over 3 qubits and the CCNOT (which I spent few days and couldn't figure why it was not working).

## Should I actually use this?
No.
