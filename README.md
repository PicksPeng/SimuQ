# Requirements

Tested SimuQ environment :

Python 3.9

NetworkX 2.8.7

NumPy 1.23.4

SciPy 1.9.3

Qiskit 0.39.0

Qiskit-Terra 0.22.0



# Project Structure

simuq/: The core compiler and language implementation of SimuQ.

systems/: Hamiltonian evolution implemented in Hamiltonian Modeling Language.

aais/: AAIS of many backends in AAPSI Specification Language.

backends/: The transpilation to API languages of different machine backends.

Analog_Hamiltonian_Simulation/: Necessary components dealing with pulses on IBM Qiskit Pulse.

demo.py: A demo file.



# Demo

demo.py contains basic settings to generate a quantum system in HML, an AAIS in AAIS Specification Language, and send the compilation intermediate results to corresponding backend. Running 'python demo.py' produces the averaged time consumption of compilation. The compilation results are stored in variables corresponding to the backends. Quantum circuits are stored in 'circ', Bloqade code is stored in 'code', and IBM pulse schedules are stored in 'program.pulse_schedule'.

qsys represents the evolution to simulate, with options MIS, Heis, and QAOA.

mode represents the machine backends, with options Circuit, Rydberg1d Bloqade, Rydberg2d Bloqade, and IBM Qiskit Pulse 

N and K represent the size of the quantum system.

Repetition represents the number of repetition the compilation runs over. The time consumption is averaged over these repetitions.

D represents the discretization number for continuous systems.

