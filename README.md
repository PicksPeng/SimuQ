# SimuQ

SimuQ is a domain-specific language for quantum simulation with analog compilation.

## Requirements

Basic environment requirements for SimuQ:

* Python 3.9

* NetworkX 2.8.7

* NumPy 1.23.4

* SciPy 1.9.3

To generate IBM machine AAIS for specific machines, programs in Qiskit Pulses or circuits in Qiskit, SimuQ also needs:

* Qiskit 0.38.0

* Qiskit-Terra 0.21.2

To run the generated Bloqade program, you will need:

* Julia 1.7.2

* Bloqade 0.1.13


## Demo

The file `demo.py` contains basic settings to generate a quantum system `qs` in HML, an AAIS `mach` in AAIS Specification Language, and sends the compilation intermediate results to corresponding backend. Running `python demo.py` produces the averaged time consumption of compilation.

The compilation results are stored in variables corresponding to the backends. Quantum circuits are stored in 'circ', Bloqade code is stored in 'code', and IBM pulse schedules are stored in 'program.pulse_schedule'.

The settings are mainly:

* `qsys`: the target system. Options: `MIS`, `Heis`, and `QAOA`, representing respectively the maximal-independent-set Hamitlonian evolution, the Heisenberg model, and the max-cut Hamiltonian evolution.

* `mode`: the machine backend. Options: `Rydberg1d Bloqade`, `IBM Qiskit Pulse`, `Circuit`, `Rydberg2d Bloqade`, and `Rydberg QuEra`, representing respectively 1d Rydberg atom arrays with local controls, IBM's transmon qubit system, quantum circuits supporting rotation along XX, YY, and ZZ, 2d Rydberg atom arrays with local controls, 2d Rydberg atom arrays without local controls (and transpilation to AWS's API for QuEra machines).

* `N` and `K`: the size of the target system. For detailed explanation please check the referred files in `aais/`.

* `Repetition`: the number of repetition that the compilation runs over. The time consumption is averaged over them.

* `D`: the discretization number for continuous systems.

* `Trot`: the Trotterization number in conflict resolve.


## Project Structure

`aais/`: AAIS of many backends in AAPSI Specification Language.

`backends/`: The transpilation to API languages of different machine backends.

`scripts/`: Utility scripts.

`simuq/`: The core compiler and language implementation of SimuQ.

`systems/`: Hamiltonian evolution implemented in Hamiltonian Modeling Language.

`demo.py`: A demo file.



