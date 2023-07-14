# SimuQ

SimuQ is a domain-specific language for quantum simulation with analog compilation.

## Requirements

Basic environment requirements for SimuQ:

* Python 3.9

* NetworkX 2.8.7

* NumPy 1.23.4

* SciPy 1.9.3

The following are optional, based on what part of SimuQ you want to use:

* Qiskit 0.43.0

* Qiskit-Terra 0.24.1

* Julia 1.7.2

* Bloqade 0.1.13

* Amazon Braket SDK 1.39.1

## Project Structure

`examples/`: Example notebooks of running SimuQ and obtaining results.

`aais/`: AAIS of many backends in AAPSI Specification Language.

`backends/`: The translation to API languages of different machine backends.

`simuq/`: The core compiler and language implementation of SimuQ.

`systems/`: Hamiltonian evolution implemented in Hamiltonian Modeling Language.

`systems/benchmark/`: A small benchmark of quantum Hamiltonian simulation as reported in the paper.



