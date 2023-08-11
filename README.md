# SimuQ

**SimuQ** is a domain-specific language for quantum simulation with analog compilation.

Our [project website](https://pickspeng.github.io/SimuQ/) includes use cases and documentation of SimuQ.

We illustrate our design and benefits in our [arXiv paper](https://arxiv.org/abs/2303.02775).

## Installation

We recommend Python 3.9 or greater since many optional dependencies have requirements for Python version.

We encourage using ``pip`` for SimuQ installation. To install the core components of SimuQ, run the following command in shell:

(Currently we haven't uploaded the repo to PyPI, so please download/clone the repository, enter the (unzipped) folder, and follow the following instructions but substituting ``.`` for ``simuq``, i.e., ``pip install .[all]``.)

```bash
pip install simuq
```

Multiple platforms are supported by SimuQ through different APIs and you may optionally install them based on your needs. You may install all optional dependencies, though not recommended, by

```bash
pip install simuq[all]
```

### Amazon Braket

SimuQ supports compilation to IonQ's trapped-ion devices and QuEra's neutral atom devices through Amazon Braket.

To enable these backends, you may install the dependencies via

```bash
pip install simuq[braket]
```

Then you need to configure AWS IAM roles, regions, and AWS credentials. Please refer to [Amazon Braket Python SDK](https://github.com/aws/amazon-braket-sdk-python).

For examples using Amazon Braket providers of SimuQ, please refer to [Simulating Ising model on QuEra devices](https://github.com/PicksPeng/SimuQ/blob/main/notebooks/ising_quera.ipynb).

### IonQ Quantum Cloud

SimuQ supports compilation to IonQ's trapped-ion devices through IonQ Quantum Cloud.

To enable these backends, you may install the dependencies via

```bash
pip install simuq[ionq]
```

To run through IonQ Quantum Cloud, you must obtain an API key from [IonQ](https://ionq.com/quantum-cloud).

When creating IonQ providers in SimuQ, you must provide the API key either through a string or a file storing the key. For further details please refer to [Simulating Ising model on IonQ devices](https://github.com/PicksPeng/SimuQ/blob/main/notebooks/ising_ionq.ipynb).

### Qiskit

SimuQ supports compilation to IBM's superconducting devices through Qiskit and IBM Quantum.

To enable these backends, you may install the dependencies via

```bash
pip install simuq[qiskit]
```

To run through IBM Quantum, you must obtain an API token from [IBM](https://quantum-computing.ibm.com/).

The examples are TODO.

### dReal solver

The default solver of SimuQ is based on [SciPy](https://scipy.org/)'s least square solver. Other SMT solvers are also supported in SimuQ compiler, like dReal.

To try dReal solver, please refer to [installation of dReal4](https://github.com/dreal/dreal4), and enable the solver by

```bash
pip install simuq[dreal]
```

## Project Structure

`notebooks/`: Example notebooks of running SimuQ and obtaining results.

`simuq/`: The core compiler and language implementation of SimuQ.

`simuq/aais/`: AAIS of many backends in AAIS Specification Language.

`simuq/backends/`: The translation to API languages of different machine backends.

`simuq/systems/`: Hamiltonian evolution implemented in Hamiltonian Modeling Language.

`simuq/systems/benchmark/`: A small benchmark of quantum Hamiltonian simulation as reported in our arXiv paper.


