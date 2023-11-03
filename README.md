# SimuQ

[![Latest Version](https://img.shields.io/pypi/v/simuq.svg)](https://pypi.python.org/pypi/simuq)
[![Supported Python Versions](https://img.shields.io/pypi/pyversions/simuq.svg)](https://pypi.python.org/pypi/simuq)

**SimuQ** is a domain-specific language for quantum simulation with analog compilation.

Our [project website](https://pickspeng.github.io/SimuQ/) includes use cases and documentation of SimuQ.

We illustrate our design and benefits in our [arXiv paper](https://arxiv.org/abs/2303.02775).

## Installation

We recommend Python 3.9 or greater since many optional dependencies have a minimum Python version requirement.

We encourage using `pip` for SimuQ installation. To install the core components of SimuQ, run the following command in shell:

```bash
pip install simuq
```

Multiple platforms are supported by SimuQ through different APIs and you may optionally install them based on your needs. You may install all optional dependencies, though not recommended, by

```bash
pip install "simuq[all]"
```

### Amazon Braket

SimuQ supports compilation to IonQ's trapped-ion devices and QuEra's neutral atom devices through Amazon Braket. Install the Amazon Braket provider and its dependencies by running

```bash
pip install "simuq[braket]"
```

If running on QPUs, make sure that your AWS account is onboarded to Amazon Braket, as per the instructions [here](https://github.com/amazon-braket/amazon-braket-sdk-python#prerequisites) (this isn't necessary for running on the local simulator).

### IonQ Quantum Cloud

SimuQ supports compilation to IonQ's trapped-ion devices through IonQ Quantum Cloud. Install the IonQ provider and its dependencies by running

```bash
pip install "simuq[ionq]"
```

To run through IonQ Quantum Cloud, you must obtain an API key from [IonQ](https://ionq.com/quantum-cloud). When creating IonQ providers in SimuQ, you must provide the API key either through a string or a file storing the key.

### Qiskit

SimuQ supports compilation to IBM's superconducting devices through Qiskit and IBM Quantum. Install the Qiskit provider and its dependencies by running

```bash
pip install "simuq[qiskit]"
```

To run through IBM Quantum, you must obtain an API token from [IBM](https://quantum-computing.ibm.com/).

### QuTiP

You can simulate the dynamics of a SimuQ Hamiltonian with QuTiP. Install the QuTiP provider and its dependencies by running

```bash
pip install "simuq[qutip]"
```

### Installing from source

You can also install from source by cloning this repository and running a pip install command in the root directory of the repository:

```bash
git clone git@github.com:PicksPeng/SimuQ.git
cd SimuQ
pip install .
```

If installing extras, remember to include quotes around the target (for example, `pip install ".[dev]"`). To install in editable mode, run `pip install` with the `-e` flag: `pip install -e ".[dev]"`.

## Examples

For examples of SimuQ usage, refer to the notebooks in the [tutorials folder](https://github.com/PicksPeng/SimuQ/tree/main/notebooks/tutorials).

## Project Structure

`notebooks/`: Example notebooks of running SimuQ and obtaining results.

`notebooks/tutorials`: Tutorial notebooks on usage of SimuQ.

`notebooks/artifact_evaluation`: AE notebooks to reproduce the data in our paper.

`src/simuq/`: The core compiler and language implementation of SimuQ.

`src/simuq/aais/`: AAIS of many backends in AAIS Specification Language.

`src/simuq/backends/`: The translation to API languages of different machine backends.

`src/simuq/systems/`: Hamiltonian systems implemented in Hamiltonian Modeling Language.

`src/simuq/systems/benchmark/`: A small benchmark of quantum Hamiltonian simulation as reported in our paper.

## Citations

If you use SimuQ in your work, please cite our paper.

```
@article{peng2023simuq,
  title    = {SimuQ: A Domain-Specific Language For Quantum Simulation With Analog Compilation},
  author   = {Peng, Yuxiang and Young, Jacob and Liu, Pengyu and Wu, Xiaodi},
  journal  = {arXiv preprint arXiv:2303.02775},
  year     = {2023}
}
```