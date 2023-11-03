# SimuQ

This file contains necessary information for the artifact evaluation of paper "SimuQ: a Framework for Programming Quantum Hamiltonian Simulation with Analog Compilation".

**Disclaimer:** Our work help program and deploy quantum Hamiltonian simulation on cloud-based quantum computing devices. Most of the demonstrations related to artifact evaluation in this repository are deployed on classical simulators (either locally or through cloud-based services) for the convenience of the reviewers, since they are free. In our paper, we demonstrate many executions on real cloud-based quantum devices. This can be easily achieved by setting `on_simulator=False` when running the jobs. Be careful when you decide to run them on real quantum devices: (1) it inevitably requires your access credentials to the cloud service providers; (2) it costs real money; (3) and it is expensive (normally take dozens of dollars for an experiment).

For (mathematical) background knowledge of quantum Hamiltonian simulation, please refer to Section 2.1 of our paper. Detailed backgrounds are in Section 1.

## List of Claims

We propse SimuQ, a framework for programming quantum Hamiltonian simulation with analog compilation. Its design and implementation enable us to conduct the following case studies in Section 5.
- We can program and deploy the chain or cycle Ising model simulation on QuEra, IonQ, and IBM devices, as in Section 5.1.
- We can exploit native instructions on IBM devices to generate shorter pulse schedules for a specific quantum system, as in Section 5.2.
- We can exploit interaction-based gates on IonQ and IBM devices to reduce errors in an implementation of quantum approximate optimization algorithm, as in Section 5.3.
- We establish a quantum simulation benchmark and test SimuQ compiler on it, as in Section 5.4.

## Environment Setup

SimuQ requires Python 3 to run. We recommend Python 3.9 or greater since many optional dependencies have a minimum Python version requirement. We use `pip` for installation in this repository. In the following, we assume that you have a proper Python environment with `pip`. 

For easy reproduction of our case studies in the paper, we prepare notebooks for each of them. We recommend a local setup of Jupyter Notebook by

```bash
pip install notebook
```

Alternatively, you may choose online Jupyter Notebook services, like Google Colab and Amazon Web Services. For environment setup of these services, please refer to their documentation.

### Installation

To install SimuQ, please enter `SimuQ/` folder, where there is a `src/` folder containing the source code. The core of SimuQ may be installed by the following command, which normally takes a minute.

```bash
pip install .
```

SimuQ supports multiple quantum device platforms and cloud service providers. To run the case study notebooks, you may install all the related providers by the following command, which takes a few minutes to install.

```bash
pip install ".[all]"
```

If you only want to install a specific provider, please refer to `SimuQ/README.md` for detailed instructions.

### Sanity Check

For a simply sanity check, you may run the following Python code to check if you have successfully installed SimuQ for the artifact evaluation. The following code programs a simply quantum system and deploy it on a classical simulator based on package [QuTiP](https://qutip.org/).

```python
from simuq import QSystem, Qubit
from simuq.qutip import QuTiPProvider

qs = QSystem()
q0, q1 = Qubit(qs), Qubit(qs)
h = q0.X * q1.X
qs.add_evolution(h, 1.0)

qtp = QuTiPProvider()
qtp.compile(qs)
qtp.run()
print(qtp.results())
```

It takes a few seconds to run the code. The output should be

```bash
Compiled.
Solved.
{'00': 0.2919499925618451, '01': 0.0, '10': 0.0, '11': 0.7080500074381548}
```

It represents the ideal measurement results after evolving the quantum Hamiltonian system `qs`.

### Cloud Service Providers

Currently SimuQ supports three cloud-based quantum service providers: Amazon Braket, IonQ Quantum Cloud, and IBM Quantum Experience. Running case study notebooks does not require you to pay, since most of the tasks are deployed on simulators, which are totally free. However, you still need accounts from these providers to run the code (except for Amazon Braket provider, whose simulator runs locally).

The user interfaces are organized as `providers` with similar APIs. To run tasks on these providers and their devices, you might need their accounts and credits. We list the details below.

For Amazon Braket, SimuQ supports its access to IonQ's trapped ion devices and QuEra's neutral atom devices. Running on Braket classical simulators (for these platforms) does not require an account and can be executed locally. If you want to submit jobs to Amazon Braket accessing real devices, make sure that your AWS account is onboarded to Amazon Braket, as per the instructions [here](https://github.com/amazon-braket/amazon-braket-sdk-python#prerequisites).

For IonQ Quantum Cloud, SimuQ supports its access to IonQ's trapped ion devices. Running on IonQ's ideal and noisy classical simulators requires an [IonQ Quantum Cloud](https://ionq.com/quantum-cloud) account and is free. When creating IonQ providers in SimuQ, you must provide the API key corresponding your account from IonQ either through a string or a file storing the key. To run on real devices, please contact IonQ for access permission.

For IBM Quantum Experience, SimuQ supports its access to IBM's superconducting transmon devices. Running on IBM's ideal and noisy classical simulators requires an [IBM Quantum Experience](https://quantum-computing.ibm.com/) account and is free. Similarly, when creating IonQ providers in SimuQ, you must provide the API token corresponding your account from IBM either through a string or a file storing the key. IBM's open plan grants you access to small real devices for free. For the case studies, many of them require larger devices. Please refer to [IBM Quantum Experience](https://quantum-computing.ibm.com/) for more information on pricing if you want to run on large real devices.

## Evaluation Instructions

We have prepared a notebook for each case study. You may directly execute the notebooks line by line to evaluate our claims. The notebooks are in folder `notebooks/artifact_evaluation`. 

Note that for the case study on a quantum simulation benchmark, compiling large systems may take around an hour, depending on your computing system. It is expected that the compilation time is similar to our report in the paper. We expect that the examples in Table 3 can be compiled within an hour, except for possible time outs on your device for `ising_chain` of size 96 on QuEra and IBM, `ising_cycle` of size 64 on QuEra and IBM.

## Additional Artifact Description

You are welcome to read SimuQ tutorial notebooks in folder `notebooks/tutorials`, where we have prepared 4 notebooks from basical to advanced usage of SimuQ with detailed illustrations. These notebooks may also help you evaluate the reusability of our artifact.

We encourage you to interact with SimuQ providers for a smooth front-end user experience. You are also welcome to examine and explore the detailed implementation behind the providers, according to the following description of the artifact structure.

### Artifact Structure

`notebooks/`: Example notebooks of running SimuQ and obtaining results.

`src/simuq/`: The core compiler and language implementation of SimuQ.

`src/simuq/aais/`: AAIS of many backends in AAIS Specification Language.

`src/simuq/backends/`: The translation to API languages of different machine backends.

`src/simuq/systems/`: Hamiltonian evolution implemented in Hamiltonian Modeling Language.

`src/simuq/systems/benchmark/`: A small benchmark of quantum Hamiltonian simulation as reported in our paper.

`src/simuq/qutip/`: The QuTiP provider of SimuQ.

`src/simuq/braket/`: The Amazon Braket provider of SimuQ.

`src/simuq/ionq/`: The IonQ Quantum Cloud provider of SimuQ.

`src/simuq/ibm/`: The IBM Quantum Experience provider of SimuQ.

