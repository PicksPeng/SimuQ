---
layout: page
menubar: docs_menu
title: Installation
subtitle: Prerequisites and steps to install SimuQ
show_sidebar: false
toc: true
---

## Download

Download SimuQ from [GitHub](https://github.com/PicksPeng/SimuQ) or [.zip](https://github.com/PicksPeng/SimuQ/archive/refs/heads/main.zip).

## Prerequisites

We have tested our implementation with `Python 3.9` and the following packages:
* `NetworkX 2.8.7`
* `NumPy 1.23.4`
* `SciPy 1.9.3`

To generate quantum circuits in Qiskit, IBM machine AAIS for specific machines, or programs in Qiskit Pulses, the following packages are also required:
* `Qiskit 0.38.0`
* `Qiskit-Terra 0.21.2`

The generated Bloqade programs are tested in `Julia 1.7.2` with `Bloqade 0.1.13`.

## Installation

Since Qiskit 0.38.0 requires Python version >= 3.9, if you want to use qiskit-related features, please update your Python to at least 3.9.

### Linux

To install the prerequisites, run in terminal:
```
pip install networkx numpy scipy
```

Furthermore, to run qiskit-related code, run in terminal:
```
pip install qiskit==0.38.0 qiskit-terra
```

We recommend Qiskit version 0.38.0 because there are known bugs in Qiskit 0.39.0.

For the installation of Julia, please refer to [Julialang website](https://julialang.org/downloads/platform/).

To install Bloqade, input `]` in Julia to enter the package manager, and run `add Bloqade`.

### Other platforms

To be added.