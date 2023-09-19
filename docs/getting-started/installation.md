---
layout: page
menubar: docs_menu
title: Installation
subtitle: Prerequisites and steps to install SimuQ
show_sidebar: false
toc: true
---

## Download

You may install SimuQ directly via `pip`:

```python
pip install simuq
```

## Providers

Currently, most functionalities of SimuQ can be accessed via providers, interfaces for different cloud services. To enable a specific provider, for example Amazon Braket, you may install via

```python
pip install "simuq[braket]"
```

SimuQ now supports Amazon Braket `[braket]`, IonQ Quantum Cloud `[ionq]`, and IBM Quantum Experience `[ibm]`. For numerical simulation, it also supports QuTiP `[qutip]`.


## Installation from source

To install from source code, download SimuQ from [GitHub](https://github.com/PicksPeng/SimuQ) or [.zip](https://github.com/PicksPeng/SimuQ/archive/refs/heads/main.zip).

Enter the SimuQ folder which contains a `setup.py`, and execute

```python
pip install .
```

For specific backends and providers, you may optionally install other dependencies according to your need.
