---
layout: page
menubar: docs_menu
title: Simulate your quantum system
subtitle: Compile and run via SimuQ's QuTiP provider
show_sidebar: false
toc: true
---


Following [programming an Ising model](/SimuQ/docs/qsystem), we have `qs` storing the target system to simulate. We numerically simulate it via SimuQ's QuTiP provider, which employs a translation from SimuQ to QuTiP and solve the Schrodinger equation.

## Install dependencies

To enable QuTiP providers, we need to install [QuTiP](https://qutip.org/) first. This can be easily done by running

```bash
pip install simuq[qutip]
```

## Create a QuTiP provider

A provider is a user interface for convenient manipulations of functionalities of SimuQ. We use QuTiP provider as a basic example on how to use providers to deploy quantum simulation problems on devices and obtain results.

We can create a QuTiP provider via the following code

```python
from simuq.providers import QuTiPProvider
qpp = QuTiPProvider()
```

## Compilation in provider

To simulate a quantum system `qs` programmed in HML via SimuQ, we need three major steps of a provider: `compile`, `run`, `results`.

We call the `compile` function of the provider to process the system into a runnable executable. For QuTiP provider, we can execute

```python
qpp.compile(qs)
```

QuTiP provider processes the quantum system `qs` and translate it into a Hamiltonian in QuTiP. When succeeds, a message `Compiled.` will show up.

For other providers, `compile` command may specify the backend device, AAIS, and compiler specifications.

When compilation succeeds, the job will be recorded in the provider.

## Run and obtain results from providers

Running a job will send the compilation results to backend devices to execute. For QuTiP provider, we execute

```python
qpp.run()
```

When successfully executed, the command will show a message `Submitted`, meaning that the executables are sent to devices to run.

For other providers, you can choose whether to run on classical simulators, and how many shots to run on the real device or simulators. Note that these classical simulators are mostly from hardware providers and executed on the cloud.

To retrieve the results, we can execute

```python
qpp.results()
```

If the task is still in the queue on the device, an exception will be raised. When the retrieval succeeds, a dictionary is returned, which contains the frequencies of obtaining a measurement array (encoded as a 0/1 string). A bit in the string corresponds to a site of the quantum system. We can call `qpp.print_sites()` to show the order of the sites in the measurement output.

The resulting dictionary of the 3-site Ising model is

```python
{'000': 0.6972053600815633,
 '001': 0.0,
 '010': 0.0,
 '011': 0.05973358916427098,
 '100': 0.0,
 '101': 0.18332746158989466,
 '110': 0.05973358916427098,
 '111': 0.0}
```
