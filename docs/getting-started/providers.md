---
layout: page
menubar: docs_menu
title: Compile and run quantum simulation on devices
subtitle: Use SimuQ's providers for easy deployments
show_sidebar: false
toc: true
---

We demonstrate how to employ SimuQ to simulate the 3-site Ising model on IonQ device through IonQ cloud. A more detailed discussion is in [Simulate Ising model on multiple platforms](/SimuQ/examples/2023/08/05/ising/).

## Install dependencies

To activate IonQ Quantum Cloud in SimuQ, we need to install several optional dependencies. This can be done by executing

```bash
pip install simuq[ionq]
```

## Declare IonQ providers

Accesses to IonQ Quantum Cloud require an API key from [IonQ](https://ionq.com/). It directly gives you access to IonQ's noisy and ideal simulators for free. 

Assuming that you have a string `API_key` storing it, we can declare an IonQ provider in SimuQ by

```python
from simuq.ionq import IonQProvider
iqp = IonQProvider(API_key)
```

Alternatively, you can also store the key in a file and pass the path to the file in the delcaration of the provider:

```python
iqp = IonQProvider(from_file = './API_key')
```

## Compile, run, and retrieve results

To compile and run the 3-site Ising system `qs` in [programming an Ising model](/SimuQ/docs/qsystem) on IonQ devices, we need to specify a device and an AAIS. By default, the device is Aria-1 and the AAIS in compilation is Heisenberg-AAIS. One can also pass detailed parameters in the compilation like error tolerance bound and number of Trotterization here. For simplicity, we execute with the default compiler configuration by

```python
iqp.compile(qs)
```

We then run the experiment on simulator, which executes

```python
iqp.run(on_simulator = True)
```

This will submit a job to IonQ Quantum Cloud using the compiled results, and print the job id and status. The simulator will account for the noise model of Aria-1. To retrieve the results, we wait for a short period of time (normally several seconds) for the cloud server to queue and process the simulation, and then run

```python
iqp.results()
```

This returns a dictionary, whose results may depend on calibration data changing by date. In our experiments, it returns

```python
{'000': 0.609375,
 '011': 0.1042480469,
 '101': 0.2021484375,
 '110': 0.08422851562}
 ```

To run on real devices, we can execute

```python
iqp.run(shots = 4096)
```

Here `shots` represent the number of repetition of the experiment. After queuing and executing, you may retrieve the results by calling `iqp.results()`.