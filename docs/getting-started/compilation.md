---
layout: page
menubar: docs_menu
title: Compile to Qiskit Pulse Schedules
subtitle: Compile your target evolution to IBM Qiskit pulse schedules
show_sidebar: false
toc: true
---

Assume that we have `qs` storing the target system to simulate. We want to compile it to an IBM machine (`ibmq_jakarta` in this case).

## Import the machine object

We first import the IBM machine object from the `FakeJakarta` machine backend.

```python
from qiskit.providers.fake_provider import FakeJakarta
from aais.ibm import get_mach
mach = get_mach(FakeJakarta())
```

## Synthesize an abstract schedule

We then apply the SimuQ solver to solve and store the intermediate representation.

```python
from simuq.solver import generate_as
as = generate_as(qs, mach)
```

Here `as` is a tuple storing the abstract schedule, containing `(Mapping, ValuationGVars, InstructionBoxes, TemporalRelation)`.

## Generate pulse schedule

We then apply a transpiler to generate the pulse schedule.

```python
from backends.qiskit_pulse_ibm import transpile
schedule = transpile(backend, *as)
```

To submit this schedule to real IBM machines, please refer to [Qiskit Pulse](https://qiskit.org/documentation/apidoc/pulse.html).
