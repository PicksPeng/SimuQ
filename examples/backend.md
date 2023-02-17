---
title: Backend
subtitle: Build backend supports for AAIS to generate pulse-level schedules
layout: page
show_sidebar: false
---

We illustrate implementations of several backends, translating schedules used in SimuQ solver to pulse schedules with knowledges of how machines implement AAIS.

## Abstract Schedules

SimuQ solver generates abstract schedules that reproduce the target Hamiltonian. An abstract schedule contains a valuation of the global variables and a temporal graph of instruction blocks. An instruction block is consisted by several instructions, their local variable valuations and an evolution time $$t$$. It represents these instructions are applied simultaneously for time $$t$$. The temporal graph shows the temporal relation between two instruction blocks: an instruction block should start executing after the finish of the other instruction block execution.

## Concrete Schedules and Pulse Schedules

The first component of a backend is a scheduler that settles down instruction executions' start and end times according to the durations of their implementations. It specifies the sequence of instructions executed on a signal line (including idle intervals). This is done by importing the implementation details from the machine and running a topological sort of the temporal graph. The resulting schedule is called concrete schedule.

Then a pulse scheduler transpiles the concrete schedule into a pulse schedule by transpiling each instruction call with the pulse implementation of the instruction.

## Supported Backends

* [QuEra backend]({{ 'examples/backend/quera/' | relative_url }}): using Rydberg AAIS with global control, through Amazon Braket or QuEra Bloqade.
* [IBM backend]({{ 'examples/backend/ibm/' | relative_url }}): using Heisenberg-like AAIS with IBM machines' topology, through IBM Quantum.