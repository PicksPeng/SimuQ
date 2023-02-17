---
title: AAIS
subtitle: Programming abstract analog instruction set in AAIS Specification Language
layout: page
show_sidebar: false
---

We program several interesting AAISs in AAIS Specification Language. 

## Abstract Analog Instruction Set

We illustrate the structure of AAIS here.

The top structure of AAIS contains a list of signal lines. They model the signal carriers of the physical devices.

For each signal line, there are instructions belonging to it. These instructions correspond to engineered pulses delivered through the signal carriers. Instructions belonging to the same signal line should not be invoked at the same time, since they are deployed via the same signal carrier exclusively.

Instructions may have different properties. A native instruction should describe the effect of its parameters. Formally, a native instruction can be applied simultaneously with other native instructions, and the total effect of simultaneous applications are the summation of their Hamiltonian. A derived instruction, on the other hand, may be characterized by the effective Hamiltonian of an engineered pulse on several sites. Our assumption here is that, a derived instruction should not borrow other sites to realize its effect. Derived instructions are similar to quantum gates, and cannot be simultaneously applied with other instructions that share effected sites with it. 

Each instruction may take several _local_ variables as input parameters, and have its effect Hamiltonian expressed in terms of these variables.

Besides, there are also _global_ variables belonging to the AAIS, corresponding to those configurable parameters before the experiment runs. One can also set system Hamiltonian expressed by these global variables to an AAIS, representing the always-on effects of the system.

## Design Principles of AAIS

AAISs are designed specifically for different systems. Nevertheless, there are shared principles behind their designs:
* The instructions should be expressive enough for most target Hamiltonians.
* The instructions should be closely related to the physics effect of the pulses.
* Rarely used instructions should be avoided.

## Supported Devices and Their AAIS

Currently we designed two AAIS:
* [Rydberg AAIS]({{ 'examples/aais/rydberg/' | relative_url }}): for QuEra's Rydberg atom arrays.
* [Heisenberg-like AAIS]({{ 'examples/aais/rydberg/' | relative_url }}): for IBM's superconducting qubit systems and IonQ's trapped ion systems.

SimuQ generates code directly runnable for these devices according to the API provided by hardware providers:
* QuEra's machines: through Amazon Braket or QuEra Bloqade.
* IBM's machines: through IBM Quantum.
* IonQ's machines: through Amazon Braket (currently no pulse-level control).

We are working on the supports of Rigetti's superconducting devices and Oxford Quantum Circuits through AWS Braket.

Current AAISs model the programmability of machines we have access to. Nevertheless, the abstraction of AAIS readily captures programmability of most quantum devices in the near future. With better devices in the future and more complicated target system to simulate, one can design more complex AAIS within the scope of expresiveness of AAIS specification language. 
