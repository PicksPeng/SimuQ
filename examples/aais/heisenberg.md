---
title: Heisenberg-like AAIS
subtitle: Programming the Heisenberg-like AAIS for IBM machines and IonQ machines
layout: page
show_sidebar: false
---

## IBM transmon systems

The IBM superconducting machines have native single qubit X-Y plane evolution and virtual (free) Z evolution. Besides, there are ZX evolution implemented via echoing cross-resonance pulses ([details]({{ '_posts/2023-1-18-design-aais-2' | relative_url }}). By adding single qubit rotations, it effectively realizes XX, YY, and ZZ evolution, which appear more in target Hamiltonian. Since the cross-resonance pulses are for two neibouring qubits in a 2-D topology, the topology in the AAIS is the same as the machine topology.

## IonQ trapped ion systems

Similar to IBM transmon systems, IonQ has trapped ion systems with native single qubit X-Y plane evolution and virtual (free) Z evolution. (Effectively) slightly different from the IBM machines, trapped ion machines have XX evolution implemented via precisely controlled micro movements. With single qubit rotations, we can also build YY and ZZ evolution from it. The difference to IBM machines is on the topology. The XX evolution can happen between an arbitrary pair of qubits.

## Heisenberg-like AAIS

Because of the similarities between the effects of these systems, we design a Heisenberg-like AAIS for both of them. Let $$E$$ represent the topology graph of the device. For IBM machine, $$E$$ is the machine topology. For IonQ machines, $$E$$ is a complete graph.

For each edge $$(i, j)\in E$$, we create a signal line. Three derived instructions implementing XX, YY, ZZ evolution respectively are added to this signal line. These instructions cannot be applied simultaneously with other instructions utilizing these sites. Since their implementations include a sequence of operations on the drive lines of the two sites $$i, j$$ and eventually realize effective Hamiltonian $$XX, YY, ZZ$, blending these impelementations with other pulses will result in wrong effective Hamiltonian.

Besides, we create a signal line for each qubit dealing with the single qubit operations. The X-Y plane evolution is native to both systems in the rotating frame, hence a native instruction implementing $$\alpha(\cos(\phi)X+\sin(\phi)Y)$$ where $$\alpha, \phi$$ are local variables. The Z evolution is realized by free Z gates, where phase shifts are recorded in classical pulse generators. We model it as a derived instruction implementing $$\alpha Z$$ where $$\alpha$$ is a local variable.

