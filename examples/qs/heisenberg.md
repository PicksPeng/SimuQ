---
title: Heisenberg model
subtitle: Programming Heisenberg model in HML
layout: page
show_sidebar: false
---

The quantum Heisenberg model, developed by Werner Heisenberg, is a statistical mechanical model used in the study of critical points and phase transitions of magnetic systems, in which the spins of the magnetic systems are treated quantum mechanically. 

The Hamiltonian of a Heisenberg model with n qubits can be described as 
$$H=-\frac{1}{2} \sum_{j=1}^N\left(J_x \sigma_j^x \sigma_{j+1}^x+J_y \sigma_j^y \sigma_{j+1}^y+J_z \sigma_j^z \sigma_{j+1}^z+h \sigma_j^z\right)$$
where $$J_x, J_y$$ and $$J_z$$ are real-valued coupling constants. This Hamiltonian can be described by the following HML code:

```python
import numpy as np
from simuq.environment import qubit
from simuq.qsystem import QSystem


def Heisenberg_model(n_qubits, Jx, Jy, Jz, h, T):
    qs = QSystem()
    ql = [qubit(qs) for i in range(n_qubits)]
    link = [(i, i + 1) for i in range(n_qubits - 1)]
    ham = 0
    for (q0, q1) in link:
        ham -= Jx * ql[q0].X * ql[q1].X / 2
        ham -= Jy * ql[q0].Y * ql[q1].Y / 2
        ham -= Jz * ql[q0].Z * ql[q1].Z / 2
    for i in range(n):
        ham -= ql[i].X / 2
    qs.add_evolution(ham, T)
    return qs
```