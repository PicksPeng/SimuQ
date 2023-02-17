---
title: Jaynes-Cumming model
subtitle: Programming Jaynes-Cumming model in HML
layout: page
show_sidebar: false
---

The Jaynes–Cummings model (sometimes abbreviated JCM) is a theoretical model in quantum optics. It describes the system of a two-level atom interacting with a quantized mode of an optical cavity (or a bosonic field), with or without the presence of light (in the form of a bath of electromagnetic radiation that can cause spontaneous emission and absorption). 

The Hamiltonian of a Jaynes–Cummings model with one qubits can be described as 
$$\hat{H}_{\mathrm{JC}}=\hbar \omega_c \hat{a}^{\dagger} \hat{a}+\hbar \omega_a \frac{\hat{\sigma}_z}{2}+\frac{\hbar \Omega}{2}\left(\hat{a} \hat{\sigma}_{+}+\hat{a}^{\dagger} \hat{\sigma}_{-}\right)$$
where $$\omega_c, \omega_a$$ and $$\Omega$$ are real-valued coupling constants. This Hamiltonian can be described by the following HML code:

```python
import numpy as np
from simuq.environment import qubit, boson
from simuq.qsystem import QSystem


def JC_model(omega_c, omega_a, Omega, T):
    qs = QSystem()
		q, bos = qubit(qs), boson(qs)

		ham = omega_c * bos.c * bos.a + omega_a / 2 * q.Z + Omega / 2 * (bos.a * (q.X - 1j * q.Y) + bos.c * (q.X + 1j * q.Y))

    qs.add_evolution(ham, T)
    return qs
```