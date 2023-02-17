---
title: Rydberg AAIS
subtitle: Programming the AAIS for Rydberg machines
layout: page
show_sidebar: false
---

The design principles of AAIS of Rydberg atom arrays are illustrated in a [tutorial]({{ 'tutorial/2023/01/15/design-aais-1/' | relative_url }}).

For QuEra Bloqade backend, the AAIS considers an 1-D or 2-D Rydberg atom array with local laser control. Each atom is targeted by a laser beam whose detuning, amplitude and phase can be arbitrarily set as functions of time.

For real QuEra machines through Bloqade or AWS, the AAIS contains a 2-D Rydberg atom array with global laser control. A laser with fully configurable detuning, amplitude and phase targets the atoms uniformly.

We provide classes to generate machine object of Rydberg machines in git-repo:aais/rydberg*.py, where you can specify the dimension and number of atoms.