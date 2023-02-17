---
title: IBM backend implementation
subtitle: How to implement IBM AAIS on IBM devices
layout: page
show_sidebar: false
---

One can easily translate instruction calls of AAIS into single qubit and two qubit rotation gates in Qiskit.

Since the compiled abstract schedule may contain a large amount of short evolution fragments, the Qiskit compiler makes use of two CNOT gates to realize one rotation gate because of calibration convenience. We make use of the [pulse-efficient circuit transpilation](https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.3.043088) to shorten the total duration of generated pulse 