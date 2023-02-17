---
title: IBM backend implementation
subtitle: How to implement IBM AAIS on IBM devices
layout: page
show_sidebar: false
---

The AAIS for Rydberg atom arrays have native instructions directly corresponding to the pulse channels of the devices. The valuation of the local variables then can be read as piecewise constant pulses. The transpiler for QuEra's machine creates these pulses using Amazon Braket API of Analog Hamiltonian Simulation. Users can then send these pulses to the QuEra devices via AWS.