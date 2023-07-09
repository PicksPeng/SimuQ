# Currently this file is not usable:
# AWS Braket does not support arbitrary-angle MS gate for now
# Will be tested later
# If you want to use IonQ's device through braket,
# try qiskit_iontrap.py which generates a circuit in bk_c.

import networkx as nx
import numpy as np
from qiskit import QuantumCircuit

# import gates
# import utils
from qiskit_ionq import GPI2Gate, GPIGate, IonQProvider, MSGate


def clean_as(n, boxes, edges):
    link = [(i, j) for i in range(n) for j in range(i + 1, n)]
    # link = [(i, (i + 1) % 12) for i in range(n)]
    DG = nx.DiGraph()
    DG.add_nodes_from([i for i in range(len(boxes))])
    DG.add_edges_from(edges)
    topo_order = list(nx.topological_sort(DG))
    circ = QuantumCircuit(n, n)
    accum_phase = [0 for i in range(n)]

    # Generate unitary GPI2(a)*GPI(b)*GPI2(c)

    for i in range(len(boxes)):
        idx = topo_order[i]
        t = boxes[idx][1]
        for (line, ins), params in boxes[idx][0]:
            if line < n:
                if ins == 0:
                    q = line
                    rot = params[0] * t
                    phi = params[1]
                    # expm(-i*rot*(cos(phi)X+sin(phi)Y))=I-i*rot*(cos(phi)X+sin(phi)Y)
                    circ.rz(-phi, q)
                    circ.rx(2 * rot, q)
                    circ.rz(phi, q)

                else:
                    q = line
                    # Rz(q, 2 * params[0] * t)
                    rot = params[0] * t
                    circ.rz(2 * rot, q)

            else:
                (q0, q1) = link[line - n]
                theta = 2 * params[0] * t
                if ins == 0:
                    # R_XX(theta)
                    if abs(theta) > 1e-5:
                        circ.rxx(theta, q0, q1)
                elif ins == 1:
                    # R_YY(theta)
                    if abs(theta) > 1e-5:
                        # Hadamard on q0
                        circ.ryy(theta, q0, q1)
                else:
                    # R_ZZ(theta)
                    if abs(theta) > 1e-5:
                        # Hadamard on q0
                        circ.rzz(theta, q0, q1)
    for i in range(12):
        circ.h(i)
    circ.measure(np.arange(n), np.arange(n))
    return circ, accum_phase


def transpile(alignment, sol_gvars, boxes, edges):
    n = len(alignment)
    circ, accum_phase = clean_as(n, boxes, edges)
    return circ
