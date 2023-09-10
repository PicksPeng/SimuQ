# Currently this file is not usable:
# AWS Braket does not support arbitrary-angle MS gate for now
# Will be tested later
# If you want to use IonQ's device through braket,
# try qiskit_iontrap.py which generates a circuit in bk_c.

from math import atan2, cos, floor, pi, sin, sqrt

import networkx as nx
import numpy as np
from qiskit import QuantumCircuit

# import gates
# import utils
from qiskit_ionq import GPI2Gate, GPIGate, IonQProvider, MSGate


def clean_as(n, boxes, edges):
    link = [(i, j) for i in range(n) for j in range(i + 1, n)]
    DG = nx.DiGraph()
    DG.add_nodes_from([i for i in range(len(boxes))])
    DG.add_edges_from(edges)
    topo_order = list(nx.topological_sort(DG))
    circ = QuantumCircuit(n, n)
    accum_phase = [0 for i in range(n)]

    def add_gpi2gpigpi2(q, a, b, c):
        circ.append(GPI2Gate((c + accum_phase[q]) / (2 * np.pi)), [q])
        circ.append(GPIGate((b + accum_phase[q]) / (2 * np.pi)), [q])
        circ.append(GPI2Gate((a + accum_phase[q]) / (2 * np.pi)), [q])

    def add_hadamard(q):
        add_gpi2gpigpi2(q, np.pi, np.pi / 4, 0)

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
                    # circ.rz(-phi, q)
                    # circ.rx(2 * rot, q)
                    # circ.rz(phi, q)
                    if abs(rot) > 1e-5:
                        add_gpi2gpigpi2(q, np.pi / 2 - phi, rot - np.pi / 2 - phi, np.pi / 2 - phi)
                else:
                    q = line
                    # Rz(q, 2 * params[0] * t)
                    accum_phase[q] += 2 * params[0] * t
            elif line < n + len(link):
                (q0, q1) = link[line - n]
                theta = 2 * params[0] * t
                if ins == 0:
                    # R_YX(theta)
                    if abs(theta) > 1e-5:
                        # Hadamard on q0
                        add_hadamard(q0)
                        # S on q0
                        accum_phase[q0] += np.pi / 2

                        circ.append(
                            MSGate(
                                accum_phase[q0] / (2 * np.pi),
                                accum_phase[q1] / (2 * np.pi),
                                theta / (2 * np.pi),
                            ),
                            [q0, q1],
                        )

                        # Sdg on q0
                        accum_phase[q0] -= np.pi / 2
                        # Hadamard on q0
                        add_hadamard(q0)
            else:
                t = params[0] * t
                if t < 0 or t >= 2 * pi:
                    t -= floor(t / (2 * pi)) * 2 * pi
                repeat = floor(t / (pi / 2)) + 1
                t = t / repeat
                t1 = atan2(-sqrt(sin(2 * t)), 1) / 2
                t2 = atan2(+sqrt(sin(2 * t)), cos(t) - sin(t)) / 2

                def eYY(T, q0, q1):
                    # Hadamard on q0
                    add_hadamard(q0)
                    # S on q0
                    accum_phase[q0] += np.pi / 2
                    # Hadamard on q1
                    add_hadamard(q1)
                    # S on q1
                    accum_phase[q1] += np.pi / 2

                    circ.append(
                        MSGate(
                            accum_phase[q0] / (2 * np.pi),
                            accum_phase[q1] / (2 * np.pi),
                            T / np.pi,
                        ),
                        [q0, q1],
                    )

                    # Sdg on q0
                    accum_phase[q0] -= np.pi / 2
                    # Hadamard on q0
                    add_hadamard(q0)
                    # Sdg on q1
                    accum_phase[q1] -= np.pi / 2
                    # Hadamard on q1
                    add_hadamard(q1)

                def eXX(T, q0, q1):
                    circ.append(
                        MSGate(
                            accum_phase[q0] / (2 * np.pi),
                            accum_phase[q1] / (2 * np.pi),
                            T / np.pi,
                        ),
                        [q0, q1],
                    )

                eYY(t1, 0, 1)
                eXX(t2, 1, 2)
                eYY(t2, 0, 1)
                eXX(t1, 1, 2)

    circ.measure(np.arange(n), np.arange(n))
    return circ, accum_phase


def transpile(alignment, sol_gvars, boxes, edges):
    n = len(alignment)
    circ, accum_phase = clean_as(n, boxes, edges)
    return circ


def debug_circ():
    n = 2
    circ = QuantumCircuit(n, n)
    accum_phase = [0 for i in range(n)]

    def add_gpi2gpigpi2(q, a, b, c):
        circ.append(GPI2Gate((c + accum_phase[q]) / (2 * np.pi)), [q])
        circ.append(GPIGate((b + accum_phase[q]) / (2 * np.pi)), [q])
        circ.append(GPI2Gate((a + accum_phase[q]) / (2 * np.pi)), [q])

    def add_hadamard(q):
        add_gpi2gpigpi2(q, np.pi, -np.pi / 4, 0)

    # Generate unitary GPI2(a)*GPI(b)*GPI2(c)

    add_hadamard(0)
    circ.measure(np.arange(n), np.arange(n))
    return circ
