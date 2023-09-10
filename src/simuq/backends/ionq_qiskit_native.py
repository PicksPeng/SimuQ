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
                    rot = 2 * params[0] * t
                    phi = params[1]
                    if abs(rot) > 1e-5:
                        # Rz(q, phi)
                        accum_phase[q] -= phi
                        accum_phase[q] %= 2 * np.pi
                        # Rx(q, rot)
                        if abs(rot / (2 * np.pi) - 0.25) < 1e-6:
                            circ.append(GPI2Gate(((accum_phase[q] + 0) / (2 * np.pi)) % 1, [q]))
                        elif abs(rot / (2 * np.pi) + 0.25) < 1e-6:
                            circ.append(GPI2Gate(((accum_phase[q] + np.pi) / (2 * np.pi)) % 1, [q]))
                        elif abs(rot / (2 * np.pi) - 0.5) < 1e-6:
                            circ.append(GPIGate(((accum_phase[q] + 0) / (2 * np.pi)) % 1, [q]))
                        elif abs(rot / (2 * np.pi) + 0.5) < 1e-6:
                            circ.append(GPIGate(((accum_phase[q] + np.pi) / (2 * np.pi)) % 1, [q]))
                        else:
                            circ.append(
                                GPI2Gate(((accum_phase[q] + 3 * np.pi / 2) / (2 * np.pi)) % 1, [q])
                            )
                            accum_phase[q] -= rot
                            accum_phase[q] %= 2 * np.pi
                            circ.append(
                                GPI2Gate(((accum_phase[q] + np.pi / 2) / (2 * np.pi)) % 1, [q])
                            )
                        # Rz(q, -phi)
                        accum_phase[q] += phi
                        accum_phase[q] %= 2 * np.pi
                else:
                    q = line
                    # Rz(q, 2 * params[0] * t)
                    accum_phase[q] -= 2 * params[0] * t
                    accum_phase[q] %= 2 * np.pi
            else:
                (q0, q1) = link[line - n]
                theta = 2 * params[0] * t
                if ins == 0:
                    # R_XX(theta)
                    if abs(theta) > 1e-5:
                        if (theta / (2 * np.pi)) % 1 <= 0.25 or (theta / (2 * np.pi)) % 1 >= 0.75:
                            circ.append(
                                MSGate(
                                    accum_phase[q0] / (2 * np.pi),
                                    accum_phase[q1] / (2 * np.pi),
                                    (theta / (2 * np.pi)) % 1,
                                    [q0, q1],
                                )
                            )
                        elif 0.25 <= (theta / (2 * np.pi)) % 1 <= 0.5:
                            circ.append(GPIGate(((accum_phase[q0]) / (2 * np.pi)) % 1, [q0]))
                            circ.append(GPIGate(((accum_phase[q1]) / (2 * np.pi)) % 1, [q1]))
                            circ.append(
                                MSGate(
                                    (accum_phase[q0] / (2 * np.pi) + 0.5) % 1,
                                    accum_phase[q1] / (2 * np.pi),
                                    0.5 - ((theta / (2 * np.pi)) % 1),
                                    [q0, q1],
                                )
                            )
                        elif 0.5 <= (theta / (2 * np.pi)) % 1 <= 0.75:
                            circ.append(GPIGate(((accum_phase[q0]) / (2 * np.pi)) % 1, [q0]))
                            circ.append(GPIGate(((accum_phase[q1]) / (2 * np.pi)) % 1, [q1]))
                            circ.append(
                                MSGate(
                                    accum_phase[q0] / (2 * np.pi),
                                    accum_phase[q1] / (2 * np.pi),
                                    ((theta / (2 * np.pi)) % 1) - 0.5,
                                    [q0, q1],
                                )
                            )
                        else:
                            raise ValueError(
                                f"Rotation angle is {theta}, should be between 0 and 2*pi"
                            )
                elif ins == 1:
                    # R_YY(theta)
                    if abs(theta) > 1e-5:
                        if (theta / (2 * np.pi)) % 1 <= 0.25 or (theta / (2 * np.pi)) % 1 >= 0.75:
                            circ.append(
                                MSGate(
                                    (accum_phase[q0] / (2 * np.pi) + 0.25) % 1,
                                    (accum_phase[q1] / (2 * np.pi) + 0.25) % 1,
                                    (theta / (2 * np.pi)) % 1,
                                    [q0, q1],
                                )
                            )
                        elif 0.25 <= (theta / (2 * np.pi)) % 1 <= 0.5:
                            circ.append(GPIGate(((accum_phase[q0]) / (2 * np.pi)) % 1, [q0]))
                            circ.append(GPIGate(((accum_phase[q1]) / (2 * np.pi)) % 1, [q1]))
                            circ.append(
                                MSGate(
                                    (accum_phase[q0] / (2 * np.pi) + 0.75) % 1,
                                    (accum_phase[q1] / (2 * np.pi) + 0.25) % 1,
                                    0.5 - ((theta / (2 * np.pi)) % 1),
                                    [q0, q1],
                                )
                            )
                        elif 0.5 <= (theta / (2 * np.pi)) % 1 <= 0.75:
                            circ.append(GPIGate(((accum_phase[q0]) / (2 * np.pi)) % 1, [q0]))
                            circ.append(GPIGate(((accum_phase[q1]) / (2 * np.pi)) % 1, [q1]))
                            circ.append(
                                MSGate(
                                    (accum_phase[q0] / (2 * np.pi) + 0.25) % 1,
                                    (accum_phase[q1] / (2 * np.pi) + 0.25) % 1,
                                    ((theta / (2 * np.pi)) % 1) - 0.5,
                                    [q0, q1],
                                )
                            )
                        else:
                            raise ValueError(
                                f"Rotation angle is {theta}, should be between 0 and 2*pi"
                            )
                else:
                    # R_ZZ(theta)
                    if abs(theta) > 1e-5:
                        # R_X(-pi/2)
                        circ.append(GPI2Gate(((accum_phase[q0] + np.pi) / (2 * np.pi)) % 1, [q0]))
                        circ.append(GPI2Gate(((accum_phase[q1] + np.pi) / (2 * np.pi)) % 1, [q1]))
                        # R_YY(theta)
                        if (theta / (2 * np.pi)) % 1 <= 0.25 or (theta / (2 * np.pi)) % 1 >= 0.75:
                            circ.append(
                                MSGate(
                                    (accum_phase[q0] / (2 * np.pi) + 0.25) % 1,
                                    (accum_phase[q1] / (2 * np.pi) + 0.25) % 1,
                                    (theta / (2 * np.pi)) % 1,
                                    [q0, q1],
                                )
                            )
                        elif 0.25 <= (theta / (2 * np.pi)) % 1 <= 0.5:
                            circ.append(GPIGate(((accum_phase[q0]) / (2 * np.pi)) % 1, [q0]))
                            circ.append(GPIGate(((accum_phase[q1]) / (2 * np.pi)) % 1, [q1]))
                            circ.append(
                                MSGate(
                                    (accum_phase[q0] / (2 * np.pi) + 0.75) % 1,
                                    (accum_phase[q1] / (2 * np.pi) + 0.25) % 1,
                                    0.5 - ((theta / (2 * np.pi)) % 1),
                                    [q0, q1],
                                )
                            )
                        elif 0.5 <= (theta / (2 * np.pi)) % 1 <= 0.75:
                            circ.append(GPIGate(((accum_phase[q0]) / (2 * np.pi)) % 1, [q0]))
                            circ.append(GPIGate(((accum_phase[q1]) / (2 * np.pi)) % 1, [q1]))
                            circ.append(
                                MSGate(
                                    (accum_phase[q0] / (2 * np.pi) + 0.25) % 1,
                                    (accum_phase[q1] / (2 * np.pi) + 0.25) % 1,
                                    ((theta / (2 * np.pi)) % 1) - 0.5,
                                    [q0, q1],
                                )
                            )
                        else:
                            raise ValueError(
                                f"Rotation angle is {theta}, should be between 0 and 2*pi"
                            )
                        # R_X(pi/2)
                        circ.append(GPI2Gate(((accum_phase[q0]) / (2 * np.pi)) % 1, [q0]))
                        circ.append(GPI2Gate(((accum_phase[q1]) / (2 * np.pi)) % 1, [q1]))
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
