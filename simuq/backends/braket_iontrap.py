# Currently this file is not usable:
# AWS Braket does not support arbitrary-angle MS gate for now
# Will be tested later
# If you want to use IonQ's device through braket,
# try qiskit_iontrap.py which generates a circuit in bk_c.

import networkx as nx
import numpy as np
from braket.circuits import Circuit

def clean_as(n, boxes, edges) :
    link = [(i, j) for i in range(n) for j in range(i + 1, n)]
    DG = nx.DiGraph()
    DG.add_nodes_from([i for i in range(len(boxes))])
    DG.add_edges_from(edges)
    topo_order = list(nx.topological_sort(DG))
    circ = Circuit()
    accum_phase = [0 for i in range(n)]
    
    for i in range(len(boxes)) :
        idx = topo_order[i]
        t = boxes[idx][1]
        for ((line, ins), params) in boxes[idx][0] :
            if line < n:
                if ins == 0:
                    q = line
                    rot = 2 * params[0] * t
                    phi = params[1]
                    if abs(rot) > 1e-5:
                        # Rz(q, phi)
                        accum_phase[q] -= phi
                        accum_phase[q] %= (2 * np.pi)
                        # Rx(q, rot)
                        if abs(rot / (2 * np.pi) - 0.25) < 1e-6:
                            circ.gpi2(q, (accum_phase[q] + 0) % (2 * np.pi))
                        elif abs(rot / (2 * np.pi) + 0.25) < 1e-6:
                            circ.gpi2(q, (accum_phase[q] + np.pi) % (2 * np.pi))
                        elif abs(rot / (2 * np.pi) - 0.5) < 1e-6:
                            circ.gpi(q, (accum_phase[q] + 0) % (2 * np.pi))
                        elif abs(rot / (2 * np.pi) + 0.5) < 1e-6:
                            circ.gpi(q, (accum_phase[q] + np.pi) % (2 * np.pi))
                        else:
                            circ.gpi2(q, (accum_phase[q] + 3 * np.pi / 2) % (2 * np.pi))
                            accum_phase[q] -= rot
                            accum_phase[q] %= (2 * np.pi)
                            circ.gpi2(q, (accum_phase[q] + np.pi / 2) % (2 * np.pi))
                        # Rz(q, -phi)
                        accum_phase[q] += phi
                        accum_phase[q] %= (2 * np.pi)
                else :
                    q = line
                    # Rz(q, 2 * params[0] * t)
                    accum_phase[q] -= 2 * params[0] * t
                    accum_phase[q] %= (2 * np.pi)
            else:
                (q0, q1) = link[line - n]
                theta = 2 * params[0] * t
                if ins == 0:
                    # R_XX(theta)
                    if abs(theta) > 1e-5:
                        if (theta / (2 * np.pi)) % 1 <= 0.25 or (theta / (2 * np.pi)) % 1 >= 0.75:
                            circ.ms(q0, q1, accum_phase[q0], accum_phase[q1], theta % (2 * np.pi))
                        elif 0.25 <= (theta / (2 * np.pi)) % 1 <= 0.5:
                            circ.ms(q0, q1, accum_phase[q0], accum_phase[q1], (theta % (2 * np.pi)) / 2)
                            circ.ms(q0, q1, accum_phase[q0], accum_phase[q1], (theta % (2 * np.pi)) / 2)
                        elif 0.5 <= (theta / (2 * np.pi)) % 1 <= 0.75:
                            circ.ms(q0, q1, accum_phase[q0], accum_phase[q1], ((theta % (2 * np.pi)) / 2 - np.pi) % (2 * np.pi))
                            circ.ms(q0, q1, accum_phase[q0], accum_phase[q1], ((theta % (2 * np.pi)) / 2 - np.pi) % (2 * np.pi))
                        else:
                            raise ValueError(f"Rotation angle is {theta}, should be between 0 and 2*pi")
                elif ins == 1:
                    # R_YY(theta)
                    if abs(theta) > 1e-5:
                        if (theta / (2 * np.pi)) % 1 <= 0.25 or (theta / (2 * np.pi)) % 1 >= 0.75:
                            circ.ms(q0, q1, (accum_phase[q0] + np.pi / 2) % (2 * np.pi), (accum_phase[q1] + np.pi / 2) % (2 * np.pi), theta % (2 * np.pi))
                        elif 0.25 <= (theta / (2 * np.pi)) % 1 <= 0.5:
                            circ.ms(q0, q1, (accum_phase[q0] + np.pi / 2) % (2 * np.pi), (accum_phase[q1] + np.pi / 2) % (2 * np.pi), (theta % (2 * np.pi)) / 2)
                            circ.ms(q0, q1, (accum_phase[q0] + np.pi / 2) % (2 * np.pi), (accum_phase[q1] + np.pi / 2) % (2 * np.pi), (theta % (2 * np.pi)) / 2)
                        elif 0.5 <= (theta / (2 * np.pi)) % 1 <= 0.75:
                            circ.ms(q0, q1, (accum_phase[q0] + np.pi / 2) % (2 * np.pi), (accum_phase[q1] + np.pi / 2) % (2 * np.pi), ((theta % (2 * np.pi)) / 2 - np.pi) % (2 * np.pi))
                            circ.ms(q0, q1, (accum_phase[q0] + np.pi / 2) % (2 * np.pi), (accum_phase[q1] + np.pi / 2) % (2 * np.pi), ((theta % (2 * np.pi)) / 2 - np.pi) % (2 * np.pi))
                        else:
                            raise ValueError(f"Rotation angle is {theta}, should be between 0 and 2*pi")
                else:
                    # R_ZZ(theta)
                    if abs(theta) > 1e-5:
                        # R_X(-pi/2)
                        circ.gpi2(q0, (accum_phase[q0] + np.pi) % (2 * np.pi))
                        circ.gpi2(q1, (accum_phase[q1] + np.pi) % (2 * np.pi))
                        # R_YY(theta)
                        if (theta / (2 * np.pi)) % 1 <= 0.25 or (theta / (2 * np.pi)) % 1 >= 0.75:
                            circ.ms(q0, q1, (accum_phase[q0] + np.pi / 2) % (2 * np.pi), (accum_phase[q1] + np.pi / 2) % (2 * np.pi), theta % (2 * np.pi))
                        elif 0.25 <= (theta / (2 * np.pi)) % 1 <= 0.5:
                            circ.ms(q0, q1, (accum_phase[q0] + np.pi / 2) % (2 * np.pi), (accum_phase[q1] + np.pi / 2) % (2 * np.pi), (theta % (2 * np.pi)) / 2)
                            circ.ms(q0, q1, (accum_phase[q0] + np.pi / 2) % (2 * np.pi), (accum_phase[q1] + np.pi / 2) % (2 * np.pi), (theta % (2 * np.pi)) / 2)
                        elif 0.5 <= (theta / (2 * np.pi)) % 1 <= 0.75:
                            circ.ms(q0, q1, (accum_phase[q0] + np.pi / 2) % (2 * np.pi), (accum_phase[q1] + np.pi / 2) % (2 * np.pi), ((theta % (2 * np.pi)) / 2 - np.pi) % (2 * np.pi))
                            circ.ms(q0, q1, (accum_phase[q0] + np.pi / 2) % (2 * np.pi), (accum_phase[q1] + np.pi / 2) % (2 * np.pi), ((theta % (2 * np.pi)) / 2 - np.pi) % (2 * np.pi))
                        else:
                            raise ValueError(f"Rotation angle is {theta}, should be between 0 and 2*pi")
                        # R_X(pi/2)
                        circ.gpi2(q0, (accum_phase[q0]) % (2 * np.pi))
                        circ.gpi2(q1, (accum_phase[q1]) % (2 * np.pi))

    circ = Circuit().add_verbatim_box(circ)
    return circ, accum_phase


def transpile(alignment, sol_gvars, boxes, edges) :
    n = len(alignment)
    circ, accum_phase = clean_as(n, boxes, edges)
    return circ
