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

    # Generate unitary GPI2(a)*GPI(b)*GPI2(c)
    def add_gpi2gpigpi2(q, a, b, c) :
        circ.gpi2(q, c + accum_phase[q])
        circ.gpi(q, b + accum_phase[q])
        circ.gpi2(q, a + accum_phase[q])

    def add_hadamard(q) :
        add_gpi2gpigpi2(q, np.pi, -np.pi / 4, 0)
    
    for i in range(len(boxes)) :
        idx = topo_order[i]
        t = boxes[idx][1]
        for ((line, ins), params) in boxes[idx][0] :
            if line < n :
                if ins == 0 :
                    q = line
                    rot = params[0] * t
                    phi = params[1]
                    # expm(-i*rot*(cos(phi)X+sin(phi)Y))
                    add_gpi2gpigpi2(q, np.pi / 2 - phi,
                                    rot - np.pi / 2 - phi, np.pi / 2 - phi)
                else :
                    q = line
                    # Rz(q, 2 * params[0] * t)
                    accum_phase[q] += 2 * params[0] * t
            else :
                (q0, q1) = link[line - n]
                theta = 2 * params[0] * t
                if ins == 0 :
                    # R_XX(theta)
                    if abs(theta) > 1e-5 :
                        circ.ms(q0, q1, accum_phase[q0], accum_phase[q1], theta)
                elif ins == 1 :
                    # R_YY(theta)
                    if abs(theta) > 1e-5 :
                        # Hadamard on q0
                        add_hadamard(q0)
                        # S on q0
                        accum_phase[q0] += np.pi / 2
                        # Hadamard on q1
                        add_hadamard(q1)
                        # S on q1
                        accum_phase[q1] += np.pi / 2
                        
                        circ.ms(q0, q1, accum_phase[q0], accum_phase[q1], theta)

                        # Sdg on q0
                        accum_phase[q0] -= np.pi / 2
                        # Hadamard on q0
                        add_hadamard(q0)
                        # Sdg on q1
                        accum_phase[q1] -= np.pi / 2
                        # Hadamard on q1
                        add_hadamard(q1)
                else :
                    # R_ZZ(theta)
                    if abs(theta) > 1e-5 :
                        # Hadamard on q0
                        add_hadamard(q0)
                        # Hadamard on q1
                        add_hadamard(q1)
                        
                        circ.ms(q0, q1, accum_phase[q0], accum_phase[q1], theta)

                        # Hadamard on q0
                        add_hadamard(q0)
                        # Hadamard on q1
                        add_hadamard(q1)
    circ = Circuit().add_verbatim_box(circ)
    return circ, accum_phase


def transpile(alignment, sol_gvars, boxes, edges) :
    n = len(alignment)
    circ, accum_phase = clean_as(n, boxes, edges)
    return circ
