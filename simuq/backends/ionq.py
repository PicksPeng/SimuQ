import networkx as nx
import numpy as np
from ionq_circuit import Circuit

def clean_as(n, boxes, edges, backend="simulator", noise_model=None):
    link = [(i, j) for i in range(n) for j in range(i + 1, n)]
    DG = nx.DiGraph()
    DG.add_nodes_from([i for i in range(len(boxes))])
    DG.add_edges_from(edges)
    topo_order = list(nx.topological_sort(DG))
    circ = Circuit("test", n, backend, noise_model)
    
    for i in range(len(boxes)):
        idx = topo_order[i]
        t = boxes[idx][1]
        for (line, ins), params in boxes[idx][0]:
            if line < n:
                if ins == 0:
                    q = line
                    # expm(-i*(rot/2)*(cos(phi)X+sin(phi)Y))
                    rot = 2 * params[0] * t
                    phi = params[1]
                    if abs(rot) > 1e-5:
                        circ.rz(q, phi)
                        
                        # Rx(q, rot)
                        turns = rot / (2 * np.pi)
                        if abs(turns - 0.25) < 1e-6:
                            circ.gpi2(q, 0)
                        elif abs(turns - 0.75) < 1e-6:
                            circ.gpi2(q, np.pi)
                        elif abs(turns - 0.5) < 1e-6:
                            circ.gpi(q, 0)
                        else:
                            circ.gpi2(q, 3 * np.pi / 2)
                            circ.rz(q, rot)
                            circ.gpi2(q, np.pi / 2)
                        
                        circ.rz(q, -phi)

                else:
                    q = line
                    # expm(-i*(rot/2)*Z)
                    rot = 2 * params[0] * t
                    circ.rz(q, rot)
            else:
                (q0, q1) = link[line - n]
                theta = 2 * params[0] * t
                if ins == 0:
                    # R_XX(theta)
                    if abs(theta) > 1e-5:
                        circ.ms(q0, q1, 0, 0, theta)
                elif ins == 1:
                    # R_YY(theta)
                    if abs(theta) > 1e-5:
                        circ.ms(q0, q1, np.pi / 2, np.pi / 2, theta)
                else:
                    # R_ZZ(theta)
                    if abs(theta) > 1e-5:
                        # R_X(-pi/2)
                        circ.gpi2(q0, np.pi)
                        circ.gpi2(q1, np.pi)
                        # R_YY(theta)
                        circ.ms(q0, q1, np.pi / 2, np.pi / 2, theta)
                        # R_X(pi/2)
                        circ.gpi2(q0, 0)
                        circ.gpi2(q1, 0)
    return circ


def transpile(n, sol_gvars, boxes, edges, backend="simulator", noise_model=None):
    if isinstance(n, list):
        n = len(n)
    circ = clean_as(n, boxes, edges, backend, noise_model)
    return circ.job
