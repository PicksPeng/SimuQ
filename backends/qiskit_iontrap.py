import numpy as np
import networkx as nx
from qiskit import QuantumCircuit


def clean_as(n, boxes, edges) :
    link = [(i, j) for i in range(n) for j in range(i + 1, n)]
    print(link)
    circ = QuantumCircuit(n)
    DG = nx.DiGraph()
    DG.add_nodes_from([i for i in range(len(boxes))])
    DG.add_edges_from(edges)
    topo_order = list(nx.topological_sort(DG))
    for i in range(len(boxes)) :
        idx = topo_order[i]
        t = boxes[idx][1]
        for ((line, ins), params) in boxes[idx][0] :
            if line < n :
                if ins == 0 :
                    q = line
                    theta = 2 * params[0] * t
                    lamb = - params[1]
                    phi = params[1] + np.pi
                    circ.u(theta, lamb, phi, q)
                else :
                    q = line
                    lamb = 2 * params[0] * t
                    circ.rz(lamb, q)
            else :
                (q0, q1) = link[line - n]
                theta = 2 * params[0] * t
                if ins == 0 :
                    circ.rxx(theta, q0, q1)
                elif ins == 1 :
                    circ.ryy(theta, q0, q1)
                else :
                    circ.rzz(theta, q0, q1)
    return circ


def transpile(alignment, sol_gvars, boxes, edges) :
    circ = clean_as(5, boxes, edges)
    return circ
