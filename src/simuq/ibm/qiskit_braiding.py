import networkx as nx
from qiskit import QuantumCircuit

from simuq.ibm.qiskit_transpiler import RYXGate, eYZXGate, get_pm

n = 3


def clean_as(n, boxes, edges):
    link = [(0, 1), (1, 2)]
    print(link)
    circ = QuantumCircuit(3)
    DG = nx.DiGraph()
    DG.add_nodes_from([i for i in range(len(boxes))])
    DG.add_edges_from(edges)
    topo_order = list(nx.topological_sort(DG))
    for i in range(len(boxes)):
        idx = topo_order[i]
        t = boxes[idx][1]
        for (line, ins), params in boxes[idx][0]:
            if line < n:
                if ins == 0:
                    q = line
                    # circ.rgate(2 * params[0] * t, params[1], q)
                    circ.rz(-params[1], q)
                    circ.rx(2 * params[0] * t, q)
                    circ.rz(params[1], q)
                else:
                    q = line
                    lamb = 2 * params[0] * t
                    circ.rz(lamb, q)
            elif line < n + len(link):
                (q0, q1) = link[line - n]
                theta = 2 * params[0] * t
                circ.append(RYXGate(theta), [q0, q1], [])
            else:
                theta = params[0] * t
                circ.append(eYZXGate(theta), [0, 1, 2], [])

    # circ.measure_all()
    return circ


def transpile(backend, alignment, sol_gvars, boxes, edges):
    circ = clean_as(n, boxes, edges)
    pm = get_pm(backend, for_braiding=True)
    return pm.run(circ)
