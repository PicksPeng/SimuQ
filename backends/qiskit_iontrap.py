import networkx as nx
import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import RGate

def clean_as(n, boxes, edges) :
    link = [(i, j) for i in range(n) for j in range(i + 1, n)]
    #print(link)
    circ = QuantumCircuit(n)
    DG = nx.DiGraph()
    DG.add_nodes_from([i for i in range(len(boxes))])
    DG.add_edges_from(edges)
    topo_order = list(nx.topological_sort(DG))
    bk_c = '''from braket.circuits import Circuit
circ = Circuit()'''
    for i in range(len(boxes)) :
        idx = topo_order[i]
        t = boxes[idx][1]
        for ((line, ins), params) in boxes[idx][0] :
            if line < n :
                if ins == 0 :
                    q = line
                    circ.append(RGate(2 * params[0] * t, params[1]), [q])
                    circ.barrier()
                    #circ.rz(-params[1],q)
                    #circ.rx(2 * params[0] * t,q)
                    #circ.rz(params[1],q)
                    bk_c += f'.rz({q}, {-params[1]})'
                    bk_c += f'.rx({q}, {2 * params[0] * t})'
                    bk_c += f'.rz({q}, {params[1]})'
                else :
                    q = line
                    lamb = 2 * params[0] * t
                    circ.rz(lamb, q)
                    bk_c += f'.rz({q}, {lamb})'
            else :
                #print(line)
                (q0, q1) = link[line - n]
                theta = 2 * params[0] * t
                if ins == 0 :
                    circ.rxx(theta, q0, q1)
                    circ.barrier()
                    bk_c += f'.xx({q0}, {q1}, {theta})'
                elif ins == 1 :
                    circ.ryy(theta, q0, q1)
                    circ.barrier()
                    bk_c += f'.yy({q0}, {q1}, {theta})'
                else :
                    circ.rzz(theta, q0, q1)
                    circ.barrier()
                    bk_c += f'.zz({q0}, {q1}, {theta})'
    with open('ionq.braket', 'w') as f :
        print(bk_c, file = f)
    circ.measure_all()
    return circ


def transpile(n, alignment, sol_gvars, boxes, edges) :
    circ = clean_as(n, boxes, edges)
    return circ
