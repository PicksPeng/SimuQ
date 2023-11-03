from abc import ABC, abstractmethod

import networkx as nx
import numpy as np
import scipy as sp

from simuq.backends.ionq_circuit import IonQCircuit
from simuq.transpiler import Transpiler


def randomized_topo_sort(G):
    n = len(G)
    ret = []
    cand = []
    for i in range(n):
        if G.in_degree(i) == 0:
            cand.append(i)
    for _ in range(n):
        j = cand[np.random.randint(len(cand))]
        cand.remove(j)
        ret.append(j)
        nxts = list(G.neighbors(j))
        for k in nxts:
            if G.in_degree(k) == 0:
                continue
            G.remove_edge(j, k)
            if G.in_degree(k) == 0:
                cand.append(k)
    return ret


class IonQTranspiler(Transpiler, ABC):
    @abstractmethod
    def generate_circuit(self, n: int, *args, **kwargs) -> IonQCircuit:
        pass

    @staticmethod
    def clean_as(n, boxes, edges, circ: IonQCircuit, randomized=False) -> IonQCircuit:
        link = [(i, j) for i in range(n) for j in range(i + 1, n)]
        dg = nx.DiGraph()
        dg.add_nodes_from([i for i in range(len(boxes))])
        dg.add_edges_from(edges)
        if randomized:
            topo_order = randomized_topo_sort(dg)
        else:
            topo_order = list(nx.topological_sort(dg))

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
                            cosphi = np.cos(phi)
                            sinphi = np.sin(phi)
                            U = sp.linalg.expm(-1j * (rot / 2) * np.array([[0, cosphi - 1j * sinphi], [cosphi + 1j * sinphi, 0]]))
                            circ._add_unitary(q, U)

                            """
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
                            """
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
        return circ.optimize()

    def transpile(
        self, n, sol_gvars, boxes, edges, *generate_circuit_args, **generate_circuit_kwargs
    ):
        n = len(n) if isinstance(n, list) else n
        if "randomized" in generate_circuit_kwargs:
            randomized = generate_circuit_kwargs["randomized"]
            del generate_circuit_kwargs["randomized"]
        else:
            randomized = False
        circ = self.generate_circuit(n, *generate_circuit_args, **generate_circuit_kwargs)
        return self.clean_as(n, boxes, edges, circ, randomized)
