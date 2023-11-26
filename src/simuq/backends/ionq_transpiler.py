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

def decompose_4x4_hermitian(H):
    import numpy as np
    import torch
    H=torch.tensor(H,dtype=torch.complex128)
    X=torch.tensor([[0,1],[1,0]],dtype=torch.complex128)
    Y=torch.tensor([[0,-1j],[1j,0]],dtype=torch.complex128)
    Z=torch.tensor([[1,0],[0,-1]],dtype=torch.complex128)
    I=torch.tensor([[1,0],[0,1]],dtype=torch.complex128)

    def create_2x2_hermitian(x,y,z,i):
        return x*X+y*Y+z*Z+i*I

    def empty_4x4_hermitian():
        return torch.zeros(4,4,dtype=torch.complex128)

    def cal_commutator(A,B):
        return A@B-B@A
    for n in range(1,4):
        params=torch.rand((2,4,n),dtype=torch.float,requires_grad=True)
        optimizer = torch.optim.Adam([params], lr=0.1)
        for i in range(1000):
            optimizer.zero_grad()
            guess=empty_4x4_hermitian()
            commutator=0
            Hamiltonian_terms=[]
            for j in range(n):
                A = create_2x2_hermitian(params[0,0,j],params[0,1,j],params[0,2,j],params[0,3,j])
                B = create_2x2_hermitian(params[1,0,j],params[1,1,j],params[1,2,j],params[1,3,j])
                Hamiltonian_terms.append(torch.kron(A, B))
                C = torch.kron(A, B)
            
                guess+=C
            for j in range(n):
                for k in range(j+1,n):
                    commutator+=cal_commutator(Hamiltonian_terms[j],Hamiltonian_terms[k]).norm()
            l=(guess-H).norm()+0.01*torch.sum(abs(params))+0.01*commutator
            l.backward()
            optimizer.step()
            # weight decay

        if (guess-H).norm()<0.01:
            print("decomposable with",n,"terms")
            # print(params)
            # print(guess)
            break
    return params