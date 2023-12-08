from abc import ABC, abstractmethod

import networkx as nx
import numpy as np
import scipy as sp

from simuq.backends.ionq_circuit import IonQCircuit
from simuq.transpiler import Transpiler
import torch
from qiskit import QuantumCircuit


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
                            U = sp.linalg.expm(
                                -1j
                                * (rot / 2)
                                * np.array([[0, cosphi - 1j * sinphi], [cosphi + 1j * sinphi, 0]])
                            )
                            circ._add_unitary(q, U)
                    else:
                        q = line
                        # expm(-i*(rot/2)*Z)
                        rot = 2 * params[0] * t
                        circ.rz(q, rot)
                else:
                    (q0, q1) = link[line - n]
                    params = 2 * np.array(params) * t
                    decomposed_params = decompose_ham(params)
                    n_terms = decomposed_params.shape[2]
                    # TODO add trotter
                    for term_idx in range(n_terms):
                        u0, theta0 = rotate_to_x(decomposed_params[0, :, term_idx])
                        circ._add_unitary(q0, u0)
                        u1, theta1 = rotate_to_x(decomposed_params[1, :, term_idx])
                        circ._add_unitary(q1, u1)

                        circ.ms(
                            q0,
                            q1,
                            0,
                            0,
                            theta0 * theta1,
                        )
                        circ._add_unitary(q0, u0.conj().transpose())
                        circ._add_unitary(q1, u1.conj().transpose())
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


def rotate_to_x(params):
    theta = np.sqrt((params**2).sum())
    params = params / theta
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.array([[1, 0], [0, -1]])
    Hadmard = 1 / np.sqrt(2) * np.array([[1, 1], [1, -1]])
    ham = params[0] * X + params[1] * Y + params[2] * Z
    eigvals, eigvecs = np.linalg.eigh(ham)
    if eigvals[0] > eigvals[1]:
        unitary = eigvecs
    else:
        unitary = eigvecs[:, ::-1]
    unitary = unitary @ Hadmard
    return unitary, theta


def decompose_ham(thetas):
    X = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)
    Y = torch.tensor([[0, -1j], [1j, 0]], dtype=torch.complex128)
    Z = torch.tensor([[1, 0], [0, -1]], dtype=torch.complex128)
    I = torch.tensor([[1, 0], [0, 1]], dtype=torch.complex128)
    H = (
        thetas[0] * torch.kron(X, X)
        + thetas[1] * torch.kron(X, Y)
        + thetas[2] * torch.kron(X, Z)
        + thetas[3] * torch.kron(Y, X)
        + thetas[4] * torch.kron(Y, Y)
        + thetas[5] * torch.kron(Y, Z)
        + thetas[6] * torch.kron(Z, X)
        + thetas[7] * torch.kron(Z, Y)
        + thetas[8] * torch.kron(Z, Z)
    )

    def create_local_ham(x, y, z):
        return x * X + y * Y + z * Z

    def create_zero_ham():
        return torch.zeros(4, 4, dtype=torch.complex128)

    def cal_commutator(A, B):
        return A @ B - B @ A

    for n in range(1, 4):
        params = torch.rand((2, 3, n), dtype=torch.float, requires_grad=True)
        optimizer = torch.optim.Adam([params], lr=0.01)
        for i in range(2000):
            optimizer.zero_grad()
            guess = create_zero_ham()
            commutator = 0
            Hamiltonian_terms = []
            for j in range(n):
                A = create_local_ham(params[0, 0, j], params[0, 1, j], params[0, 2, j])
                B = create_local_ham(params[1, 0, j], params[1, 1, j], params[1, 2, j])
                Hamiltonian_terms.append(torch.kron(A, B))
                C = torch.kron(A, B)

                guess += C
            for j in range(n):
                for k in range(j + 1, n):
                    commutator += cal_commutator(Hamiltonian_terms[j], Hamiltonian_terms[k]).norm()
            l = (guess - H).norm() + 0.001 * torch.sum(abs(params))
            # + 0.01 * commutator
            l.backward()
            optimizer.step()

        if (guess - H).norm() < 0.01:
            print("decomposable with", n, "terms")
            # print(params)
            # print(guess)
            break
    return params.detach().numpy()
