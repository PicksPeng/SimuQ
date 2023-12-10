from abc import ABC, abstractmethod

import networkx as nx
import numpy as np
import scipy as sp

from simuq.backends.ionq_circuit import IonQCircuit
from simuq.transpiler import Transpiler
from qiskit import QuantumCircuit
from scipy.optimize import dual_annealing


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
    def clean_as(n, boxes, edges, circ: IonQCircuit, trotter_args, randomized=False) -> IonQCircuit:
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
                    if trotter_args["order"] == 1 and trotter_args["sequential"]:
                        trotter_mode = "first_order"
                    elif trotter_args["order"] == 2 and trotter_args["sequential"]:
                        trotter_mode = "second_order"
                    else:
                        trotter_mode = "random"
                    trotter_num = trotter_args["num"]

                    (q0, q1) = link[line - n]
                    params = 2 * np.array(params) * t
                    decomposed_params = decompose_ham(params)
                    n_terms = decomposed_params.shape[0]
                    if trotter_mode == "first_order" or n_terms == 1:
                        for _ in range(trotter_num):
                            for term_idx in range(n_terms):
                                u0, theta0 = rotate_to_x(decomposed_params[term_idx, 0, :])
                                circ._add_unitary(q0, u0)
                                u1, theta1 = rotate_to_x(decomposed_params[term_idx, 1, :])
                                circ._add_unitary(q1, u1)

                                circ.ms(
                                    q0,
                                    q1,
                                    0,
                                    0,
                                    theta0 * theta1 / trotter_num,
                                )
                                circ._add_unitary(q0, u0.conj().transpose())
                                circ._add_unitary(q1, u1.conj().transpose())
                    elif trotter_mode == "second_order":
                        for trotter_id in range(trotter_num):
                            for term_idx in range(n_terms):
                                if trotter_id % 2 == 0:
                                    term_idx = n_terms - 1 - term_idx
                                u0, theta0 = rotate_to_x(decomposed_params[0, :, term_idx])
                                circ._add_unitary(q0, u0)
                                u1, theta1 = rotate_to_x(decomposed_params[1, :, term_idx])
                                circ._add_unitary(q1, u1)

                                circ.ms(
                                    q0,
                                    q1,
                                    0,
                                    0,
                                    theta0 * theta1 / trotter_num,
                                )
                                circ._add_unitary(q0, u0.conj().transpose())
                                circ._add_unitary(q1, u1.conj().transpose())
                    elif trotter_mode == "random":
                        for _ in range(trotter_num):
                            for term_idx in np.random.permutation(n_terms):
                                u0, theta0 = rotate_to_x(decomposed_params[0, :, term_idx])
                                circ._add_unitary(q0, u0)
                                u1, theta1 = rotate_to_x(decomposed_params[1, :, term_idx])
                                circ._add_unitary(q1, u1)

                                circ.ms(
                                    q0,
                                    q1,
                                    0,
                                    0,
                                    theta0 * theta1 / trotter_num,
                                )
                                circ._add_unitary(q0, u0.conj().transpose())
                                circ._add_unitary(q1, u1.conj().transpose())
        return circ.optimize()

    def transpile(
        self,
        n,
        sol_gvars,
        boxes,
        edges,
        trotter_args,
        *generate_circuit_args,
        **generate_circuit_kwargs
    ):
        n = len(n) if isinstance(n, list) else n
        if "randomized" in generate_circuit_kwargs:
            randomized = generate_circuit_kwargs["randomized"]
            del generate_circuit_kwargs["randomized"]
        else:
            randomized = False
        circ = self.generate_circuit(n, *generate_circuit_args, **generate_circuit_kwargs)
        return self.clean_as(n, boxes, edges, circ, trotter_args, randomized)


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
        unitary = eigvecs @ X
    # print(unitary @ unitary.conj().T)
    # ham=unitary @ Z @ unitary^dagger
    unitary = Hadmard @ unitary.conj().T
    return unitary, theta


def decompose_ham(thetas):
    # change this to scipy
    X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    Z = np.array([[1, 0], [0, -1]], dtype=np.complex128)
    I = np.array([[1, 0], [0, 1]], dtype=np.complex128)
    H = (
        thetas[0] * np.kron(X, X)
        + thetas[1] * np.kron(X, Y)
        + thetas[2] * np.kron(X, Z)
        + thetas[3] * np.kron(Y, X)
        + thetas[4] * np.kron(Y, Y)
        + thetas[5] * np.kron(Y, Z)
        + thetas[6] * np.kron(Z, X)
        + thetas[7] * np.kron(Z, Y)
        + thetas[8] * np.kron(Z, Z)
    )
    total_amp = np.sqrt(np.trace(H @ H.conj().T))

    def create_tensor_ham(amp, theta1, theta2, phi1, phi2):
        H1 = amp * (
            np.cos(theta1) * X
            + np.sin(theta1) * np.cos(theta2) * Y
            + np.sin(theta1) * np.sin(theta2)
        )
        H2 = amp * (
            np.cos(phi1) * X + np.sin(phi1) * np.cos(phi2) * Y + np.sin(phi1) * np.sin(phi2)
        )
        return np.kron(H1, H2)

    def create_zero_ham():
        return np.zeros((4, 4), dtype=np.complex128)

    def calc_loss(params):
        params = params.reshape(-1, 5)
        n_terms = params.shape[0]
        guess = create_zero_ham()
        for i in range(n_terms):
            guess += create_tensor_ham(
                params[i, 0], params[i, 1], params[i, 2], params[i, 3], params[i, 4]
            )
        return np.linalg.norm(guess - H)

    def transform_params(params):
        params = params.reshape(-1, 5)
        n_terms = params.shape[0]
        new_params = np.zeros((n_terms, 2, 3))
        for i in range(n_terms):
            new_params[i, 0, 0] = params[i, 0] * np.cos(params[i, 1])
            new_params[i, 0, 1] = params[i, 0] * np.sin(params[i, 1]) * np.cos(params[i, 2])
            new_params[i, 0, 2] = params[i, 0] * np.sin(params[i, 1]) * np.sin(params[i, 2])
            new_params[i, 1, 0] = params[i, 0] * np.cos(params[i, 3])
            new_params[i, 1, 1] = params[i, 0] * np.sin(params[i, 3]) * np.cos(params[i, 4])
            new_params[i, 1, 2] = params[i, 0] * np.sin(params[i, 3]) * np.sin(params[i, 4])

        return new_params

    # def cal_commutator(A, B):
    #     return A @ B - B @ A
    for i in range(1, 4):
        bounds = np.array(
            [(0, total_amp), (0, np.pi), (0, np.pi), (0, np.pi), (0, np.pi)] * i, dtype=np.float64
        )
        x0 = (
            0.5
            * np.random.rand(5 * i)
            * np.array([total_amp, np.pi, np.pi, np.pi, np.pi] * i, dtype=np.float64)
        )
        res = dual_annealing(
            calc_loss,
            bounds=bounds,
            x0=x0,
        )
        if res.fun < 0.01:
            print("success with ", i, " terms")
            break
    new_params = transform_params(res.x)
    print(new_params)
    return transform_params(res.x)
