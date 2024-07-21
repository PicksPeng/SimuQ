from abc import ABC, abstractmethod

import networkx as nx
import numpy as np
import scipy as sp
from scipy.optimize import dual_annealing

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


class IonQTranspiler_2Pauli(Transpiler, ABC):
    @abstractmethod
    def generate_circuit(self, n: int, *args, **kwargs) -> IonQCircuit:
        pass

    @staticmethod
    def clean_as(n, boxes, edges, circ: IonQCircuit, trotter_args) -> IonQCircuit:
        link = [(i, j) for i in range(n) for j in range(i + 1, n)]
        dg = nx.DiGraph()
        dg.add_nodes_from([i for i in range(len(boxes))])
        dg.add_edges_from(edges)
        seen_params = {}
        precision = 1e-5
        if trotter_args["randomized"]:
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
                    # TODO add optimization, if already decomposed, use a dict to store the decomposed params
                    key = (params / precision).astype(int).tostring()
                    if key in seen_params:
                        decomposed_params = seen_params[key]
                    else:
                        decomposed_params = decompose_ham(params)
                        seen_params[key] = decomposed_params
                    n_terms = decomposed_params.shape[0]

                    def apply_linear_tensor(decomposed_params_in_term):
                        u0, theta0 = rotate_to_x(decomposed_params_in_term[0, :])
                        u1, theta1 = rotate_to_x(decomposed_params_in_term[1, :])
                        circ._add_unitary(q0, u0)
                        circ._add_unitary(q1, u1)

                        circ.ms(q0, q1, 0, 0, theta0 * theta1 / trotter_num)

                        circ._add_unitary(q0, u0.conj().transpose())
                        circ._add_unitary(q1, u1.conj().transpose())

                    if trotter_mode == "first_order" or n_terms == 1:
                        for _ in range(trotter_num):
                            for term_idx in range(n_terms):
                                apply_linear_tensor(decomposed_params[term_idx])

                    elif trotter_mode == "second_order":
                        for trotter_id in range(trotter_num):
                            for term_idx in range(n_terms):
                                if trotter_id % 2 == 0:
                                    term_idx = n_terms - 1 - term_idx
                                apply_linear_tensor(decomposed_params[term_idx])
                    elif trotter_mode == "random":
                        for _ in range(trotter_num):
                            for term_idx in np.random.permutation(n_terms):
                                apply_linear_tensor(decomposed_params[term_idx])
        return circ.optimize()

    def transpile(
        self, n, sol_gvars, boxes, edges, *generate_circuit_args, **generate_circuit_kwargs
    ):
        n = len(n) if isinstance(n, list) else n
        if "trotter_args" in generate_circuit_kwargs:
            trotter_args = generate_circuit_kwargs["trotter_args"]
            del generate_circuit_kwargs["trotter_args"]
        else:
            trotter_args = {"order": 1, "mode": 1, "sequential": True, "randomized": False}
        circ = self.generate_circuit(n, *generate_circuit_args, **generate_circuit_kwargs)
        return self.clean_as(n, boxes, edges, circ, trotter_args)


def rotate_to_x(params):
    theta = np.sqrt((params**2).sum())
    params = params / theta
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.array([[1, 0], [0, -1]])
    ham=params[0]*X+params[1]*Y+params[2]*Z
    def get_unitary(theta,phi,lbd):
        return np.array([[np.cos(theta),-np.exp(1j*lbd)*np.sin(theta)],[np.exp(1j*phi)*np.sin(theta),np.exp(1j*(phi+lbd))*np.cos(theta)]])
    def loss(angles):
        unitary=get_unitary(angles[0],angles[1],angles[2])
        u1=unitary.conj().T@sp.linalg.expm(-1j*X)@unitary
        u2=sp.linalg.expm(-1j*ham)
        loss=np.linalg.norm(u1-u2)
        return loss
                            
    res = sp.optimize.minimize(loss, [theta,0,0], bounds=[(0,2*np.pi),(0,2*np.pi),(0,2*np.pi)])
    return get_unitary(res.x[0],res.x[1],res.x[2]),theta

def decompose_ham(thetas, verbose=0):
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
    total_amp = np.sqrt(np.trace(H @ H.conj().T).real)

    def create_tensor_ham(amp, theta1, theta2, phi1, phi2):
        H1 = amp * (
            np.cos(theta1) * X
            + np.sin(theta1) * np.cos(theta2) * Y
            + np.sin(theta1) * np.sin(theta2) * Z
        )
        H2 = amp * (
            np.cos(phi1) * X + np.sin(phi1) * np.cos(phi2) * Y + np.sin(phi1) * np.sin(phi2) * Z
        )
        return np.kron(H1, H2)

    def create_zero_ham():
        return np.zeros((4, 4), dtype=np.complex128)

    def calc_loss(params):
        params = params.reshape(-1, 5)
        n_terms = params.shape[0]
        guess = create_zero_ham()
        commutator_loss = 0
        H_list = []
        for i in range(n_terms):
            H_list.append(
                create_tensor_ham(
                    params[i, 0], params[i, 1], params[i, 2], params[i, 3], params[i, 4]
                )
            )
        for i in range(n_terms):
            for j in range(n_terms):
                if i == j:
                    continue
                else:
                    commutator_loss += np.linalg.norm(H_list[i] @ H_list[j] - H_list[j] @ H_list[i])
            guess += H_list[i]
        return np.linalg.norm(guess - H) + sum(params[:, 0] ** 2) * 0.01 + commutator_loss * 0.01

    def calc_loss_no_penalty(params):
        params = params.reshape(-1, 5)
        n_terms = params.shape[0]
        H_list = []
        guess = create_zero_ham()
        for i in range(n_terms):
            H_list.append(
                create_tensor_ham(
                    params[i, 0], params[i, 1], params[i, 2], params[i, 3], params[i, 4]
                )
            )
        for i in range(n_terms):
            guess += H_list[i]
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

    success = False
    for i in range(1, 4):
        bounds = np.array(
            [(0, 5), (0, np.pi), (0, np.pi), (0, np.pi), (0, np.pi)] * i, dtype=np.float64
        )
        for _ in range(10):
            x0 = np.random.rand(5 * i) * np.array([5, np.pi, np.pi, np.pi, np.pi] * i, dtype=np.float64)
            res = dual_annealing(
                calc_loss,
                bounds=bounds,
                x0=x0,
            )
            if calc_loss_no_penalty(res.x) < 0.01:
                if verbose > 1:
                    print("success with ", i, " terms")
                success = True
                break
        if success:
            break
    if not success:
        raise Exception("decomposition failed")
    return transform_params(res.x)
