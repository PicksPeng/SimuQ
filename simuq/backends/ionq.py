# Currently this file is not usable:
# AWS Braket does not support arbitrary-angle MS gate for now
# Will be tested later
# If you want to use IonQ's device through braket,
# try qiskit_iontrap.py which generates a circuit in bk_c.

import networkx as nx
import numpy as np


class Circuit:
    def __init__(self, name, n_qubits, backend="simulator", noise_model=None):
        self.job = {
            "lang": "json",
            "shots": 4096,
            "name": name,
            "target": "simulator",
            "body": {
                "gateset": "native",
                "qubits": n_qubits,
                "circuit": [],
            },
        }
        self.job["target"] = backend
        if noise_model is not None:
            self.job["noise"] = {"model": noise_model}

    def gpi2(self, q, phi):
        self.job["body"]["circuit"].append({})
        gate = self.job["body"]["circuit"][-1]
        gate["gate"] = "gpi2"
        gate["target"] = q
        gate["phase"] = (phi / (2 * np.pi)) % 1

    def gpi(self, q, phi):
        self.job["body"]["circuit"].append({})
        gate = self.job["body"]["circuit"][-1]
        gate["gate"] = "gpi"
        gate["target"] = q
        gate["phase"] = (phi / (2 * np.pi)) % 1

    def ms(self, q0, q1, phi0, phi1, theta):
        # assume angle is between 0 and 1
        if 0 <= theta <= np.pi / 2:
            self.job["body"]["circuit"].append({})
            gate = self.job["body"]["circuit"][-1]
            gate["gate"] = "ms"
            gate["targets"] = [q0, q1]
            gate["phases"] = [(phi0 / (2 * np.pi)) % 1, (phi1 / (2 * np.pi)) % 1]
            gate["angle"] = theta / (2 * np.pi)
        elif 3 * np.pi / 2 <= theta <= 2 * np.pi:
            self.job["body"]["circuit"].append({})
            gate = self.job["body"]["circuit"][-1]
            gate["gate"] = "ms"
            gate["targets"] = [q0, q1]
            gate["phases"] = [(phi0 / (2 * np.pi)) % 1, ((phi1 + np.pi) / (2 * np.pi)) % 1]
            gate["angle"] = 1 - theta / (2 * np.pi)
        else:
            raise ValueError(f"Parameter theta is {theta}, must be between 0 and pi/2 or 3*pi/2 and 2*pi (use two gates instead)")
        


def clean_as(n, boxes, edges, backend="simulator", noise_model=None):
    link = [(i, j) for i in range(n) for j in range(i + 1, n)]
    DG = nx.DiGraph()
    DG.add_nodes_from([i for i in range(len(boxes))])
    DG.add_edges_from(edges)
    topo_order = list(nx.topological_sort(DG))
    circ = Circuit("test", n, backend, noise_model)
    accum_phase = [0 for i in range(n)]


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
                        # Rz(q, phi)
                        accum_phase[q] -= phi
                        accum_phase[q] %= (2 * np.pi)
                        # Rx(q, rot)
                        if abs(rot / (2 * np.pi) - 0.25) < 1e-6:
                            circ.gpi2(q, (accum_phase[q] + 0) % (2 * np.pi))
                        elif abs(rot / (2 * np.pi) + 0.25) < 1e-6:
                            circ.gpi2(q, (accum_phase[q] + np.pi) % (2 * np.pi))
                        elif abs(rot / (2 * np.pi) - 0.5) < 1e-6:
                            circ.gpi(q, (accum_phase[q] + 0) % (2 * np.pi))
                        elif abs(rot / (2 * np.pi) + 0.5) < 1e-6:
                            circ.gpi(q, (accum_phase[q] + np.pi) % (2 * np.pi))
                        else:
                            circ.gpi2(q, (accum_phase[q] + 3 * np.pi / 2) % (2 * np.pi))
                            accum_phase[q] -= rot
                            accum_phase[q] %= (2 * np.pi)
                            circ.gpi2(q, (accum_phase[q] + np.pi / 2) % (2 * np.pi))
                        # Rz(q, -phi)
                        accum_phase[q] += phi
                        accum_phase[q] %= (2 * np.pi)

                else:
                    q = line
                    # Rz(q, 2 * params[0] * t)
                    accum_phase[q] -= 2 * params[0] * t
                    accum_phase[q] %= (2 * np.pi)
            else:
                (q0, q1) = link[line - n]
                theta = 2 * params[0] * t
                if ins == 0:
                    # R_XX(theta)
                    if abs(theta) > 1e-5:
                        if (theta / (2 * np.pi)) % 1 <= 0.25 or (theta / (2 * np.pi)) % 1 >= 0.75:
                            circ.ms(q0, q1, accum_phase[q0], accum_phase[q1], theta % (2 * np.pi))
                        elif 0.25 <= (theta / (2 * np.pi)) % 1 <= 0.5:
                            circ.gpi(q0, accum_phase[q0])
                            circ.gpi(q1, accum_phase[q1])
                            circ.ms(q0, q1, (accum_phase[q0] + np.pi) % (2 * np.pi), accum_phase[q1], np.pi - (theta % (2 * np.pi)))
                        elif 0.5 <= (theta / (2 * np.pi)) % 1 <= 0.75:
                            circ.gpi(q0, accum_phase[q0])
                            circ.gpi(q1, accum_phase[q1])
                            circ.ms(q0, q1, accum_phase[q0], accum_phase[q1], (theta % (2 * np.pi)) - np.pi)
                        else:
                            raise ValueError(f"Rotation angle is {theta}, should be between 0 and 2*pi")
                elif ins == 1:
                    # R_YY(theta)
                    if abs(theta) > 1e-5:
                        if (theta / (2 * np.pi)) % 1 <= 0.25 or (theta / (2 * np.pi)) % 1 >= 0.75:
                            circ.ms(q0, q1, (accum_phase[q0] + np.pi / 2) % (2 * np.pi), (accum_phase[q1] + np.pi / 2) % (2 * np.pi), theta % (2 * np.pi))
                        elif 0.25 <= (theta / (2 * np.pi)) % 1 <= 0.5:
                            circ.gpi(q0, accum_phase[q0])
                            circ.gpi(q1, accum_phase[q1])
                            circ.ms(q0, q1, (accum_phase[q0] + np.pi / 2) % (2 * np.pi), (accum_phase[q1] + np.pi / 2) % (2 * np.pi), np.pi - (theta % (2 * np.pi)))
                        elif 0.5 <= (theta / (2 * np.pi)) % 1 <= 0.75:
                            circ.gpi(q0, accum_phase[q0])
                            circ.gpi(q1, accum_phase[q1])
                            circ.ms(q0, q1, (accum_phase[q0] + np.pi / 2) % (2 * np.pi), (accum_phase[q1] + np.pi / 2) % (2 * np.pi), (theta % (2 * np.pi)) - np.pi)
                        else:
                            raise ValueError(f"Rotation angle is {theta}, should be between 0 and 2*pi")
                else:
                    # R_ZZ(theta)
                    if abs(theta) > 1e-5:
                        # R_X(-pi/2)
                        circ.gpi2(q0, (accum_phase[q0] + np.pi) % (2 * np.pi))
                        circ.gpi2(q1, (accum_phase[q1] + np.pi) % (2 * np.pi))
                        # R_YY(theta)
                        if (theta / (2 * np.pi)) % 1 <= 0.25 or (theta / (2 * np.pi)) % 1 >= 0.75:
                            circ.ms(q0, q1, (accum_phase[q0] + np.pi / 2) % (2 * np.pi), (accum_phase[q1] + np.pi / 2) % (2 * np.pi), theta % (2 * np.pi))
                        elif 0.25 <= (theta / (2 * np.pi)) % 1 <= 0.5:
                            circ.gpi(q0, accum_phase[q0])
                            circ.gpi(q1, accum_phase[q1])
                            circ.ms(q0, q1, (accum_phase[q0] + np.pi / 2) % (2 * np.pi), (accum_phase[q1] + np.pi / 2) % (2 * np.pi), np.pi - (theta % (2 * np.pi)))
                        elif 0.5 <= (theta / (2 * np.pi)) % 1 <= 0.75:
                            circ.gpi(q0, accum_phase[q0])
                            circ.gpi(q1, accum_phase[q1])
                            circ.ms(q0, q1, (accum_phase[q0] + np.pi / 2) % (2 * np.pi), (accum_phase[q1] + np.pi / 2) % (2 * np.pi), (theta % (2 * np.pi)) - np.pi)
                        else:
                            raise ValueError(f"Rotation angle is {theta}, should be between 0 and 2*pi")
                        # R_X(pi/2)
                        circ.gpi2(q0, (accum_phase[q0]) % (2 * np.pi))
                        circ.gpi2(q1, (accum_phase[q1]) % (2 * np.pi))
    # circ = Circuit().add_verbatim_box(circ)
    return circ.job, accum_phase


def transpile(n, sol_gvars, boxes, edges, backend="simulator", noise_model=None):
    if isinstance(n, list) :
        n = len(n)
    circ, accum_phase = clean_as(n, boxes, edges, backend, noise_model)
    return circ
