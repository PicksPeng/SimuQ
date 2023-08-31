import numpy as np


def to_turns(phi) :
    return (phi / (2 * np.pi)) % 1

def is_close(a, b, tol = 1e-6) :
    return abs(a - b) < tol

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
        self.accum_phase = [0 for _ in range(n_qubits)]

    def gpi2(self, q, phi):
        self.job["body"]["circuit"].append({})
        gate = self.job["body"]["circuit"][-1]
        gate["gate"] = "gpi2"
        gate["target"] = q
        gate["phase"] = to_turns(phi + self.accum_phase[q])

    def gpi(self, q, phi):
        self.job["body"]["circuit"].append({})
        gate = self.job["body"]["circuit"][-1]
        gate["gate"] = "gpi"
        gate["target"] = q
        gate["phase"] = to_turns(phi + self.accum_phase[q])
        
    def rz(self, q, phi) :
        self.accum_phase[q] -= phi

    def ms(self, q0, q1, phi0, phi1, theta):
        theta = theta % (2 * np.pi)
        if 0 <= theta <= np.pi / 2:
            self.job["body"]["circuit"].append({})
            gate = self.job["body"]["circuit"][-1]
            gate["gate"] = "ms"
            gate["targets"] = [q0, q1]
            gate["phases"] = [to_turns(phi0 + self.accum_phase[q0]), to_turns(phi1 + self.accum_phase[q1])]
            gate["angle"] = to_turns(theta / (2 * np.pi))
        elif np.pi / 2 < theta <= np.pi:
            self.gpi(q0, phi0)
            self.gpi(q1, phi1)
            self.ms(q0, q1, phi0 + np.pi, phi1, np.pi - theta)
        elif np.pi < theta <= 3 * np.pi / 2:
            self.gpi(q0, phi0)
            self.gpi(q1, phi1)
            self.ms(q0, q1, phi0, phi1, theta - np.pi)
        elif 3 * np.pi / 2 <= theta <= 2 * np.pi:
            self.job["body"]["circuit"].append({})
            gate = self.job["body"]["circuit"][-1]
            gate["gate"] = "ms"
            gate["targets"] = [q0, q1]
            gate["phases"] = [to_turns(phi0 + self.accum_phase[q0]), to_turns(phi1 + np.pi + self.accum_phase[q1])]
            gate["angle"] = to_turns(2 * pi - theta)

    @staticmethod
    def gpi_mat(turns):
        phi = turns * 2 * np.pi
        return np.array([[0, np.exp(-1j * phi)], [np.exp(1j * phi), 0]])

    @staticmethod
    def gpi2_mat(turns):
        phi = turns * 2 * np.pi
        return np.array([[1, -1j * np.exp(-1j * phi)], [-1j * np.exp(1j * phi), 1]]) / np.sqrt(2)

    @staticmethod
    def decomp_rz(U) :
        if is_close(abs(U[0, 0]), 1) and is_close(abs(U[1, 1]), 1) and is_close(U[1, 0], 0) and is_close(U[0, 1], 0):
            return np.angle(U[1, 1] / U[0, 0])
        return None

    @staticmethod
    # _GPi(theta)_Rz(phi*2)_ = GPi(theta + phi)
    def decomp_gpirz(U) :
        if is_close(U[0, 0], 0) and is_close(U[1, 1], 0) and is_close(U[0, 1], np.conjugate(U[1, 0])):
            return np.angle(U[1, 0]), 0
        return None
    
    @staticmethod
    def decomp_gpi2rz(U) :
        invsqrt2 = 1 / np.sqrt(2)
        for i in range(2) :
            for j in range(2) :
                if not is_close(abs(U[0, 0]), invsqrt2) :
                    return None
        lamb = np.angle(U[1, 1] / U[0, 0])
        ph = lamb / 2
        U *= np.exp(1j * (ph - np.angle(U[1, 1])))
        if is_close(U[0, 0], invsqrt2) and is_close(U[1, 1], invsqrt2) and is_close(U[0, 1] / -1j, np.conjugate(U[1, 0] / -1j)):
            return np.angle(U[1, 0] / -1j) - ph, lamb
        return None

    @staticmethod
    # _GPi(theta)_GPi2(-psi+theta)_ = _GPi2(psi+theta)_GPi(theta)_
    # _GPi2(psi+theta)_GPi(theta)_Rz(phi*2)_ = _GPi2(psi+theta)_GPi(theta+phi)_
    def decomp_gpi2gpirz(U) :
        invsqrt2 = 1 / np.sqrt(2)
        for i in range(2) :
            for j in range(2) :
                if not is_close(abs(U[0, 0]), invsqrt2) :
                    return None
        phi = np.angle(U[1, 0] / U[0, 1]) / 2
        U *= np.exp(1j * (phi - np.angle(U[1, 0])))
        if is_close(U[0, 0] / -1j, np.conjugate(U[1, 1] / -1j)):
            psi = np.angle(U[0, 0] / U[1, 1]) / 2
            return psi + phi, phi, 0
        return None

    @staticmethod
    # _GPi2(psi)_GPi2(phi)_Rz(lamb)_
    def decomp_gpi2gpi2rz(U) :
        if is_close(abs(U[0, 0]), 0) :
            # In this case, the decomp has psi=phi, absorbed by GPi(theta)
            return None
        d = np.angle(U[1, 1] / U[0, 0]) / 2
        U *= np.exp(1j * (d - np.angle(U[1, 1])))
        if not (is_close(U[0, 0], np.conjugate(U[1, 1])) and is_close(U[0, 1] / -1j, np.conjugate(U[1, 0] / -1j))):
            return None
        delta = np.arccos((2 - abs(U[0, 0]) ** 2) / 2)
        ph = np.angle(U[1, 1] / (1 - np.exp(-1j * delta)))
        lamb = ph * 2
        mid = np.angle(U[1, 0] / -1j / np.exp(1j * ph))
        psi = mid + delta / 2
        phi = mid - delta / 2
        if not (is_close(U[1, 0], -1j * np.exp(1j * ph) * (np.exp(1j * psi) + np.exp(1j * phi)))):
            return None
        return psi, phi, lamb

    def decomp_gpi2gpigpi2(U) :
        if is_close(abs(U[0, 0]), 0) :
            # In this case, the decomp has psi=phi, absorbed by GPi(theta)
            return None
        U=U/(U[0][0]/abs(U[0][0]))
        alpha=np.arccos(U[0][0])
        beta=np.angle(U[1][1])
        gamma=np.angle(U[1][0])-beta/2

        theta1=(beta+2*gamma)/2
        theta3=(2*gamma-beta)/2
        theta2=gamma-alpha
        return theta1,theta2,theta3
        # GPi2(theta1)@GPi(theta2)@GPi2(theta3)
    
    @staticmethod
    def decomp_unitary(q, U, circ) :
        res = decomp_rz(U)
        if res != None :
            theta = res
            circ.rz(q, theta)
            return
        res = decomp_gpirz(U)
        if res != None :
            theta, lamb = res
            circ.gpi(q, theta)
            circ.rz(q, lamb)
            return
        res = decomp_gpi2rz(U)
        if res != None :
            theta, lamb = res
            circ.gpi2(q, theta)
            circ.rz(q, lamb)
            return
        res = decomp_gpi2gpirz(U)
        if res != None :
            x, y, z = res
            circ.gpi2(q, x)
            circ.gpi(q, y)
            circ.rz(q, z)
            return
        res = decomp_gpi2gpi2rz(U)
        if res != None :
            psi, phi, lamb = res
            circ.gpi2(q, psi)
            circ.gpi2(q, phi)
            circ.rz(q, lamb)
            return

    def optimize(self) :
        new_circ = Circuit(
            self.job["name"],
            self.job["body"]["qubits"],
            self.job["target"],
            self.job["noise"]["model"],
        )
        
        identity = np.arrays([[1, 0], [0, 1]])
        n = self.job["body"]["qubits"]

        unitaries = [identity.copy() for _ in range(n)]

        for ind, gate in enumerate(self.job["body"]["circuit"]):
            if gate["gate"] == "gpi":
                unitaries[gate["target"]] = np.matmul(
                    self.gpi_mat(gate["angle"]),
                    unitaries[gate["target"]],
                )
            elif gate["gate"] == "gpi2":
                unitaries[gate["target"]] = np.matmul(
                    self.gpi2_mat(gate["angle"]),
                    unitaries[gate["target"]],
                )
            elif gate["gate"] == "ms":
                for q in gate["targets"] :
                    self.decomp_unitary(q, unitaries[q], new_circ)
                q0, q1 = gate["targets"]
                phi0, phi1 = gate["phases"]
                phi0 *= 2 * np.pi
                phi1 *= 2 * np.pi
                angle = gate["angle"] * 2 * np.pi
                new_circ.ms(q0, q1, phi0, phi1, angle)
            else :
                raise Exception("Unknown gate:", gate["gate"])
        for q in range(n) :
            self.decomp_unitary(q, unitaries[q], new_circ)
            new_circ.accum_phase[q] += self.accum_phase[q]

        return new_circ
