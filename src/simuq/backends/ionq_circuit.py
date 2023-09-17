import cmath
from abc import ABC, abstractmethod

import numpy as np


def to_turns(phi):
    return (phi / (2 * np.pi)) % 1.0


def isclose(a, b, tol=1e-6):
    return cmath.isclose(a, b, abs_tol=tol)


class IonQCircuit(ABC):
    def __init__(self, num_qubits):
        self._accum_phases = [0] * num_qubits

    @abstractmethod
    def gpi(self, q, phi):
        pass

    @abstractmethod
    def gpi2(self, q, phi):
        pass

    def rz(self, q, phi):
        self._accum_phases[q] -= phi
        return self

    def ms(self, q0, q1, phi0, phi1, theta):
        theta = theta % (2 * np.pi)
        if 0 <= theta <= np.pi / 2:
            self.ms_quarter(q0, q1, phi0, phi1, theta)
        elif np.pi / 2 < theta <= np.pi:
            self.gpi(q0, phi0)
            self.gpi(q1, phi1)
            self.ms(q0, q1, phi0 + np.pi, phi1, np.pi - theta)
        elif np.pi < theta <= 3 * np.pi / 2:
            self.gpi(q0, phi0)
            self.gpi(q1, phi1)
            self.ms(q0, q1, phi0, phi1, theta - np.pi)
        elif 3 * np.pi / 2 <= theta <= 2 * np.pi:
            self.ms_quarter(q0, q1, phi0, phi1 + np.pi, 2 * np.pi - theta)
        return self

    @abstractmethod
    def ms_quarter(self, q0, q1, phi0, phi1, theta):
        pass

    @staticmethod
    def _gpi_mat(turns):
        phi = turns * 2 * np.pi
        return np.array([[0, np.exp(-1j * phi)], [np.exp(1j * phi), 0]])

    @staticmethod
    def _gpi2_mat(turns):
        phi = turns * 2 * np.pi
        return np.array([[1, -1j * np.exp(-1j * phi)], [-1j * np.exp(1j * phi), 1]]) / np.sqrt(2)

    @staticmethod
    def _decomp_rz(U):
        if (
            isclose(abs(U[0, 0]), 1)
            and isclose(abs(U[1, 1]), 1)
            and isclose(U[1, 0], 0)
            and isclose(U[0, 1], 0)
        ):
            return np.angle(U[1, 1] / U[0, 0])
        return None

    @staticmethod
    # _GPi(theta)_Rz(phi*2)_ = GPi(theta + phi)
    def _decomp_gpirz(U):
        if isclose(U[0, 0], 0) and isclose(U[1, 1], 0) and isclose(U[0, 1], np.conjugate(U[1, 0])):
            return np.angle(U[1, 0]), 0
        return None

    @staticmethod
    def _decomp_gpi2rz(U):
        invsqrt2 = 1 / np.sqrt(2)
        for i in range(2):
            for j in range(2):
                if not isclose(abs(U[0, 0]), invsqrt2):
                    return None
        lamb = np.angle(U[1, 1] / U[0, 0])
        ph = lamb / 2
        U *= np.exp(1j * (ph - np.angle(U[1, 1])))
        e_iph = np.exp(1j * ph)
        U[0, :] *= e_iph
        U[1, :] *= e_iph.conj()
        if (
            isclose(U[0, 0], invsqrt2)
            and isclose(U[1, 1], invsqrt2)
            and isclose(U[0, 1] / -1j, np.conjugate(U[1, 0] / -1j))
        ):
            return np.angle(U[1, 0] / -1j), lamb
        return None

    @staticmethod
    # _GPi2(psi)_GPi2(phi)_Rz(lamb)_
    def _decomp_gpi2gpi2rz(U):
        if isclose(abs(U[0, 0]), 0):
            # In this case, the decomp has psi=phi, absorbed by GPi(theta)
            return None
        d = np.angle(U[1, 1] / U[0, 0]) / 2
        U *= np.exp(1j * (d - np.angle(U[1, 1])))
        if not (
            isclose(U[0, 0], np.conjugate(U[1, 1]))
            and isclose(U[0, 1] / -1j, np.conjugate(U[1, 0] / -1j))
        ):
            return None
        delta = np.arccos((2 - abs(U[0, 0] * 2) ** 2) / 2)
        ph = np.angle(U[1, 1] / (1 - np.exp(-1j * delta)))
        lamb = ph * 2
        mid = np.angle(U[1, 0] / -1j / np.exp(1j * ph))
        psi = mid + delta / 2
        phi = mid - delta / 2
        if not (
            isclose(U[1, 0], -1j * np.exp(1j * ph) * (np.exp(1j * psi) + np.exp(1j * phi)) / 2)
        ):
            return None
        return psi, phi, lamb

    @staticmethod
    def _decomp_gpi2gpigpi2(U):
        if isclose(abs(U[0, 0]), 0):
            # In this case, the decomp has psi=phi, absorbed by GPi(theta)
            return None
        U = U / (U[0][0] / abs(U[0][0]))
        alpha = np.arccos(U[0][0].real)
        beta = -np.angle(U[1][1])
        gamma = np.angle(U[1][0]) + beta / 2

        theta1 = (beta + 2 * gamma) / 2
        theta3 = (2 * gamma - beta) / 2
        theta2 = gamma - alpha
        return theta1, theta2, theta3
        # GPi2(theta1)@GPi(theta2)@GPi2(theta3)

    def _add_unitary(self, q, U):
        res = self._decomp_rz(U.copy())
        if res is not None:
            theta = res
            self.rz(q, theta)
            return
        res = self._decomp_gpirz(U.copy())
        if res is not None:
            theta, lamb = res
            self.gpi(q, theta)
            self.rz(q, lamb)
            return
        res = self._decomp_gpi2rz(U.copy())
        if res is not None:
            theta, lamb = res
            self.gpi2(q, theta)
            self.rz(q, lamb)
            return
        res = self._decomp_gpi2gpi2rz(U.copy())
        if res is not None:
            psi, phi, lamb = res
            self.gpi2(q, psi)
            self.gpi2(q, phi)
            self.rz(q, lamb)
            return
        res = self._decomp_gpi2gpigpi2(U.copy())
        if res is not None:
            theta1, theta2, theta3 = res
            self.gpi2(q, theta1)
            self.gpi(q, theta2)
            self.gpi2(q, theta3)
            return
        raise Exception("Decomposition failed.")

    @abstractmethod
    def optimize(self):
        pass

    @abstractmethod
    def add(self):
        pass

    @abstractmethod
    def copy(self):
        pass

    # @staticmethod
    # # _GPi(theta)_GPi2(-psi+theta)_ = _GPi2(psi+theta)_GPi(theta)_
    # # _GPi2(psi+theta)_GPi(theta)_Rz(phi*2)_ = _GPi2(psi+theta)_GPi(theta+phi)_
    # # Notice _GPi2(psi+phi)_GPi(phi)_ = _GPi2(phi-psi)_Rz(-psi*2)_  # coefficients may be wrong
    # # The following two can be passed, absorbed by gpi2rz
    # def _decomp_gpi2gpirz(U):
    #     """
    #     invsqrt2 = 1 / np.sqrt(2)
    #     for i in range(2) :
    #         for j in range(2) :
    #             if not is_close(abs(U[0, 0]), invsqrt2) :
    #                 return None
    #     phi = np.angle(U[1, 0] / U[0, 1]) / 2
    #     U *= np.exp(1j * (phi - np.angle(U[1, 0])))
    #     if is_close(U[0, 0] / -1j, np.conjugate(U[1, 1] / -1j)):
    #         psi = np.angle(U[0, 0] / U[1, 1]) / 2
    #         return psi + phi, phi, 0
    #     """
    #     return None
    #
    # @staticmethod
    # def decomp_gpigpi2rz(U):
    #     return None


"""
    def to_matrix(self):
        from qiskit import QuantumCircuit
        from qiskit.quantum_info import Operator

        n = self.job["body"]["qubits"]
        circ = QuantumCircuit(n)

        from qiskit.circuit.library import RZGate
        from qiskit_ionq import GPI2Gate, GPIGate, IonQProvider, MSGate

        for ind, gate in enumerate(self.job["body"]["circuit"]):
            if gate["gate"] == "gpi":
                circ.append(GPIGate(gate["phase"]), [n - 1 - gate["target"]])
            elif gate["gate"] == "gpi2":
                circ.append(GPI2Gate(gate["phase"]), [n - 1 - gate["target"]])
            elif gate["gate"] == "ms":
                circ.append(MSGate(gate["phases"][1], gate["phases"][0], gate["angle"]), list(n - 1 - np.array(gate["targets"])))
            else:
                raise Exception("Unknown gate:", gate["gate"])
        for q in range(self.job["body"]["qubits"]):
            circ.append(RZGate(-self.accum_phase[q]), [n - 1 - q])

        return Operator(circ).to_matrix()

def random_circ(n, m) :
    if n <= 1 :
        raise Exception("n must be > 1.")
    circ = Circuit("random circ", n)
    for _ in range(m) :
        g = np.random.randint(4)
        if g == 0 :
            q = np.random.randint(n)
            theta = np.random.uniform(-3 * np.pi, 3 * np.pi)
            circ.rz(q, theta)
        elif g == 1 :
            q = np.random.randint(n)
            theta = np.random.uniform(-3 * np.pi, 3 * np.pi)
            circ.gpi(q, theta)
        elif g == 2 :
            q = np.random.randint(n)
            theta = np.random.uniform(-3 * np.pi, 3 * np.pi)
            circ.gpi2(q, theta)
        else :
            q1 = np.random.randint(n)
            q2 = np.random.randint(n)
            while q1 == q2 :
                q2 = np.random.randint(n)
            phi1 = np.random.uniform(-3 * np.pi, 3 * np.pi)
            phi2 = np.random.uniform(-3 * np.pi, 3 * np.pi)
            ang = np.random.uniform(-3 * np.pi, 3 * np.pi)
            circ.ms(q1, q2, phi1, phi2, ang)
    return circ
        

TODO: convert to unit tests


def test_1q(T):
    from scipy.stats import unitary_group

    for _ in range(T):
        U = unitary_group.rvs(2)
        circ = Circuit("test", 1)
        circ.add_unitary(0, U)
        Uc = circ.to_matrix()
        Udag_Uc = np.matmul(U.conj().transpose(), Uc)
        Utar = np.array([[1, 0], [0, 1]]) * Udag_Uc[0, 0]
        if not is_close(np.linalg.norm(Udag_Uc - Utar), 0):
            print("Test failed")
            print(U)
            print(Uc)
            print(Udag_Uc)
            return circ

    print("Test passed")

def test_rand_circ(n, m, T):
    for _ in range(T):
        circ1 = random_circ(n, m)
        U1 = circ1.to_matrix()
        circ2 = circ1.optimize()
        U2 = circ2.to_matrix()
        identity = np.zeros((1<<n, 1<<n))
        for i in range(1<<n):
            identity[i, i] = 1
        
        U1dag_U2 = np.matmul(U1.conj().transpose(), U2)
        Utar = np.array(identity * U1dag_U2[0, 0])
        if not is_close(np.linalg.norm(U1dag_U2 - Utar), 0, 1e-6 * (1<<n)):
            print("Test failed")
            return circ1

    print("Test passed")

if __name__ == "__main__":
    test_1q(100)
    test_rand_circ(3, 30, 100)
"""
