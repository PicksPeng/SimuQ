import numpy as np
from braket.circuits import Circuit, gates

from simuq.backends.ionq_circuit import IonQCircuit, to_turns


class BraketIonQCircuit(IonQCircuit):
    def __init__(self, num_qubits):
        super().__init__(num_qubits)
        self._circuit = Circuit()

    def gpi(self, q, phi):
        self._circuit.gpi(q, to_turns(phi + self._accum_phases[q]))

    def gpi2(self, q, phi):
        self._circuit.gpi2(q, to_turns(phi + self._accum_phases[q]))

    def ms_quarter(self, q0, q1, phi0, phi1, theta):
        self._circuit.ms(
            q0,
            q1,
            to_turns(phi0 + self._accum_phases[q0]),
            to_turns(phi0 + self._accum_phases[q1]),
            to_turns(theta),
        )

    def optimize(self):
        """
        Optimize 1q gate sequences
        """
        new_circ = BraketIonQCircuit(self._circuit.qubit_count)

        n = len(self._accum_phases)

        unitaries = [np.eye(2)] * n

        for instruction in self._circuit.instructions:
            gate = instruction.operator
            target = instruction.target
            if isinstance(gate, gates.GPi):
                qubit = target[0]
                unitaries[qubit] = np.matmul(
                    self._gpi_mat(gate.angle),
                    unitaries[qubit],
                )
            elif isinstance(gate, gates.GPi2):
                qubit = target[0]
                unitaries[qubit] = np.matmul(
                    self._gpi2_mat(gate.angle),
                    unitaries[qubit],
                )
            elif isinstance(gate, gates.MS):
                for q in target:
                    new_circ._add_unitary(q, unitaries[q])
                    # Reset accumulated matrix
                    unitaries[q] = np.eye(2)
                q0, q1 = target
                phi0 = gate.angle_1 * 2 * np.pi
                phi1 = gate.angle_2 * 2 * np.pi
                angle = gate.angle_3 * 2 * np.pi
                new_circ.ms(q0, q1, phi0, phi1, angle)
            else:
                raise Exception("Unknown gate:", gate["gate"])

        for q in range(n):
            new_circ._add_unitary(q, unitaries[q])
            new_circ.rz(q, -self._accum_phases[q])  # TODO: is this really necessary?

        return new_circ
