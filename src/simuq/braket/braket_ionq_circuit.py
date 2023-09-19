import numpy as np
from braket.circuits import Circuit, gates

from simuq.backends.ionq_circuit import IonQCircuit, to_turns


class BraketIonQCircuit(IonQCircuit):
    def __init__(self, num_qubits):
        super().__init__(num_qubits)
        self._circuit = Circuit()

    def gpi(self, q, phi):
        self._circuit.gpi(q, phi + self._accum_phases[q])
        return self

    def gpi2(self, q, phi):
        self._circuit.gpi2(q, phi + self._accum_phases[q])
        return self

    def ms_quarter(self, q0, q1, phi0, phi1, theta):
        self._circuit.ms(
            q0,
            q1,
            phi0 + self._accum_phases[q0],
            phi1 + self._accum_phases[q1],
            theta,
        )
        return self

    def optimize(self):
        """
        Optimize 1q gate sequences
        """
        # Not self._circuit.qubit_count because sometimes that's smaller
        new_circ = BraketIonQCircuit(len(self._accum_phases))

        n = len(self._accum_phases)

        unitaries = [np.eye(2)] * n

        for instruction in self._circuit.instructions:
            gate = instruction.operator
            target = instruction.target
            if isinstance(gate, gates.GPi):
                qubit = target[0]
                unitaries[qubit] = np.matmul(
                    self._gpi_mat(to_turns(gate.angle)),
                    unitaries[qubit],
                )
            elif isinstance(gate, gates.GPi2):
                qubit = target[0]
                unitaries[qubit] = np.matmul(
                    self._gpi2_mat(to_turns(gate.angle)),
                    unitaries[qubit],
                )
            elif isinstance(gate, gates.MS):
                for q in target:
                    new_circ._add_unitary(q, unitaries[q])
                    # Reset accumulated matrix
                    unitaries[q] = np.eye(2)
                q0, q1 = target
                phi0 = gate.angle_1
                phi1 = gate.angle_2
                angle = gate.angle_3
                new_circ.ms(q0, q1, phi0, phi1, angle)
            else:
                raise Exception("Unknown gate:", gate["gate"])

        for q in range(n):
            new_circ._add_unitary(q, unitaries[q])
            new_circ.rz(q, -self._accum_phases[q])

        return new_circ

    def add(self, circ):
        """
        Append a circuit behind self
        """

        if len(circ._accum_phases) != len(self._accum_phases):
            raise Exception("Circuit sizes are different.")

        for instruction in circ._circuit.instructions:
            gate = instruction.operator
            target = instruction.target
            if isinstance(gate, gates.GPi):
                qubit = target[0]
                angle = self._accum_phases[qubit] + gate.angle
                self.gpi(qubit, angle)
            elif isinstance(gate, gates.GPi2):
                qubit = target[0]
                angle = self._accum_phases[qubit] + gate.angle
                self.gpi2(qubit, angle)
            elif isinstance(gate, gates.MS):
                q0, q1 = target
                phi0 = self._accum_phases[q0] + gate.angle_1
                phi1 = self._accum_phases[q1] + gate.angle_2
                angle = gate.angle_3
                self.ms(q0, q1, phi0, phi1, angle)
            else:
                raise Exception("Unknown gate:", gate["gate"])

        for q, phi in enumerate(circ._accum_phases):
            self._accum_phases[q] += phi

        return self

    def copy(self):
        """
        Generate a copy of self
        """
        circ = BraketIonQCircuit(len(self._accum_phases))
        circ.add(self)
        return circ

    @property
    def braket_circuit(self):
        return self._circuit
