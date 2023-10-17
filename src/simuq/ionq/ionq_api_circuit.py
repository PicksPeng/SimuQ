import numpy as np

from simuq.backends.ionq_circuit import IonQCircuit, to_turns


class IonQAPICircuit(IonQCircuit):
    def __init__(self, n_qubits, name="circuit", backend="simulator", noise_model=None):
        super().__init__(n_qubits)
        # Circuit are stored in turns, accum_phases are stored in rads.
        self.job = {
            "lang": "json",
            "name": name,
            "target": backend,
            "input": {
                "gateset": "native",
                "qubits": n_qubits,
                "circuit": [],
            },
        }
        if noise_model:
            self.job["noise"] = {"model": noise_model}

    def gpi(self, q, phi):
        self.job["input"]["circuit"].append(
            {"gate": "gpi", "target": q, "phase": to_turns(phi + self._accum_phases[q])}
        )
        return self

    def gpi2(self, q, phi):
        self.job["input"]["circuit"].append(
            {"gate": "gpi2", "target": q, "phase": to_turns(phi + self._accum_phases[q])}
        )
        return self

    def ms_quarter(self, q0, q1, phi0, phi1, theta):
        self.job["input"]["circuit"].append(
            {
                "gate": "ms",
                "targets": [q0, q1],
                "phases": [
                    to_turns(phi0 + self._accum_phases[q0]),
                    to_turns(phi1 + self._accum_phases[q1]),
                ],
                "angle": to_turns(theta),
            }
        )
        return self

    def optimize(self):
        """
        Optimize 1q gate sequences
        """
        if "noise" in self.job:
            new_circ = IonQAPICircuit(
                self.job["input"]["qubits"],
                self.job["name"],
                self.job["target"],
                self.job["noise"]["model"],
            )
        else:
            new_circ = IonQAPICircuit(
                self.job["input"]["qubits"],
                self.job["name"],
                self.job["target"],
            )

        n = len(self._accum_phases)

        unitaries = [np.eye(2)] * n

        for gate in self.job["input"]["circuit"]:
            if gate["gate"] == "gpi":
                unitaries[gate["target"]] = np.matmul(
                    self._gpi_mat(gate["phase"]),
                    unitaries[gate["target"]],
                )
            elif gate["gate"] == "gpi2":
                unitaries[gate["target"]] = np.matmul(
                    self._gpi2_mat(gate["phase"]),
                    unitaries[gate["target"]],
                )
            elif gate["gate"] == "ms":
                for q in gate["targets"]:
                    new_circ._add_unitary(q, unitaries[q])
                    # Reset accumulated matrix
                    unitaries[q] = np.eye(2)
                q0, q1 = gate["targets"]
                phi0, phi1 = np.array(gate["phases"]) * 2 * np.pi
                angle = gate["angle"] * 2 * np.pi
                new_circ.ms(q0, q1, phi0, phi1, angle)
            else:
                raise Exception("Unknown gate:", gate["gate"])

        for q in range(n):
            new_circ._add_unitary(q, unitaries[q])
            new_circ.rz(q, -self._accum_phases[q])

        return new_circ

    def add(self, circ, inherit_from_back=False):
        """
        Append a circuit behind self
        """
        if len(circ._accum_phases) != len(self._accum_phases):
            raise Exception("Circuit sizes are different.")

        for gate in circ.job["input"]["circuit"]:
            if gate["gate"] == "gpi":
                qubit = gate["target"]
                angle = gate["phase"] * 2 * np.pi
                self.gpi(qubit, angle)
            elif gate["gate"] == "gpi2":
                qubit = gate["target"]
                angle = gate["phase"] * 2 * np.pi
                self.gpi2(qubit, angle)
            elif gate["gate"] == "ms":
                q0, q1 = gate["targets"]
                phi0 = gate["phases"][0] * 2 * np.pi
                phi1 = gate["phases"][1] * 2 * np.pi
                angle = gate["angle"] * 2 * np.pi
                self.ms(q0, q1, phi0, phi1, angle)
            else:
                raise Exception("Unknown gate:", gate["gate"])

        for q, phi in enumerate(circ._accum_phases):
            self._accum_phases[q] += phi
        if inherit_from_back:
            self.job["target"] = circ.job["target"]
            if "noise" in circ.job:
                self.job["noise"] = circ.job["noise"]
        return self

    def copy(self):
        """
        Generate a copy of self
        """
        circ = IonQAPICircuit(len(self._accum_phases))
        circ.add(self)
        return circ
