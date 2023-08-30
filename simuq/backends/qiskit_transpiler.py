from math import atan2, cos, floor, pi, sin, sqrt

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit.gate import Gate
from qiskit.circuit.library.standard_gates import (
    HGate,
    RXGate,
    RXXGate,
    RYYGate,
    RZGate,
    RZXGate,
    RZZGate,
)
from qiskit.dagcircuit import DAGCircuit
from qiskit.pulse import Drag, DriveChannel, GaussianSquare, Play, Schedule, ShiftPhase
from qiskit.pulse.instruction_schedule_map import CalibrationPublisher
from qiskit.transpiler import PassManager
from qiskit.transpiler.basepasses import TransformationPass
from qiskit.transpiler.passes import CommutativeCancellation, RZXCalibrationBuilder
from qiskit.transpiler.passes.calibration.base_builder import CalibrationBuilder


class eYZXGate(Gate):
    def __init__(self, t):
        """Create new eYZX gate."""
        super().__init__("eyzx", 3, [t])

    def _define(self):
        """
        gate eyzx(t) a, b, c
        Evolve under YZX for time t
        """
        # pylint: disable=cyclic-import

        t = self.params[0]
        if t < 0 or t >= 2 * pi:
            t -= floor(t / (2 * pi)) * 2 * pi
        repeat = floor(t / (pi / 2)) + 1
        t = t / repeat
        t1 = atan2(-sqrt(sin(2 * t)), 1) / 2
        t2 = atan2(+sqrt(sin(2 * t)), cos(t) - sin(t)) / 2
        q = QuantumRegister(3, "q")
        qc = QuantumCircuit(q, name=self.name)
        rules = sum(
            [
                [
                    (RYYGate(t1 * 2), [q[0], q[1]], []),
                    (RXXGate(t2 * 2), [q[1], q[2]], []),
                    (RYYGate(t2 * 2), [q[0], q[1]], []),
                    (RXXGate(t1 * 2), [q[1], q[2]], []),
                ]
                for i in range(repeat)
            ],
            [],
        )
        qc._data = rules
        self.definition = qc

    def inverse(self):
        """Return inverse RZX gate (i.e. with the negative rotation angle)."""
        return eYZXGate(-self.params[0])


class RYXGate(Gate):
    def __init__(self, theta):
        """Create new RYX gate."""
        super().__init__("ryx", 2, [theta])

    def _define(self):
        """
        gate ryx(theta) a, b
        Rotate along YX for theta Ryx(theta) = exp(-i theta/2 YX)
        """
        # pylint: disable=cyclic-import

        theta = self.params[0]
        q = QuantumRegister(2, "q")
        qc = QuantumCircuit(q, name=self.name)
        rules = [
            (RZGate(-np.pi / 2), [q[0]], []),
            (HGate(), [q[0]], []),
            (RZXGate(theta), [q[0], q[1]], []),
            (HGate(), [q[0]], []),
            (RZGate(np.pi / 2), [q[0]], []),
        ]
        qc._data = rules
        self.definition = qc

    def inverse(self):
        """Return inverse RZX gate (i.e. with the negative rotation angle)."""
        return RYXGate(-self.params[0])


class RXCalibrationBuilder(CalibrationBuilder):
    # A transpile pass for RX gates
    def __init__(
        self,
        instruction_schedule_map,
        qubit_channel_mapping,
    ):
        super().__init__()
        self._inst_map = instruction_schedule_map
        self._channel_map = qubit_channel_mapping

    def supported(self, node_op, qubits):
        return isinstance(node_op, RXGate)

    @staticmethod
    def rescale_gaussian_inst(instruction, theta):
        pulse_ = instruction.pulse
        amp_scale = (1 + theta / np.pi) % 2 - 1
        return Play(
            Drag(
                amp=pulse_.amp * amp_scale * 2,
                sigma=pulse_.sigma,
                duration=pulse_.duration,
                beta=pulse_.beta,
            ),
            channel=instruction.channel,
        )

    def get_calibration(self, node_op, qubits):
        theta = node_op.params[0]
        qubit = qubits[0]
        x_sched = self._inst_map.get("sx", qubits=(qubit))
        rx_theta = Schedule(name="rx(%.3f)" % theta)
        rx_theta.metadata["publisher"] = CalibrationPublisher.QISKIT

        if theta == 0.0:
            return rx_theta
        inst = x_sched.instructions[0][1]
        x1 = self.rescale_gaussian_inst(inst, theta)
        return rx_theta.append(x1)


class RXXCalibrationBuilder(TransformationPass):
    def __init__(
        self,
        instruction_schedule_map,
        qubit_channel_mapping,
    ):
        super().__init__()
        self._inst_map = instruction_schedule_map
        self._channel_map = qubit_channel_mapping

    def supported(self, node_op) -> bool:
        return isinstance(node_op, RXXGate)

    def get_control_qubit(self, q1, q2):  # Control performs Z
        cx_sched = self._inst_map.get("cx", qubits=(q1, q2))

        for time, inst in cx_sched.instructions:
            if isinstance(inst.channel, DriveChannel) and not isinstance(inst, ShiftPhase):
                if isinstance(inst.pulse, GaussianSquare):
                    target = inst.channel.index
                    control = q1 if target == q2 else q2
        return control

    def run(self, dag: DAGCircuit) -> DAGCircuit:
        qubit_map = {qubit: i for i, qubit in enumerate(dag.qubits)}
        daglist = dag.gate_nodes()
        daglist.reverse()
        for node in daglist:
            if self.supported(node.op):
                qubits = [qubit_map[q] for q in node.qargs]
                control = self.get_control_qubit(qubits[0], qubits[1])
                if control == qubits[0]:
                    control = 0
                else:
                    control = 1
                mini_dag = DAGCircuit()
                p = QuantumRegister(2, "p")
                mini_dag.add_qreg(p)
                mini_dag.apply_operation_back(HGate(), qargs=[p[0]])

                mini_dag.apply_operation_back(RZXGate(node.op.params[0]), qargs=[p[0], p[1]])
                mini_dag.apply_operation_back(HGate(), qargs=[p[0]])

                # substitute the cx node with the above mini-dag
                cx_node = dag.op_nodes(op=RXXGate).pop()
                dag.substitute_node_with_dag(
                    node=cx_node, input_dag=mini_dag, wires=[p[control], p[1 - control]]
                )

        return dag


class RZZCalibrationBuilder(TransformationPass):
    def __init__(
        self,
        instruction_schedule_map,
        qubit_channel_mapping,
    ):
        super().__init__()
        self._inst_map = instruction_schedule_map
        self._channel_map = qubit_channel_mapping

    def supported(self, node_op) -> bool:
        return isinstance(node_op, RZZGate)

    def get_control_qubit(self, q1, q2):  # Control performs Z
        cx_sched = self._inst_map.get("cx", qubits=(q1, q2))

        for time, inst in cx_sched.instructions:
            if isinstance(inst.channel, DriveChannel) and not isinstance(inst, ShiftPhase):
                if isinstance(inst.pulse, GaussianSquare):
                    target = inst.channel.index
                    control = q1 if target == q2 else q2
        return control

    def run(self, dag: DAGCircuit) -> DAGCircuit:
        qubit_map = {qubit: i for i, qubit in enumerate(dag.qubits)}
        daglist = dag.gate_nodes()
        daglist.reverse()
        for node in daglist:
            if self.supported(node.op):
                qubits = [qubit_map[q] for q in node.qargs]
                control = self.get_control_qubit(qubits[0], qubits[1])
                if control == qubits[0]:
                    control = 0
                else:
                    control = 1
                mini_dag = DAGCircuit()
                p = QuantumRegister(2, "p")
                mini_dag.add_qreg(p)
                mini_dag.apply_operation_back(HGate(), qargs=[p[1]])

                mini_dag.apply_operation_back(RZXGate(node.op.params[0]), qargs=[p[0], p[1]])
                mini_dag.apply_operation_back(HGate(), qargs=[p[1]])

                # substitute the cx node with the above mini-dag
                cx_node = dag.op_nodes(op=RZZGate).pop()
                dag.substitute_node_with_dag(
                    node=cx_node, input_dag=mini_dag, wires=[p[control], p[1 - control]]
                )

        return dag


class RYYCalibrationBuilder(TransformationPass):
    def __init__(
        self,
        instruction_schedule_map,
        qubit_channel_mapping,
    ):
        super().__init__()
        self._inst_map = instruction_schedule_map
        self._channel_map = qubit_channel_mapping

    def supported(self, node_op) -> bool:
        return isinstance(node_op, RYYGate)

    def get_control_qubit(self, q1, q2):  # Control performs Z
        cx_sched = self._inst_map.get("cx", qubits=(q1, q2))

        for time, inst in cx_sched.instructions:
            if isinstance(inst.channel, DriveChannel) and not isinstance(inst, ShiftPhase):
                if isinstance(inst.pulse, GaussianSquare):
                    target = inst.channel.index
                    control = q1 if target == q2 else q2
        return control

    def run(self, dag: DAGCircuit) -> DAGCircuit:
        qubit_map = {qubit: i for i, qubit in enumerate(dag.qubits)}
        daglist = dag.gate_nodes()
        daglist.reverse()
        for node in daglist:
            if self.supported(node.op):
                qubits = [qubit_map[q] for q in node.qargs]
                control = self.get_control_qubit(qubits[0], qubits[1])
                if control == qubits[0]:
                    control = 0
                else:
                    control = 1
                mini_dag = DAGCircuit()
                p = QuantumRegister(2, "p")
                mini_dag.add_qreg(p)
                mini_dag.apply_operation_back(RZGate(-np.pi / 2), qargs=[p[0]])
                mini_dag.apply_operation_back(HGate(), qargs=[p[0]])
                mini_dag.apply_operation_back(HGate(), qargs=[p[1]])
                mini_dag.apply_operation_back(RZGate(np.pi / 2), qargs=[p[1]])

                mini_dag.apply_operation_back(RZXGate(node.op.params[0]), qargs=[p[0], p[1]])
                mini_dag.apply_operation_back(RZGate(-np.pi / 2), qargs=[p[1]])
                mini_dag.apply_operation_back(HGate(), qargs=[p[1]])
                mini_dag.apply_operation_back(HGate(), qargs=[p[0]])
                mini_dag.apply_operation_back(RZGate(np.pi / 2), qargs=[p[0]])

                # substitute the cx node with the above mini-dag
                cx_node = dag.op_nodes(op=RYYGate).pop()
                dag.substitute_node_with_dag(
                    node=cx_node, input_dag=mini_dag, wires=[p[control], p[1 - control]]
                )

        return dag


class RYXCalibrationBuilder(TransformationPass):
    def __init__(
        self,
        instruction_schedule_map,
        qubit_channel_mapping,
    ):
        super().__init__()
        self._inst_map = instruction_schedule_map
        self._channel_map = qubit_channel_mapping

    def supported(self, node_op) -> bool:
        return isinstance(node_op, RYXGate)

    def run(self, dag: DAGCircuit) -> DAGCircuit:
        daglist = dag.gate_nodes()
        daglist.reverse()
        for node in daglist:
            if self.supported(node.op):
                mini_dag = DAGCircuit()
                p = QuantumRegister(2, "p")
                mini_dag.add_qreg(p)
                mini_dag.apply_operation_back(RZGate(-np.pi / 2), qargs=[p[0]])
                mini_dag.apply_operation_back(HGate(), qargs=[p[0]])
                mini_dag.apply_operation_back(RZXGate(node.op.params[0]), qargs=[p[0], p[1]])
                mini_dag.apply_operation_back(HGate(), qargs=[p[0]])
                mini_dag.apply_operation_back(RZGate(np.pi / 2), qargs=[p[0]])

                pop_node = dag.op_nodes(op=RYXGate).pop()
                dag.substitute_node_with_dag(node=pop_node, input_dag=mini_dag, wires=p)

        return dag


class eYZXCalibrationBuilder(TransformationPass):
    def __init__(
        self,
        instruction_schedule_map,
        qubit_channel_mapping,
    ):
        super().__init__()
        self._inst_map = instruction_schedule_map
        self._channel_map = qubit_channel_mapping

    def supported(self, node_op) -> bool:
        return isinstance(node_op, eYZXGate)

    def run(self, dag: DAGCircuit) -> DAGCircuit:
        daglist = dag.gate_nodes()
        daglist.reverse()
        for node in daglist:
            if self.supported(node.op):
                mini_dag = DAGCircuit()

                t = node.op.params[0]
                if t < 0 or t >= 2 * pi:
                    t -= floor(t / (2 * pi)) * 2 * pi
                repeat = floor(t / (pi / 2)) + 1
                t = t / repeat
                t1 = atan2(-sqrt(sin(2 * t)), 1) / 2
                t2 = atan2(+sqrt(sin(2 * t)), cos(t) - sin(t)) / 2

                p = QuantumRegister(3, "p")
                mini_dag.add_qreg(p)
                for i in range(repeat):
                    mini_dag.apply_operation_back(RYYGate(t1 * 2), qargs=[p[0], p[1]])
                    mini_dag.apply_operation_back(RXXGate(t2 * 2), qargs=[p[1], p[2]])
                    mini_dag.apply_operation_back(RYYGate(t2 * 2), qargs=[p[0], p[1]])
                    mini_dag.apply_operation_back(RXXGate(t1 * 2), qargs=[p[1], p[2]])

                pop_node = dag.op_nodes(op=eYZXGate).pop()
                dag.substitute_node_with_dag(node=pop_node, input_dag=mini_dag, wires=p)

        return dag


def get_pm(backend, for_braiding=False):
    configuration = backend.configuration()
    defaults = backend.defaults()
    instruction_schedule_map = defaults.instruction_schedule_map
    qubit_channel_map = configuration._qubit_channel_map

    xx_calibrater = RXXCalibrationBuilder(
        instruction_schedule_map=instruction_schedule_map,
        qubit_channel_mapping=qubit_channel_map,
    )

    zz_calibrater = RZZCalibrationBuilder(
        instruction_schedule_map=instruction_schedule_map,
        qubit_channel_mapping=qubit_channel_map,
    )

    yy_calibrater = RYYCalibrationBuilder(
        instruction_schedule_map=instruction_schedule_map,
        qubit_channel_mapping=qubit_channel_map,
    )

    zx_calibrater = RZXCalibrationBuilder(
        instruction_schedule_map=instruction_schedule_map,
        qubit_channel_mapping=qubit_channel_map,
    )

    x_calibrater = RXCalibrationBuilder(
        instruction_schedule_map=instruction_schedule_map,
        qubit_channel_mapping=qubit_channel_map,
    )

    yzx_calibrater = eYZXCalibrationBuilder(
        instruction_schedule_map=instruction_schedule_map,
        qubit_channel_mapping=qubit_channel_map,
    )

    yx_calibrater = RYXCalibrationBuilder(
        instruction_schedule_map=instruction_schedule_map,
        qubit_channel_mapping=qubit_channel_map,
    )

    if for_braiding:
        pm = PassManager(
            [
                yzx_calibrater,
                yx_calibrater,
                zz_calibrater,
                xx_calibrater,
                yy_calibrater,
                CommutativeCancellation(),
                zx_calibrater,
                CommutativeCancellation(),
                # CommutativeCancellation(),
                x_calibrater,
            ]
        )
    else:
        pm = PassManager(
            [
                zz_calibrater,
                xx_calibrater,
                yy_calibrater,
                # CommutativeCancellation(),
                zx_calibrater,
                # CommutativeCancellation(),
                x_calibrater,
            ]
        )
    return pm
