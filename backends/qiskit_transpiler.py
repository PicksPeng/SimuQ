import matplotlib.pyplot as plt
import numpy as np
from qiskit import IBMQ, Aer, QuantumCircuit, QuantumRegister, schedule, transpile
from qiskit.circuit.library.standard_gates import (
    HGate,
    RXGate,
    RXXGate,
    RYYGate,
    RZXGate,
    RZZGate,
    SdgGate,
    SGate,
)
from qiskit.compiler import schedule, transpile
from qiskit.dagcircuit import DAGCircuit
from qiskit.providers.fake_provider import FakeGuadalupe, FakeJakarta
from qiskit.pulse import (
    ControlChannel,
    Drag,
    DriveChannel,
    GaussianSquare,
    Play,
    Schedule,
    ShiftPhase,
)
from qiskit.pulse.instruction_schedule_map import CalibrationPublisher
from qiskit.tools.monitor import job_monitor
from qiskit.transpiler import PassManager
from qiskit.transpiler.basepasses import TransformationPass
from qiskit.transpiler.passes import RZXCalibrationBuilder
from qiskit.transpiler.passes.calibration.builders import CalibrationBuilder


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
            if isinstance(inst.channel, DriveChannel) and not isinstance(
                inst, ShiftPhase
            ):
                if isinstance(inst.pulse, GaussianSquare):
                    target = inst.channel.index
                    control = q1 if target == q2 else q2
        return control

    def run(self, dag: DAGCircuit) -> DAGCircuit:
        qubit_map = {qubit: i for i, qubit in enumerate(dag.qubits)}
        params_list = []
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

                mini_dag.apply_operation_back(
                    RZXGate(node.op.params[0]), qargs=[p[0], p[1]]
                )
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
            if isinstance(inst.channel, DriveChannel) and not isinstance(
                inst, ShiftPhase
            ):
                if isinstance(inst.pulse, GaussianSquare):
                    target = inst.channel.index
                    control = q1 if target == q2 else q2
        return control

    def run(self, dag: DAGCircuit) -> DAGCircuit:
        qubit_map = {qubit: i for i, qubit in enumerate(dag.qubits)}
        params_list = []
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

                mini_dag.apply_operation_back(
                    RZXGate(node.op.params[0]), qargs=[p[0], p[1]]
                )
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
            if isinstance(inst.channel, DriveChannel) and not isinstance(
                inst, ShiftPhase
            ):
                if isinstance(inst.pulse, GaussianSquare):
                    target = inst.channel.index
                    control = q1 if target == q2 else q2
        return control

    def run(self, dag: DAGCircuit) -> DAGCircuit:
        qubit_map = {qubit: i for i, qubit in enumerate(dag.qubits)}
        params_list = []
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
                mini_dag.apply_operation_back(SdgGate(), qargs=[p[0]])
                mini_dag.apply_operation_back(HGate(), qargs=[p[0]])
                mini_dag.apply_operation_back(HGate(), qargs=[p[1]])
                mini_dag.apply_operation_back(SGate(), qargs=[p[1]])

                mini_dag.apply_operation_back(
                    RZXGate(node.op.params[0]), qargs=[p[0], p[1]]
                )
                mini_dag.apply_operation_back(SdgGate(), qargs=[p[1]])
                mini_dag.apply_operation_back(HGate(), qargs=[p[1]])
                mini_dag.apply_operation_back(HGate(), qargs=[p[0]])
                mini_dag.apply_operation_back(SGate(), qargs=[p[0]])

                # substitute the cx node with the above mini-dag
                cx_node = dag.op_nodes(op=RYYGate).pop()
                dag.substitute_node_with_dag(
                    node=cx_node, input_dag=mini_dag, wires=[p[control], p[1 - control]]
                )

        return dag


def get_pm(backend):

    configuration = backend.configuration()
    properties = backend.properties()
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

    pm = PassManager(
        [zz_calibrater, xx_calibrater, yy_calibrater, zx_calibrater, x_calibrater]
    )
    return pm
