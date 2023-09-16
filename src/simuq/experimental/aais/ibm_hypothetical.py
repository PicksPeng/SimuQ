from qiskit.pulse import DriveChannel, GaussianSquare, ShiftPhase

from simuq.environment import Qubit
from simuq.qmachine import *


def get_mach(backend):
    configuration = backend.configuration()
    defaults = backend.defaults()
    n = configuration.num_qubits
    mach = QMachine()

    connected_pairs_cnot = configuration.coupling_map
    instruction_schedule_map = defaults.instruction_schedule_map

    def get_control_qubit(q1, q2):  # Control performs Z
        cx_sched = instruction_schedule_map.get("cx", qubits=(q1, q2))

        for time, inst in cx_sched.instructions:
            if isinstance(inst.channel, DriveChannel) and not isinstance(inst, ShiftPhase):
                if isinstance(inst.pulse, GaussianSquare):
                    target = inst.channel.index
                    control = q1 if target == q2 else q2
        return control

    ql = [Qubit(mach) for i in range(n)]
    link = []
    for item in connected_pairs_cnot:
        if get_control_qubit(item[0], item[1]) == item[0]:
            link.append(item)

    for i in range(n):
        L = SignalLine(mach)

        ins1 = Instruction(L, "native", "L{}_X_Y".format(i))
        amp = LocalVar(ins1)
        phase = LocalVar(ins1)
        ins1.set_ham(amp * (Expression.cos(phase) * ql[i].X + Expression.sin(phase) * ql[i].Y))

        ins2 = Instruction(L, "derived", "L{}_Z".format(i))
        amp = LocalVar(ins2)
        ins2.set_ham(amp * ql[i].Z)

    for q0, q1 in link:
        L = SignalLine(mach)

        ins = Instruction(L, "native", "L{}{}_ZX".format(q0, q1))
        amp = LocalVar(ins)
        ins.set_ham(amp * ql[q0].Z * ql[q1].X)

        ins = Instruction(L, "native", "L{}{}_XX".format(q0, q1))
        amp = LocalVar(ins)
        ins.set_ham(amp * ql[q0].X * ql[q1].X)

        ins = Instruction(L, "native", "L{}{}_YY".format(q0, q1))
        amp = LocalVar(ins)
        ins.set_ham(amp * ql[q0].Y * ql[q1].Y)

        ins = Instruction(L, "native", "L{}{}_ZZ".format(q0, q1))
        amp = LocalVar(ins)
        ins.set_ham(amp * ql[q0].Z * ql[q1].Z)
    return mach
