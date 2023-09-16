from qiskit.pulse import DriveChannel, GaussianSquare, ShiftPhase

from simuq.environment import Qubit
from simuq.expression import Expression
from simuq.qmachine import QMachine


def generate_qmachine(backend):
    configuration = backend.configuration()
    defaults = backend.defaults()
    n = configuration.num_qubits
    mach = QMachine()

    connected_pairs_cnot = configuration.coupling_map
    instruction_schedule_map = defaults.instruction_schedule_map

    def get_control_qubit(q1, q2):  # Control performs Z
        cx_sched = instruction_schedule_map.get("cx", qubits=(q1, q2))
        supported = False
        for time, inst in cx_sched.instructions:
            if isinstance(inst.channel, DriveChannel) and not isinstance(inst, ShiftPhase):
                if isinstance(inst.pulse, GaussianSquare):
                    target = inst.channel.index
                    control = q1 if target == q2 else q2
                    supported = True
        if not supported:
            raise ValueError("This machine is not supported!")
        return control

    ql = [Qubit(mach) for i in range(n)]
    link = []
    for item in connected_pairs_cnot:
        if get_control_qubit(item[0], item[1]) == item[0]:
            link.append(item)

    for i in range(n):
        L = mach.add_signal_line()

        ins1 = L.add_instruction("native", "L{}_X_Y".format(i))
        amp = ins1.add_local_variable()
        phase = ins1.add_local_variable()
        ins1.set_ham(amp * (Expression.cos(phase) * ql[i].X + Expression.sin(phase) * ql[i].Y))

        ins2 = L.add_instruction("derived", "L{}_Z".format(i))
        amp = ins2.add_local_variable()
        ins2.set_ham(amp * ql[i].Z)

    for q0, q1 in link:
        L = mach.add_signal_line()

        ins = L.add_instruction("derived", "L{}{}_ZX".format(q0, q1))
        amp = ins.add_local_variable()
        ins.set_ham(amp * ql[q0].Z * ql[q1].X)

        ins = L.add_instruction("derived", "L{}{}_XX".format(q0, q1))
        amp = ins.add_local_variable()
        ins.set_ham(amp * ql[q0].X * ql[q1].X)

        ins = L.add_instruction("derived", "L{}{}_YY".format(q0, q1))
        amp = ins.add_local_variable()
        ins.set_ham(amp * ql[q0].Y * ql[q1].Y)

        ins = L.add_instruction("derived", "L{}{}_ZZ".format(q0, q1))
        amp = ins.add_local_variable()
        ins.set_ham(amp * ql[q0].Z * ql[q1].Z)
    return mach
