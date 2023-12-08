from simuq.environment import Qubit
from simuq.expression import Expression
from simuq.qmachine import QMachine


def generate_qmachine(n=3):
    mach = QMachine()
    ql = [Qubit(mach) for i in range(n)]
    link = [(i, j) for i in range(n) for j in range(i + 1, n)]

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

        ins = L.add_instruction("derived", "L{}{}_YX".format(q0, q1))
        amps = [ins.add_local_variable() for i in range(9)]

        ins.set_ham(
            amps[0] * ql[q0].X * ql[q1].X
            + amps[1] * ql[q0].X * ql[q1].Y
            + amps[2] * ql[q0].X * ql[q1].Z
            + amps[3] * ql[q0].Y * ql[q1].X
            + amps[4] * ql[q0].Y * ql[q1].Y
            + amps[5] * ql[q0].Y * ql[q1].Z
            + amps[6] * ql[q0].Z * ql[q1].X
            + amps[7] * ql[q0].Z * ql[q1].Y
            + amps[8] * ql[q0].Z * ql[q1].Z
        )

    return mach
