from simuq.environment import Qubit
from simuq.expression import Expression
from simuq.qmachine import QMachine


def generate_qmachine():
    mach = QMachine()
    rings_per_row = 5
    qubits_per_ring = 8
    row_count = 2
    qubits = {}
    for row in range(row_count):
        for r in range(rings_per_row):
            for q in range(qubits_per_ring):
                qubit_index = row * 100 + r * 10 + q
                qubits[qubit_index] = Qubit(mach)

    links = []
    for row in range(row_count):
        for r in range(rings_per_row):
            for q in range(qubits_per_ring):
                qubit_index = row * 100 + r * 10 + q
                qubit_neighbor_index = row * 100 + r * 10 + ((q - 1) % qubits_per_ring)
                L = mach.add_signal_line()
                ins1 = L.add_instruction("native", "L{}_X_Y".format(qubit_index))
                amp = ins1.add_local_variable()
                phase = ins1.add_local_variable()
                ins1.set_ham(
                    amp
                    * (
                        Expression.cos(phase) * qubits[qubit_index].X
                        + Expression.sin(phase) * qubits[qubit_index].Y
                    )
                )

                ins2 = L.add_instruction("derived", "L{}_Z".format(qubit_index))
                amp = ins2.add_local_variable()
                ins2.set_ham(amp * qubits[qubit_index].Z)

                # All eight links between qubits in a ring
                links.append((qubit_index, qubit_neighbor_index))

        for r in range(rings_per_row - 1):
            # Links between adjacent rings
            # Note we match indices 1-6 and 2-5, and this is hardcoded, as in an eight-qubit ring topology
            links.append((row * 100 + r * 10 + 1, row * 100 + r * 10 + 16))
            links.append((row * 100 + r * 10 + 2, row * 100 + r * 10 + 15))

    for row in range(row_count - 1):
        for r in range(rings_per_row):
            # Links between adjacent rows
            # Note we match indices 7-4 and 0-3, and this is hardcoded, as in an eight-qubit ring topology
            links.append((row * 100 + r * 10 + 7, row * 100 + r * 10 + 104))
            links.append((row * 100 + r * 10, row * 100 + r * 10 + 103))

    for q0, q1 in links:
        L = mach.add_signal_line()

        ins = L.add_instruction("derived", "L{}{}_XX_YY".format(q0, q1))
        amp = ins.add_local_variable()
        ins.set_ham(amp * (qubits[q0].X * qubits[q1].X + qubits[q0].Y + qubits[q1].Y))

        ins = L.add_instruction("derived", "L{}{}_ZZ".format(q0, q1))
        amp = ins.add_local_variable()
        ins.set_ham(amp * qubits[q0].Z * qubits[q1].Z)

    return mach
