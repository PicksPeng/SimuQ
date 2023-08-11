from simuq.aais.qmachine_factory import QMachineFactory
from simuq.qmachine import *
from simuq.expression import Expression


class RigettiQMachineFactory(QMachineFactory):
    @staticmethod
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
                    qubits[qubit_index] = qubit(mach)

        links = []
        for row in range(row_count):
            for r in range(rings_per_row):
                for q in range(qubits_per_ring):
                    qubit_index = row * 100 + r * 10 + q
                    qubit_neighbor_index = row * 100 + r * 10 + ((q - 1) % qubits_per_ring)
                    L = SignalLine(mach)
                    ins1 = Instruction(L, 'native', 'L{}_X_Y'.format(qubit_index))
                    amp = LocalVar(ins1)
                    phase = LocalVar(ins1)
                    ins1.set_ham(amp * (Expression.cos(phase) * qubits[qubit_index].X + Expression.sin(phase) * qubits[
                        qubit_index].Y))

                    ins2 = Instruction(L, 'derived', 'L{}_Z'.format(qubit_index))
                    amp = LocalVar(ins2)
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

        for (q0, q1) in links:
            L = SignalLine(mach)

            ins = Instruction(L, 'derived', 'L{}{}_XX_YY'.format(q0, q1))
            amp = LocalVar(ins)
            ins.set_ham(amp * (qubits[q0].X * qubits[q1].X + qubits[q0].Y + qubits[q1].Y))

            ins = Instruction(L, 'derived', 'L{}{}_ZZ'.format(q0, q1))
            amp = LocalVar(ins)
            ins.set_ham(amp * qubits[q0].Z * qubits[q1].Z)

        return mach
