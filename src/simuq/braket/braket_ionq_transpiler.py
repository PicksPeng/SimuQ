from simuq.backends.ionq_transpiler import IonQTranspiler
from simuq.backends.ionq_transpiler_2Pauli import IonQTranspiler_2Pauli
from simuq.braket import BraketIonQCircuit


class BraketIonQTranspiler(IonQTranspiler):
    def generate_circuit(self, n: int) -> BraketIonQCircuit:
        return BraketIonQCircuit(n)

class BraketIonQTranspiler_2Pauli(IonQTranspiler_2Pauli):
    def generate_circuit(self, n: int) -> BraketIonQCircuit:
        return BraketIonQCircuit(n)
