from simuq.backends.ionq_transpiler import IonQTranspiler
from simuq.braket import BraketIonQCircuit


class BraketIonQTranspiler(IonQTranspiler):
    def generate_circuit(self, n: int) -> BraketIonQCircuit:
        return BraketIonQCircuit(n)
