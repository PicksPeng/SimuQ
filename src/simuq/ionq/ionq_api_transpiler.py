from simuq.backends.ionq_transpiler import IonQTranspiler
from simuq.backends.ionq_transpiler_2Pauli import IonQTranspiler_2Pauli
from simuq.ionq.ionq_api_circuit import IonQAPICircuit


class IonQAPITranspiler(IonQTranspiler):
    def generate_circuit(
        self, n: int, backend, noise_model=None, name: str = "test"
    ) -> IonQAPICircuit:
        return IonQAPICircuit(n, name, backend, noise_model)

class IonQAPITranspiler_2Pauli(IonQTranspiler_2Pauli):
    def generate_circuit(
        self, n: int, backend, noise_model=None, name: str = "test"
    ) -> IonQAPICircuit:
        return IonQAPICircuit(n, name, backend, noise_model)
