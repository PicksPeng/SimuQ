from simuq.backends.ionq_transpiler import IonQTranspiler
from simuq.ionq.ionq_api_circuit import IonQAPICircuit


class IonQAPITranspiler(IonQTranspiler):
    def generate_circuit(self, n: int, backend, noise_model) -> IonQAPICircuit:
        return IonQAPICircuit(n, "test", backend, noise_model)
