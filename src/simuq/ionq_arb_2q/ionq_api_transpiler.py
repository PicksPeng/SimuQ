from simuq.backends.ionq_transpiler_arb_2q import IonQTranspiler
from simuq.ionq_arb_2q.ionq_api_circuit import IonQAPICircuit


class IonQAPITranspiler(IonQTranspiler):
    def generate_circuit(
        self, n: int, backend, noise_model=None, name: str = "test"
    ) -> IonQAPICircuit:
        return IonQAPICircuit(n, name, backend, noise_model)
