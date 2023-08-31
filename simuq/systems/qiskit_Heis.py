import warnings

import numpy as np
from qiskit import IBMQ, QuantumCircuit, QuantumRegister, schedule, transpile
from qiskit.providers.aer import PulseSimulator, QasmSimulator
from qiskit.providers.aer.pulse import PulseSystemModel

warnings.filterwarnings("ignore")

if IBMQ.active_account() is None:
    IBMQ.load_account()

provider = IBMQ.get_provider(
    hub="ibm-q-community", group="ibmquantumawards", project="open-science-22"
)
backend = provider.get_backend("ibmq_jakarta")
sim_noisy_jakarta = QasmSimulator.from_backend(backend)
# sim_backend = PulseSimulator()
# backend_model = PulseSystemModel.from_backend(backend)
# sim_backend.set_options(system_model=backend_model)

n = 7
link = [(0, 1), (1, 2), (1, 3), (3, 5), (4, 5), (5, 6)]
T = np.pi
n_step = 1

qr = QuantumRegister(7)
qc = QuantumCircuit(qr)
for i in range(n_step):
    for q0, q1 in link:
        qc.rxx(T / n_step, q0, q1)
        qc.ryy(T / n_step, q0, q1)
        qc.rzz(T / n_step, q0, q1)

sch = schedule(transpile(qc, backend), backend)
