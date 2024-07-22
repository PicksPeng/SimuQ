import qiskit
from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime import QiskitRuntimeService

from simuq.provider import BaseProvider
from simuq.solver import generate_as


class IBMProvider(BaseProvider):
    def __init__(self, api_key=None, from_file=None):
        if from_file is not None:
            with open(from_file, "r") as f:
                api_key = f.readline().strip()
        self.api_key = api_key
        QiskitRuntimeService.save_account(channel="ibm_quantum", token=api_key, overwrite=True)
        self.provider = QiskitRuntimeService()
        super().__init__()

    def supported_backends(self):

        print(self.provider.backends())

    def compile(
        self,
        qs,
        backend="ibmq_jakarta",
        aais="heisenberg",
        tol=0.01,
        trotter_num=6,
        verbose=0,
        use_pulse=True,
        state_prep=None,
        with_measure=True,
    ):

        self.backend = self.provider.backend(backend)
        nsite = self.backend.num_qubits

        if qs.num_sites > nsite:
            raise Exception("Device has less sites than the target quantum system.")

        if aais == "heisenberg":
            from simuq.aais import ibm
            from simuq.ibm.qiskit_pulse_ibm import transpile

            mach = ibm.generate_qmachine(self.backend)
            comp = transpile

        layout, sol_gvars, boxes, edges = generate_as(
            qs,
            mach,
            trotter_num,
            solver="least_squares",
            solver_args={"tol": tol},
            override_layout=None,
            verbose=verbose,
        )
        self.prog = comp(
            self.backend,
            layout,
            sol_gvars,
            boxes,
            edges,
            use_pulse=use_pulse,
            with_measure=with_measure,
        )
        from qiskit import transpile as transpile_qiskit

        self.prog = transpile_qiskit(self.prog, backend=self.backend)
        self.layout = layout
        self.qs_names = qs.print_sites()
        if state_prep is not None:
            self.prog = self.prog.compose(state_prep, qubits=layout, front=True)

    def run(self, shots=4096, on_simulator=False, verbose=0):
        if on_simulator:
            self.backend_sim = AerSimulator(method="statevector")
            job = self.backend_sim.run(self.prog, shots=shots)
        else:
            job = self.backend.run(self.prog, shots=shots)
        self.task = job
        if verbose >= 0:
            print(self.task)

    def results(self, job_id=None, on_simulator=False):
        if job_id is None:
            if self.task is not None:
                job_id = self.task.job_id()
            else:
                raise Exception("No submitted job in record.")
        if on_simulator:
            job = self.task
        else:
            job = self.provider.job(job_id)
        status = job.status()
        if status.name == "QUEUED":
            print("Job is not completed")
            return

        def layout_rev(res):
            n = len(self.layout)
            # print(self.layout)
            b = res
            ret = ""
            for i in range(n):
                ret += b[-1 - self.layout[i]]
            return ret

        def results_from_data(data):
            ret = dict()
            for key in data.keys():
                new_key = layout_rev(key)
                if new_key in ret:
                    ret[new_key] += data[key] / n_shots
                else:
                    ret[new_key] = data[key] / n_shots
            return ret

        count = job.result().get_counts()
        n_shots = sum(count.values())
        return results_from_data(count)
