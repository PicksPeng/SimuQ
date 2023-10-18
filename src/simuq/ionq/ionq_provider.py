from simuq.aais import heisenberg
from simuq.ionq.ionq_api_transpiler import IonQAPITranspiler
from simuq.provider import BaseProvider
from simuq.solver import generate_as


class IonQProvider(BaseProvider):
    def __init__(self, api_key=None, from_file=None):
        if api_key is None:
            if from_file is None:
                raise Exception("No API_key provided.")
            else:
                with open(from_file, "r") as f:
                    api_key = f.readline().strip()
        self.API_key = api_key
        self.all_backends = ["harmony", "aria-1", "aria-2", "forte"]

        super().__init__()

    def supported_backends(self):
        print(self.all_backends)

    def compile(
        self,
        qs,
        backend="aria-1",
        aais="heisenberg",
        tol=0.01,
        trotter_num=6,
        trotter_mode=1,
        state_prep=None,
        meas_prep=None,
        verbose=0,
    ):
        if backend == "harmony":
            nsite = 11
        elif backend == "aria-1":
            nsite = 23
        elif backend == "aria-2":
            nsite = 23
        elif backend == "forte":
            nsite = 32
        elif isinstance(backend, int):
            if backend > 0:
                nsite = backend
        else:
            raise Exception("Backend is not supported.")

        if qs.num_sites > nsite:
            raise Exception("Device has less sites than the target quantum system.")

        if aais == "heisenberg":
            mach = heisenberg.generate_qmachine(qs.num_sites, e=None)
            comp = IonQAPITranspiler().transpile

        if trotter_mode == "random":
            trotter_args = {"num": trotter_num, "order": 1, "sequential": False}
            randomized = True
        else:
            trotter_args = {"num": trotter_num, "order": trotter_mode, "sequential": True}
            randomized = False

        layout, sol_gvars, boxes, edges = generate_as(
            qs,
            mach,
            trotter_args=trotter_args,
            solver="least_squares",
            solver_args={"tol": tol},
            override_layout=[i for i in range(qs.num_sites)],
            verbose=verbose,
        )
        self.prog = comp(
            qs.num_sites,
            sol_gvars,
            boxes,
            edges,
            randomized=randomized,
            backend="qpu." + backend,
            noise_model=backend,
        )

        if state_prep is not None:
            self.prog = state_prep.copy().add(self.prog, inherit_from_back=True)

        if meas_prep is not None:
            self.prog.add(meas_prep)

        self.prog = self.prog.optimize()

        self.prog = self.prog.job

        self.layout = layout
        self.qs_names = qs.print_sites()

    def print_circuit(self):
        if self.prog is None:
            raise Exception("No compiled job in record.")
        print(self.prog["input"]["circuit"])

    def run(self, shots, on_simulator=False, with_noise=False, verbose=0):
        if self.prog is None:
            raise Exception("No compiled job in record.")

        import json

        import requests

        headers = {
            "Authorization": "apiKey " + self.API_key,
            "Content-Type": "application/json",
        }

        data = self.prog.copy()

        data["shots"] = shots
        if on_simulator:
            data["target"] = "simulator"
            if not with_noise:
                data["noise"] = {"model": "ideal"}

        # print(data)

        response = requests.post(
            "https://api.ionq.co/v0.3/jobs", headers=headers, data=json.dumps(data)
        )
        self.task = response.json()

        if verbose >= 0:
            print(self.task)

        return self.task

    def results(self, job_id=None, verbose=0):
        if job_id is None:
            if self.task is not None:
                job_id = self.task["id"]
            else:
                raise Exception("No submitted job in record.")

        import requests

        headers = {
            "Authorization": "apiKey " + self.API_key,
        }

        response = requests.get("https://api.ionq.co/v0.2/jobs/" + job_id, headers=headers)
        res = response.json()

        if res["status"] != "completed":
            if verbose >= 0:
                print("Job is not completed")
                print(res)
            return None

        def layout_rev(res):
            n = len(self.layout)
            b = IonQProvider.to_bin(res, n)
            ret = ""
            for i in range(n):
                ret += b[self.layout[i]]
            return ret

        def results_from_data(data):
            ret = dict()
            for key in data.keys():
                ret[layout_rev(int(key))[-1::-1]] = data[key]
            return ret

        return results_from_data(res["data"]["histogram"])
