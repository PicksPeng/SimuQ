import matplotlib.pyplot as plt
import numpy as np
from braket.ahs.atom_arrangement import SiteType
from braket.aws import AwsDevice, AwsQuantumTask
from braket.circuits import Circuit
from braket.devices import LocalSimulator

from simuq import _version
from simuq.aais import heisenberg, rydberg1d_global, rydberg2d_global
from simuq.braket.braket_ionq_transpiler import BraketIonQTranspiler
from simuq.braket.braket_rydberg_transpiler import BraketRydbergTranspiler
from simuq.provider import BaseProvider
from simuq.solver import generate_as


class BraketProvider(BaseProvider):
    def __init__(self):
        self.backend_aais = dict()
        self.backend_aais[("ionq", "Harmony")] = ["heisenberg"]
        self.backend_aais[("ionq", "Aria-1")] = ["heisenberg"]
        self.backend_aais[("ionq", "Aria-2")] = ["heisenberg"]
        self.backend_aais[("quera", "Aquila")] = [
            "rydberg1d_global",
            "rydberg2d_global",
        ]

        super().__init__()

    def supported_backends(self):
        for comp, dev in self.backend_aais.keys():
            print(
                f"Hardware provider: {comp}  -  Device name: {dev}  -  AAIS supports: {self.backend_aais[(comp, dev)]}"
            )

    def compile(
        self,
        qs,
        provider,
        device,
        aais,
        tol=0.01,
        trotter_num=6,
        state_prep=None,
        meas_prep=None,
        no_main_body=False,
        verbose=0,
    ):
        if (provider, device) not in self.backend_aais.keys():
            print(provider)
            print(device)
            raise Exception("Not supported hardware provider or device.")
        if aais not in self.backend_aais[(provider, device)]:
            raise Exception("Not supported AAIS on this device.")

        self.qs_names = qs.print_sites()

        self.provider = provider
        self.device = device

        if self.provider == "quera":
            nsite = qs.num_sites

            if aais == "rydberg1d_global":
                transpiler = BraketRydbergTranspiler(1)
                mach = rydberg1d_global.generate_qmachine(nsite)
            elif aais == "rydberg2d_global":
                transpiler = BraketRydbergTranspiler(2)
                mach = rydberg2d_global.generate_qmachine(nsite)
            else:
                raise NotImplementedError

            if state_prep is None:
                state_prep = {"times": [], "omega": [], "delta": [], "phi": []}

            if meas_prep is not None:
                raise Exception(
                    "Currently SimuQ does not support measurement preparation pulses for QuEra devices."
                )

            layout, sol_gvars, boxes, edges = generate_as(
                qs,
                mach,
                trotter_num=1,
                solver="least_squares",
                solver_args={"tol": tol, "time_penalty": 1},
                override_layout=[i for i in range(nsite)],
                verbose=verbose,
            )
            self.sol = [layout, sol_gvars, boxes, edges]

            # Only use this when debugging state_preparation pulses
            if no_main_body:
                for j, (b, t) in enumerate(boxes):
                    boxes[j] = (b, 0.01)
            self.ahs_prog = transpiler.transpile(
                sol_gvars, boxes, edges, state_prep=state_prep, verbose=verbose
            )
            self.prog = self.ahs_prog
            self.layout = layout

        elif self.provider == "ionq":
            if device == "Harmony":
                nsite = 11
            elif device == "Aria-1":
                nsite = 25
            elif device == "Aria-2":
                nsite = 25
            elif isinstance(device, int):
                if device > 0:
                    nsite = device
            else:
                raise Exception("Backend is not supported.")

            if qs.num_sites > nsite:
                raise Exception("Device has less sites than the target quantum system.")

            if aais == "heisenberg":
                mach = heisenberg.generate_qmachine(qs.num_sites, e=None)
                transpiler = BraketIonQTranspiler()

            layout, sol_gvars, boxes, edges = generate_as(
                qs,
                mach,
                trotter_num,
                solver="least_squares",
                solver_args={"tol": tol},
                override_layout=[i for i in range(qs.num_sites)],
                verbose=verbose,
            )

            self.prog = transpiler.transpile(qs.num_sites, sol_gvars, boxes, edges)

            if state_prep is not None:
                self.prog = state_prep.copy().add(self.prog)

            if meas_prep is not None:
                self.prog.add(meas_prep)

            self.prog = Circuit().add_verbatim_box(self.prog.braket_circuit)

            self.layout = layout
            self.qs_names = qs.print_sites()

    def visualize_quera(self):
        braket_prog = self.ahs_prog

        fig1 = _show_register(braket_prog.register)
        fig2 = _show_global_drive(braket_prog.hamiltonian)
        plt.show()

    def visualize_circuit(self):
        print(self.prog)

    def visualize(self):
        if self.prog is None:
            raise Exception("No compiled job in record.")
        if self.provider == "quera":
            self.visualize_quera()

    def run(self, shots, on_simulator=False, verbose=0):
        if self.prog is None:
            raise Exception("No compiled job in record.")

        if self.provider == "quera":
            if on_simulator:
                simulator = LocalSimulator("braket_ahs")
                self.task = simulator.run(self.ahs_prog, shots=shots)
                if verbose >= 0:
                    print("Submitted.")
            else:
                aquila_qpu = AwsDevice("arn:aws:braket:us-east-1::device/qpu/quera/" + self.device)
                user_agent = f"SimuQ/{_version.__version__}"
                aquila_qpu.aws_session.add_braket_user_agent(user_agent)
                discretized_ahs_program = self.ahs_prog.discretize(aquila_qpu)
                self.task = aquila_qpu.run(discretized_ahs_program, shots=shots)
                meta = self.task.metadata()
                if verbose >= 0:
                    print("Submitted.")
                    print("Task arn: ", meta["quantumTaskArn"])
                    print("Task status: ", meta["status"])
        elif self.provider == "ionq":
            if on_simulator:
                simulator = LocalSimulator()

                # Insert identity when a qubit is not targeted by any gates
                prog = self.prog.copy()
                used_qubits = {qubit for inst in prog.instructions for qubit in inst.target}
                max_qubit = max(used_qubits)
                for index in range(max_qubit):
                    if index not in used_qubits:
                        prog.i(index)

                self.task = simulator.run(prog, shots=shots)
                if verbose >= 0:
                    print("Submitted.")
            else:
                ionq_qpu = AwsDevice("arn:aws:braket:us-east-1::device/qpu/ionq/" + self.device)
                user_agent = f"SimuQ/{_version.__version__}"
                ionq_qpu.aws_session.add_braket_user_agent(user_agent)
                self.task = ionq_qpu.run(self.prog, shots=shots)
                meta = self.task.metadata()
                if verbose >= 0:
                    print("Submitted.")
                    print("Task arn: ", meta["quantumTaskArn"])
                    print("Task status: ", meta["status"])

    def results(self, task_arn=None, verbose=0):
        if task_arn is not None:
            task = AwsQuantumTask(arn=task_arn)
        else:
            if self.task is None:
                raise Exception("No submitted job in record.")
            task = self.task

        state = task.state()
        if state != "COMPLETED":
            if verbose >= 0:
                print("Job is not completed")
                print(task.metadata())
            return None

        result = task.result()

        if self.provider == "quera":
            N = len(self.layout)
            counts = dict([])
            for i in range(result.task_metadata.shots):
                if sum(result.measurements[i].pre_sequence) == N:
                    s = ""
                for j in range(N):
                    s += "g" if result.measurements[i].post_sequence[j] == 1 else "r"
                if s not in counts.keys():
                    counts[s] = 1
                else:
                    counts[s] += 1

            prob = _extract_prob(counts, N)
            ret = {
                BraketProvider.to_bin(i, N): prob[i] for i in range(1 << N) if abs(prob[i]) > 1e-6
            }
            return ret
        elif self.provider == "ionq":
            return dict(sorted(result.measurement_probabilities.items()))


def _show_register(register):
    filled_sites = [
        site.coordinate for site in register._sites if site.site_type == SiteType.FILLED
    ]
    empty_sites = [site.coordinate for site in register._sites if site.site_type == SiteType.VACANT]

    fig = plt.figure(figsize=(4, 4))
    if len(filled_sites) > 0:
        plt.plot(
            np.array(filled_sites)[:, 0],
            np.array(filled_sites)[:, 1],
            "r.",
            ms=15,
            label="filled",
        )
    if len(empty_sites) > 0:
        plt.plot(
            np.array(empty_sites)[:, 0],
            np.array(empty_sites)[:, 1],
            ".",
            color="k",
            ms=2,
            label="vacant",
        )
    plt.legend(bbox_to_anchor=(1.1, 1.05))
    return fig


def _show_global_drive(drive):
    data = {
        "detuning [rad/s]": drive.detuning.time_series,
        "amplitude [rad/s]": drive.amplitude.time_series,
        "phase [rad]": drive.phase.time_series,
    }

    fig, axes = plt.subplots(3, 1, figsize=(6, 4), sharex=True)
    for ax, data_name in zip(axes, data.keys()):
        if data_name == "phase [rad]":
            ax.step(data[data_name].times(), data[data_name].values(), ".-", where="post")
        else:
            ax.plot(data[data_name].times(), data[data_name].values(), ".-")
        ax.set_ylabel(data_name)
        ax.grid(ls=":")
    axes[-1].set_xlabel("time [s]")
    plt.tight_layout()

    return fig


def _extract_prob(counts, N):
    tot = sum(counts.values())
    freq = {k: counts[k] / tot for k in counts.keys()}
    exp_prob = []
    for i in range(1 << N):
        key = _int_to_gr(N, i)
        exp_prob.append(freq[key] if key in freq.keys() else 0)
    return exp_prob


def _int_to_gr(N, s):
    ret_list = ["g" if (s >> i) & 1 == 0 else "r" for i in range(N)]
    return "".join(ret_list)
