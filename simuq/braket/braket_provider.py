import matplotlib.pyplot as plt
import numpy as np
from braket.ahs.atom_arrangement import SiteType
from braket.aws import AwsDevice, AwsQuantumTask
from braket.devices import LocalSimulator

from simuq import _version
from simuq.aais import rydberg1d_global, rydberg2d_global
from simuq.braket.braket_rydberg_transpiler import BraketRydbergTranspiler
from simuq.provider import BaseProvider
from simuq.solver import generate_as


class BraketProvider(BaseProvider):
    def __init__(self):
        self.backend_aais = dict()
        self.backend_aais[("ionq", "harmony")] = ["heisenberg"]
        self.backend_aais[("ionq", "aria-1")] = ["heisenberg"]
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
        no_main_body=False,
        verbose=0,
    ):
        if (provider, device) not in self.backend_aais.keys():
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

    def visualize_quera(self):
        braket_prog = self.ahs_prog

        def show_register(register):
            filled_sites = [
                site.coordinate for site in register._sites if site.site_type == SiteType.FILLED
            ]
            empty_sites = [
                site.coordinate for site in register._sites if site.site_type == SiteType.VACANT
            ]

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

        def show_global_drive(drive):
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

        fig1 = show_register(braket_prog.register)
        fig2 = show_global_drive(braket_prog.hamiltonian)
        plt.show()

    def visualize(self):
        if self.prog is None:
            raise Exception("No compiled job in record.")
        if self.provider == "quera":
            self.visualize_quera()

    def run(self, shots=1000, on_simulator=False, verbose=0):
        if self.prog is None:
            raise Exception("No compiled job in record.")

        if self.provider == "quera":
            if on_simulator:
                simulator = LocalSimulator("braket_ahs")
                self.task = simulator.run(self.ahs_prog, shots=shots)
                meta = self.task.metadata()
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

        def int_to_gr(N, s):
            ret = ""
            for i in range(N):
                ret += "g" if (s >> i) & 1 == 0 else "r"
            return ret

        def extract_prob(counts):
            tot = sum(counts.values())
            freq = dict([])
            for k in counts.keys():
                freq[k] = counts[k] / tot
            exp_prob = []
            for i in range(1 << N):
                key = int_to_gr(N, i)
                if key not in freq.keys():
                    exp_prob.append(0)
                else:
                    exp_prob.append(freq[key])
            return exp_prob

        prob = extract_prob(counts)
        ret = dict()
        for i in range(1 << N):
            if abs(prob[i]) > 1e-6:
                ret[BraketProvider.to_bin(i, N)] = prob[i]
        return ret
