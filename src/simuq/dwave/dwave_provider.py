from simuq.provider import BaseProvider
from simuq.solver import generate_as
import numpy as np
import json
import re
from dwave.system import DWaveSampler, EmbeddingComposite, FixedEmbeddingComposite
from dimod import ising_to_qubo

from simuq.aais import ising
from simuq.dwave.dwave_transpiler import DwaveTranspiler


class DWaveProvider(BaseProvider):
    def __init__(self, api_key=None, from_file = None):
        # insert all log in details
        super().__init__()
        self._samples = None
        self.api_key = api_key
        if from_file is not None:
            with open(from_file, "r") as f:
                self.api_key = f.readline().strip()
        if self.api_key is None:
            raise Exception("No API key provided.")

    def compile(self,
                qs,
                anneal_schedule=[[0, 0], [20, 1]],
                chain_strength=0.5,
                verbose=0):
        h = [0 for _ in range(qs.num_sites)]
        J = {}
        for ham in qs.evos[0][0].ham:
            keys = list(ham[0].d.keys())
            vals = list(ham[0].d.values())
            if 'Z' in vals:
                if len(vals) == 1:
                    h[keys[0]] = ham[1]
                elif len(vals) == 2 and ham[1] != 0:
                    J[(keys[0], keys[1])] = ham[1]
        self.prog = h, J, anneal_schedule
        self.chain_strength = chain_strength
        return h, J

    def compare_qubo(self, q1, q2):
        numDifferent = 0
        if len(q1.keys()) != len(q2.keys()): return False
        for (d1, d2) in q1.keys():
            if abs(q1[(d1, d2)] - q2[(int(d1), int(d2))]) > 10**-3:
                numDifferent += 1
        return numDifferent


    def run(self, shots = 100, solver="Advantage_system6.4"):
        self.shots = shots
        if self.prog is None:
            raise Exception("No compiled job in record.")
        qpu = DWaveSampler(token=self.api_key, solver=solver)
        sampler = EmbeddingComposite(qpu)
        h, J, anneal_schedule = self.prog
        response = sampler.sample_ising(h, J,
                                       chain_strength=self.chain_strength,
                                       num_reads=self.shots,
                                       anneal_schedule=anneal_schedule,
                                       answer_mode="raw"
                                        )
        self.samples = list(response.samples())
        self.time_on_machine = response.info['timing']['qpu_access_time'] * 1e-6
        self.avg_qpu_time = response.info['timing']['qpu_access_time'] * 10e-6 / self.shots
        self.num_occurrences = list(response.data_vectors['num_occurrences'])
        return response

    def isingToqubo(self, h, J):
        n = len(h)
        QUBO = {}

        for i in range(n):
            s = 0
            for ii in range(n):
                if (i,ii) in J.keys():
                    s += J[(i,ii)]
                if (ii,i) in J.keys():
                    s += J[(ii,i)]

            QUBO[(i,i)] = -2 * (h[i] + s)

            for j in range(i+1, n):
                if (i,j) in J.keys() and J[i, j] != 0:
                    QUBO[(i,j)] = 4 * J[(i,j)]

        return QUBO

    def run_qubo(self):
        if self.prog is None:
            raise Exception("No compiled job in record.")
        qpu = DWaveSampler(token=self.api_key)
        sampler = EmbeddingComposite(qpu)
        h, J, anneal_schedule = self.prog
        h = {i: h[i] for i in range(len(h))}
        qubo = self.isingToqubo(h, J)
        max_interaction_qhd = np.max(np.abs(np.array(list(qubo.values()))))
        response = sampler.sample_qubo(qubo,
                                       num_reads=self.shots,
                                       anneal_schedule=anneal_schedule,
                                       chain_strength=1.1*max_interaction_qhd)
        self.samples = list(response.samples())

    def results(self):
        if self.samples == None:
            raise Exception("Job has not been run yet.")
        return self.samples
