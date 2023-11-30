from simuq.provider import BaseProvider
from simuq.solver import generate_as
import numpy as np
import json
import re
from dwave.system import DWaveSampler, EmbeddingComposite
from dimod import ising_to_qubo

from simuq.aais import ising
from simuq.dwave.dwave_transpiler import DwaveTranspiler


class DWaveProvider(BaseProvider):
    def __init__(self, api_key, numruns=100):
        # insert all log in details
        super().__init__()
        self._samples = None
        self.api_key = api_key
        self.numruns = numruns

    def compile(self,
                qs,
                tol=0.01,
                trotter_num=6,
                verbose=0):
        # mach = ising.generate_qmachine(qs.num_sites)
        # layout, sol_gvars, boxes, edges = generate_as(
        #     qs,
        #     mach,
        #     trotter_num,
        #     solver="least_squares",
        #     solver_args={"tol": tol},
        #     override_layout=None,
        #     verbose=verbose,
        # )
        h = [0 for _ in range(qs.num_sites)]
        J = {(i, j): 0 for i in range(qs.num_sites) for j in range(i + 1, qs.num_sites)}
        for ham in qs.evos[0][0].ham:
            keys = list(ham[0].d.keys())
            vals = list(ham[0].d.values())
            if 'Z' in vals:
                if len(vals) == 1:
                    h[keys[0]] = ham[1]
                elif len(vals) == 2:
                    J[(keys[0], keys[1])] = ham[1]

        self.prog = h, J, None

    def compareToQubo(self,jsonQM, isingQM):
        errors = []
        # with open("qubomethodQ.txt", 'r', encoding='utf-8') as f:
        #     jsonQM = json.loads(json.loads(f.read()))
        for key in isingQM.keys():
            if abs(jsonQM[key] - isingQM[key]) > 10 ** -10:
                errors.append(key)
        print(errors)
        return errors

    def returnOldQubo(self):
        with open("qubomethodQ.txt", 'r', encoding='utf-8') as f:
            jsonQM = json.loads(json.loads(f.read()))
        newQubo = {}
        for key in jsonQM.keys():
            a, b = self.extract_and_convert_numbers(key)
            newQubo[(a, b)] = jsonQM[key]
        return newQubo

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

    def extract_and_convert_numbers(self, input_string):
        # Regex pattern to match two positive integers separated by a comma
        pattern = r"(\d+),(\d+)"

        # Using regex to find matches
        match = re.match(pattern, input_string)

        # Extracting the numbers a and b if a match is found
        if match:
            a, b = match.groups()
            return int(a), int(b)
        else:
            return None, None

    def run(self):
        if self.prog is None:
            raise Exception("No compiled job in record.")
        qpu = DWaveSampler(token=self.api_key)
        sampler = EmbeddingComposite(qpu)
        h, J, anneal_schedule = self.prog
        h = {i: h[i] for i in range(len(h))}
        qubo = self.isingToqubo(h, J)
        dwaveQubo = ising_to_qubo(h, J)
        quboCompare = {f'{a},{b}': qubo[(a, b)] for (a, b) in qubo.keys()}
        # quboForDwave = {(str(a), str(b)): qubo[(a, b)] for key in qubo.keys()}
        # errorCount = self.compareToQubo(self.returnOldQubo(), qubo)
        max_interaction_qhd = np.max(np.abs(np.array(list(qubo.values()))))

        response = sampler.sample_qubo(qubo, chain_strength=max_interaction_qhd,
                                       num_reads=self.numruns)
        self.samples = list(response.samples())

    def results(self):
        if self._samples == None:
            raise Exception("Job has not been run yet.")
        return self._samples