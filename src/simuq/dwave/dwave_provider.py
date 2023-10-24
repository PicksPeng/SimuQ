from simuq.provider import BaseProvider
from simuq.solver import generate_as
import numpy as np
from dwave.system import DWaveSampler, EmbeddingComposite

from simuq.aais import ising
from simuq.dwave.dwave_transpiler import DwaveTranspiler


class DWaveProvider(BaseProvider):
    def __init__(self, api_key):
        # insert all log in details
        super().__init__()
        self._samples = None
        self.api_key = api_key

    def compile(self,
                qs,
                tol=0.01,
                trotter_num=6,
                verbose=0):
        mach = ising.generate_qmachine(qs.num_sites)
        layout, sol_gvars, boxes, edges = generate_as(
            qs,
            mach,
            trotter_num,
            solver="least_squares",
            solver_args={"tol": tol},
            override_layout=None,
            verbose=verbose,
        )
        self.prog = DwaveTranspiler.transpile(sol_gvars, boxes, edges, qs.num_sites)

    def run(self):
        if self.prog is None:
            raise Exception("No compiled job in record.")
        qpu = DWaveSampler()
        sampler = EmbeddingComposite(qpu)
        h, J, anneal_schedule = self.prog
        self._samples = sampler.sample_ising(h, J, anneal_schedule, return_embedding=True, token=self.api_key).samples()

    def results(self):
        if self._samples == None:
            raise Exception("Job has not been run yet.")
        return self._samples