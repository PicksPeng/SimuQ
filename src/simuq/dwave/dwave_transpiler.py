import numpy as np
from dwave.system import DWaveSampler, EmbeddingComposite
from simuq.transpiler import Transpiler


class DwaveTranspiler(Transpiler):
    def transpile(self, sol_gvars, boxes, edges, num_sites):
        h = [sol_gvars[i] for i in range(num_sites)]
        i = num_sites
        J = dict()
        for k in range(num_sites):
            for j in range(k + 1, num_sites):
                J[(j, k)] = num_sites[i]
                i += 1

        anneal_schedule = []
        for box in boxes:
            s = box[0][0][3][0]
            anneal_schedule.append([s, box[1]])

        return h, J, anneal_schedule



