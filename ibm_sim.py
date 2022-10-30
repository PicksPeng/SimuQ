from qiskit import IBMQ
from qiskit import transpile as IBMtranspile
from qiskit.providers.fake_provider import FakeJakarta
from qiskit.providers.aer import AerSimulator
import numpy as np

n = 3

def sim_run(circ, ideal = False) :
    if ideal :
        sim = AerSimulator()
    else :
        sim = AerSimulator.from_backend(FakeJakarta())
    circ_t = IBMtranspile(circ, backend = sim)
    res = sim.run(circ_t, shots=10000).result()
    rc = res.get_counts(0)
    return rc

def get_distr(rc) :
    cnt = np.zeros(n, dtype = 'float')
    total = 0
    for s in rc.keys() :
        total += rc[s]
        for i in range(n) :
            if s[i] == '1' :
                cnt[i] += rc[s]
    cnt /= total
    return cnt

def interp(cnts, mult = 10) :
    l = len(cnts)
    cnts = [[0 for i in range(n)]] + cnts
    L = l * mult
    d = np.zeros((n, L+1), dtype='float')
    for i in range(l) :
        for j in range(n) :
            for k in range(1, mult + 1) :
                d[j, i * mult + k] = cnts[i][j] + k * 1. / mult * (cnts[i + 1][j] - cnts[i][j])
    return d

from simuq.solver import generate_as
from systems.mis import qs, set_qs
from aais.iontrap import mach
from backends.qiskit_iontrap import transpile
import matplotlib.pyplot as plt

cnts = []
#m = 100
m = 10
for i in range(3, 4) :#m) :
    print('Now ', i)
    set_qs(m, i + 1)
    #circ=transpile(*generate_as(qs, mach, 1, 0.5))
    #cnts.append(get_distr(sim_run(circ, ideal = False)))
    circ=transpile(*generate_as(qs, mach, 100, 0.5))
    cnts.append(get_distr(sim_run(circ, ideal = True)))
data = interp(cnts, mult = int(100. / m))

fig, ax = plt.subplots(figsize = (6, 4))
shw = ax.imshow(data, interpolation = "nearest", aspect = "auto", extent = [0, 1, 0, 3 + 0.5], vmin=0, vmax=0.65)
ax.set_xlabel("time (μs)")
ax.set_ylabel("site")
ax.set_xticks(list(np.linspace(0, 1, 6)))
ax.set_yticks([1, 2, 3])
bar = fig.colorbar(shw)
fig.savefig("fig.png")

