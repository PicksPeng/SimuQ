from simuq.environment import Fermion
from simuq.qsystem import QSystem

# Fermi-Hubbard model on a 1D chain with open-boundary condition,
# based on arXiv:2010.07965 (but omitted the "spin-dependent local
# potentials" for simplicity).


L = 5  # number of lattice sites
# number of fermionic systems. x2 because of possibility of being
# spin up/down at each lattice site. 2*i (2*i+1, resp.) will
# represent a fermion at site i with spin up (down, resp.).
N = 2 * L

qs = QSystem()
fermions = [Fermion(qs) for _ in range(N)]
nops = [fermions[k].c * fermions[k].a for k in range(N)]  # number operators

J = 0.1  # hopping integral
U = 0.3  # strength of on-site interaction
hh = 0  # hopping term
ho = 0  # on-site interaction term

for i in range(L - 1):
    for s in range(2):
        k1, k2 = 2 * i + s, 2 * (i + 1) + s
        hh += -J * (fermions[k1].c * fermions[k2].a + fermions[k2].c * fermions[k1].a)
for i in range(L):
    ho += U * nops[2 * i] * nops[2 * i + 1]

h = hh + ho
t = 0.1
qs.add_evolution(h, t)
