import numpy as np

from simuq.environment import Boson, Fermion, Qubit
from simuq.qsystem import QSystem


# The Jordan-Wigner transformation for fermionic modes
# Converts all fermions in the system into spins
def jw_transform(qs):
    n = qs.num_sites

    new_qs = QSystem()
    new_sites = []
    for i in range(n):
        if isinstance(qs.sites[i], Fermion):
            new_sites.append(Qubit(new_qs))
        elif isinstance(qs.sites[i], Qubit):
            new_sites.append(Qubit(new_qs))
        elif isinstance(qs.sites[i], Boson):
            new_sites.append(Boson(new_qs))

    for evo_ind in range(len(qs.evos)):
        (h, t) = qs.evos[evo_ind]
        new_h = 0
        for prod, c in h.ham:
            new_prod = c
            prev_Z = 1
            for i in range(n):
                if isinstance(qs.sites[i], Fermion):
                    for op in prod[i]:
                        if op == "a":
                            new_prod *= prev_Z * (new_sites[i].X + 1j * new_sites[i].Y) / 2
                        elif op == "c":
                            new_prod *= prev_Z * (new_sites[i].X - 1j * new_sites[i].Y) / 2
                    prev_Z *= new_sites[i].Z
                elif isinstance(qs.sites[i], Qubit):
                    for op in prod[i]:
                        if op == "X":
                            new_prod *= new_sites[i].X
                        elif op == "Y":
                            new_prod *= new_sites[i].Y
                        elif op == "Z":
                            new_prod *= new_sites[i].Z
                elif isinstance(qs.sites[i], Boson):
                    for op in prod[i]:
                        if op == "a":
                            new_prod *= new_sites[i].a
                        elif op == "c":
                            new_prod *= new_sites[i].c
            new_h += new_prod
        new_qs.add_evolution(new_h, t)

    return new_qs, new_sites


# The one-hot transformation for bosonic modes
# Converts all bosons in the system into spins
# Encode |n> into |11..11011..11> where the n-th bit is 0
def oh_transform(qs, truncation_levels=3):
    n = qs.num_sites
    m = truncation_levels

    new_qs = QSystem()
    new_sites = []
    label = []
    for i in range(n):
        if isinstance(qs.sites[i], Fermion):
            label[i] = len(new_sites)
            new_sites.append(Fermion(new_qs))
        elif isinstance(qs.sites[i], Qubit):
            label[i] = len(new_sites)
            new_sites.append(Qubit(new_qs))
        elif isinstance(qs.sites[i], Boson):
            label[i] = len(new_sites)
            new_sites += [Qubit(new_qs) for _ in range(m)]

    for evo_ind in range(len(qs.evos)):
        (h, t) = qs.evos[evo_ind]
        new_h = 0
        for prod, c in h.ham:
            new_prod = c
            for i in range(n):
                if isinstance(qs.sites[i], Fermion):
                    for op in prod[i]:
                        if op == "a":
                            new_prod *= new_sites[label[i]].a
                        elif op == "c":
                            new_prod *= new_sites[label[i]].c
                elif isinstance(qs.sites[i], Qubit):
                    for op in prod[i]:
                        if op == "X":
                            new_prod *= new_sites[label[i]].X
                        elif op == "Y":
                            new_prod *= new_sites[label[i]].Y
                        elif op == "Z":
                            new_prod *= new_sites[label[i]].Z
                elif isinstance(qs.sites[i], Boson):
                    for op in prod[i]:
                        sigminus = [
                            (new_sites[label[i] + j].X - 1j * new_sites[label[i] + j].Y) / 2
                            for j in range(m)
                        ]
                        sigplus = [
                            (new_sites[label[i] + j].X + 1j * new_sites[label[i] + j].Y) / 2
                            for j in range(m)
                        ]
                        if op == "a":
                            tr_a = 0
                            for j in range(m - 1):
                                tr_a += np.sqrt(j + 1) * sigplus[j] * sigminus[j + 1]
                            new_prod *= tr_a
                        elif op == "c":
                            tr_c = 0
                            for j in range(m - 1):
                                tr_c += np.sqrt(j + 1) * sigminus[j] * sigplus[j + 1]
                            new_prod *= tr_c
            new_h += new_prod
        new_qs.add_evolution(new_h, t)

    return new_qs, new_sites
