from math import sqrt
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


def tfim_3to2_transform(qs, penalty):
    """Transforming 3-local transverse-field Ising model (TFIM) to 2-local TFIM
    using 2nd order perturbation theory.
    Here 3-local TFIM means
        H = \sum_ijk J_ijk Z_i Z_j Z_k
            + \sum_ij J_ij Z_i Z_j
            + \sum_i   J_i Z_i
            + \sum_i   h_i X_i
    and 2-local TFIM means
        H = \sum_ij J_ij Z_i Z_j
            + \sum_i J_i Z_i
            + \sum_i h_i X_i
    Parameters
    ----------
    qs : QSystem
        quantum system to be transformed
    penalty : float
        penalty coefficient

    Returns
    -------
    QSystem, list[Qubit]
        transformed quantum system, along with the list of new qubits
    """
    if penalty <= 0:
        raise ValueError("penalty must be positive")
    if len(qs.evos) != 1:
        raise NotImplementedError("Only a single evolution is supported for now")
    h, t = qs.evos[0]
    n = qs.num_sites # num of computational qubits
    m = sum(1 for prod, _ in h.ham if len(prod) == 3) # num of ancillary qubits
    N = n + m # total num of physical qubits

    new_qs = QSystem()
    new_sites = [Qubit(new_qs) for _ in range(N)]

    Hpen = 0 # penalty hamiltonian
    V1, V2 = 0, 0 # perturbations

    anc_ind = n
    for prod, c in h.ham:
        inds, paulistr = list(prod.keys()), "".join(prod.values())
        if paulistr == "": # ignoring terms proportional to the identity
            continue
        elif paulistr == "X":
            i = inds[0]
            V1 += c * new_sites[i].X
        elif paulistr == "Z":
            i = inds[0]
            V1 += c * new_sites[i].Z
        elif paulistr == "ZZ":
            i0, i1 = inds[0], inds[1]
            V1 += c * new_sites[i0].Z * new_sites[i1].Z
        elif paulistr == "ZZZ":
            i0, i1, i2 = inds[0], inds[1], inds[2]
            V1 += abs(c) * (5/3) * (new_sites[i0].Z * new_sites[i1].Z +
                                        new_sites[i1].Z * new_sites[i2].Z +
                                        new_sites[i2].Z * new_sites[i0].Z)
            V1 -= c * (11/3) * (new_sites[i0].Z + new_sites[i1].Z + new_sites[i2].Z)
            V2 += sqrt(32 * abs(c)) * new_sites[anc_ind].X
            Hpen += (1 - new_sites[anc_ind].Z) / 2
            if c > 0:
                Hpen += (1 + new_sites[i0].Z) * (1 - new_sites[anc_ind].Z) / 4
                Hpen += (1 + new_sites[i1].Z) * (1 - new_sites[anc_ind].Z) / 4
                Hpen += (1 + new_sites[i2].Z) * (1 - new_sites[anc_ind].Z) / 4
            else:
                Hpen += (1 - new_sites[i0].Z) * (1 - new_sites[anc_ind].Z) / 4
                Hpen += (1 - new_sites[i1].Z) * (1 - new_sites[anc_ind].Z) / 4
                Hpen += (1 - new_sites[i2].Z) * (1 - new_sites[anc_ind].Z) / 4
            anc_ind += 1
        else:
            raise Exception("Invalid Hamiltonian")

    new_h = penalty * Hpen + V1 + sqrt(penalty) * V2
    new_qs.add_evolution(new_h, t)

    return new_qs, new_sites