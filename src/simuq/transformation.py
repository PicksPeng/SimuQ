import itertools
from collections import defaultdict
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
    r"""Transform 3-local transverse-field Ising model (TFIM) to 2-local TFIM
    using 2nd order perturbation theory, such that the former effectively emerges 
    in the low-energy subspace of the latter (provided that `penalty` is 
    sufficiently large).
    Here 3-local TFIM means
        H = \sum_ijk J_ijk Z_i Z_j Z_k
            + \sum_ij J_ij Z_i Z_j
            + \sum_i   J_i Z_i
            + \sum_i   h_i X_i
    and 2-local TFIM means
        H = \sum_ij J_ij Z_i Z_j
            + \sum_i J_i Z_i
            + \sum_i h_i X_i
    New ancillary qubits will be numbered after the original computational qubits.

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


def ising_3to2_transform(qs, variant, penalty=None, peek=None):
    r"""Transform 3-local Ising model to 2-local Ising model (i.e. quadratization of 
    pseudo-Boolean functions), such that ground space is preserved (some variants 
    require that `penalty` be sufficiently large).
    Here 3-local Ising model means
        H = \sum_ijk J_ijk Z_i Z_j Z_k
            + \sum_ij J_ij Z_i Z_j
            + \sum_i   J_i Z_i
    and 2-local Ising model means
        H = \sum_ij J_ij Z_i Z_j
            + \sum_i J_i Z_i
    New ancillary qubits will be numbered after the original computational qubits.

    Parameters
    ----------
    qs : QSystem
        quantum system to be transformed
    variant : str
        specify which variant of the transformation is called, must be one of:
        "sub", "min_sel"
    penalty : float
        penalty coefficient, must be supplied for variant "sub", will be ignored 
        for variant "min_sel".
    peek : dict
        store auxilliary information about the transformation. This could be 
        useful for debug purposes.

    Returns
    -------
    QSystem, list[Qubit]
        transformed quantum system, along with the list of new qubits
    """
    if len(qs.evos) != 1:
        raise NotImplementedError("Only a single evolution is supported for now")
    if peek is not None:
        if not isinstance(peek, dict):
            raise ValueError("peek has to be a dict")

    if variant == "sub":
        if penalty is None or penalty <= 0:
            raise ValueError("penalty must be positive")
        new_qs, new_qubits = _ising_3to2_transform_sub(qs, penalty, peek)
    elif variant == "min_sel":
        new_qs, new_qubits = _ising_3to2_transform_min_sel(qs)
    else:
        raise ValueError("Unknown value of `variant`")

    return new_qs, new_qubits


def _ising_3to2_transform_sub(qs, penalty, peek=None):
    """Transform 3-local Ising model to 2-local Ising model, using the technique of substitution,
    such that ground space is preserved when `penalty` is sufficiently large.
    New ancillary qubits will be numbered after the original computational qubits.

    Procedure: Convert all Pauli Z's to number operators `n := (1-Z) / 2`; select the pair `{n_i, n_j}` 
    that appears most often among all the 3-local terms, replace it with `n_a` acting on an ancillary 
    qubit `a`, and add a penalty term to realize constraint `n_i n_j |g> = n_a |g>` for any ground state 
    `|g>`; repeat until no more 3-local terms remain; convert number operators back to Pauli Z's.

    For the exact form of the penalty term, we use the formula from
    https://docs.dwavesys.com/docs/latest/handbook_reformulating.html#example-boolean-variables.

    Parameters
    ----------
    qs : QSystem
        quantum system to be transformed
    penalty : float
        penalty coefficient
    peek : dict
        store auxilliary information about the transformation. This could be useful for debug purposes.
        peek["pair2anc"] will be a dict that maps frozenset([i,j]) to index of the ancillary qubit used
        to substitute n_i*n_j.

    Returns
    -------
    QSystem, list[Qubit]
        transformed quantum system, along with the list of new qubits
    """
    h, t = qs.evos[0]

    # ---- convert Pauli Z to number operator ---
    coeff_nnn, coeff_nn, coeff_n = _helper_func_convert_zzz_to_nnn(h)

    # --- count frequency of each pair appearing in 3-local terms ---
    pair2freq = defaultdict(int)
    # e.g. pair2freq[frozenset([i,j])] is number of 3-local terms in which qubits i and j both appear
    for triple in coeff_nnn.keys():
        for pair in itertools.combinations(triple, 2):
            pair2freq[frozenset(pair)] += 1

    # --- reduce 3-local nnn to 2-local interactions ---
    pair2anc = dict()
    # e.g. pair2anc[frozenset([i,j])] records which ancillary qubit is used to substitute n_i*n_j
    nq = qs.num_sites # num of qubits
    while len(coeff_nnn) > 0:
        # find the pair {i,j} with largest frequency
        ij, _ = max(pair2freq.items(), key=lambda x:x[1])
        i, j = ij
        # introduce an ancillary qubit, whose number operator n_a will substitute n_i*n_j
        pair2anc[ij] = a = nq
        nq += 1
        # update coeff and pair2freq
        all_triples = list(coeff_nnn.keys())
        for ijk in all_triples:
            if ij.issubset(ijk):
                k = list(ijk - ij)[0]
                coeff_nn[frozenset([a, k])] += coeff_nnn[ijk]
                del coeff_nnn[ijk]
                ik, jk = frozenset([i, k]), frozenset([j, k])
                pair2freq[ik] -= 1
                if pair2freq[ik] == 0:
                    del pair2freq[ik]
                pair2freq[jk] -= 1
                if pair2freq[jk] == 0:
                    del pair2freq[jk]
        del pair2freq[ij]
    assert len(pair2freq) == 0

    if peek is not None:
        peek["pair2anc"] = pair2anc

    # --- add penalty terms ---
    for ij, a in pair2anc.items():
        i, j = ij
        ia, ja = frozenset([i, a]), frozenset([j, a])
        # To make sure that n_i * n_j = n_a on the ground space, add a penalty term:
        # penalty * (n_i * n_j - 2 (n_i + n_j) n_a + 3 n_a)
        coeff_nn[ij] += penalty
        coeff_nn[ia] += -2 * penalty
        coeff_nn[ja] += -2 * penalty
        coeff_n[a] += 3 * penalty

    # ---- convert number operator to Pauli Z ---
    new_qs = QSystem()
    new_qubits = [Qubit(new_qs) for _ in range(nq)]
    new_h = _helper_func_convert_nn_to_zz(new_qubits, coeff_nn, coeff_n)
    new_qs.add_evolution(new_h, t)

    return new_qs, new_qubits


def _ising_3to2_transform_min_sel(qs):
    """Transform 3-local Ising model to 2-local Ising model, using the technique of minimum selection,
    such that ground space is preserved. New ancillary qubits will be numbered after the original 
    computational qubits.

    Procedure: Convert all Pauli Z's to number operators `n := (1-Z) / 2`; replace each 3-local term 
    `n_i * n_j * n_k` by some degree-2 polynomial `p(n_i, n_j, n_k, n_a)` where `a` is an ancillary 
    qubit; convert number operators back to Pauli Z's.

    For the exact form of the degree-2 polynomial, we use the formula from
    https://docs.dwavesys.com/docs/latest/handbook_reformulating.html#polynomial-reduction-by-minimum-selection.

    Parameters
    ----------
    qs : QSystem
        quantum system to be transformed

    Returns
    -------
    QSystem, list[Qubit]
        transformed quantum system, along with the list of new qubits
    """
    h, t = qs.evos[0]

    # ---- convert Pauli Z to number operator ---
    coeff_nnn, coeff_nn, coeff_n = _helper_func_convert_zzz_to_nnn(h)

    # --- reduce 3-local nnn to 2-local interactions ---
    nq = qs.num_sites # num of qubits
    for ijk, c in coeff_nnn.items():
        i, j, k = ijk
        a = nq
        nq += 1
        if c < 0:
            coeff_n[a] += -2 * c
            coeff_nn[frozenset([i, a])] += c
            coeff_nn[frozenset([j, a])] += c
            coeff_nn[frozenset([k, a])] += c
        else:
            coeff_n[i] += -c
            coeff_n[j] += -c
            coeff_n[k] += -c
            coeff_n[a] += -c
            coeff_nn[frozenset([i, j])] += c
            coeff_nn[frozenset([i, k])] += c
            coeff_nn[frozenset([j, k])] += c
            coeff_nn[frozenset([i, a])] += c
            coeff_nn[frozenset([j, a])] += c
            coeff_nn[frozenset([k, a])] += c

    # ---- convert number operator to Pauli Z ---
    new_qs = QSystem()
    new_qubits = [Qubit(new_qs) for _ in range(nq)]
    new_h = _helper_func_convert_nn_to_zz(new_qubits, coeff_nn, coeff_n)
    new_qs.add_evolution(new_h, t)

    return new_qs, new_qubits


def _helper_func_convert_zzz_to_nnn(h):
    """Helper function used in `ising_3to2_transform()` to reexpress a 3-local 
    Ising Hamiltonian in terms of number operators.
    Caution: terms proportional to the identity will be dropped.

    Parameters
    ----------
    h : TIHamiltonian
        3-local Ising Hamiltonian

    Returns
    -------
    coeff_nnn : defaultdict(float)
        coefficients of 3-local terms,
        e.g. `coeff_nn[frozenset([i, j, k])]` is coefficient in front of `n_i * n_j * n_k`.
    coeff_nn : defaultdict(float)
        coefficients of 2-local terms,
        e.g. `coeff_nn[frozenset([i, j])]` is coefficient in front of `n_i * n_j`.
    coeff_n : defaultdict(float)
        coefficients of 1-local terms,
        e.g. `coeff_n[i]` is coefficient in front of `n_i`.
    """
    coeff_nnn, coeff_nn, coeff_n = defaultdict(float), defaultdict(float), defaultdict(float)

    for prod, c in h.ham:
        inds, paulistr = list(prod.keys()), "".join(prod.values())
        if paulistr == "": # ignoring terms proportional to the identity
            continue
        elif paulistr == "Z": # replacing Z -> 1 - 2n
            coeff_n[inds[0]] += -2 * c
        elif paulistr == "ZZ": # replacing Z_i*Z_j -> 1 - 2(n_i+n_j) + 4n_i*n_j
            for i in inds:
                coeff_n[i] += -2 * c
            coeff_nn[frozenset(inds)] += 4 * c
        elif paulistr == "ZZZ": # replacing Z_i*Z_j*Z_k -> 1 - 2(n_i+n_j+n_k) + 4(n_i*n_j+n_i*n_k+n_j*n_k) - 8n_i*n_j*n_k
            for i in inds:
                coeff_n[i] += -2 * c
            for pair in itertools.combinations(inds, 2):
                coeff_nn[frozenset(pair)] += 4 * c
            coeff_nnn[frozenset(inds)] += -8 * c
        else:
            raise ValueError("Invalid Hamiltonian")

    return coeff_nnn, coeff_nn, coeff_n


def _helper_func_convert_nn_to_zz(qubits, coeff_nn, coeff_n):
    """Helper function used in `ising_3to2_transform()` to construct a 2-local 
    Ising Hamiltonian based on coefficients of number operators.

    Parameters
    ----------
    qubits : list[Qubit]
        list of qubits on which the Hamiltonian acts on.
    coeff_nn : defaultdict(float)
        coefficients of 2-local terms,
        e.g. `coeff_nn[frozenset([i, j])]` is coefficient in front of `n_i * n_j`.
    coeff_n : defaultdict(float)
        coefficients of 1-local terms,
        e.g. `coeff_n[i]` is coefficient in front of `n_i`.

    Returns
    -------
    h : TIHamiltonian
        3-local Ising Hamiltonian
    """
    h = 0
    for i, c in coeff_n.items(): # replacing n_i -> (1 - Z_i) / 2
        h += (1/2) * c * (1 - qubits[i].Z)
    for ij, c in coeff_nn.items(): # replacing n_i * n_j -> (1 - Z_i) * (1 - Z_j) / 4
        i, j = ij
        h += (1/4) * c * (1 - qubits[i].Z) * (1 - qubits[j].Z)
    return h