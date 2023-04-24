from simuq.qsystem import QSystem
from simuq.environment import qubit, fermion


# The Jordan-Wigner transformation
# Converts all fermions in the system into spins
def JW_transform(qs) :

    n = qs.num_sites

    new_qs = QSystem()
    new_sites = []
    for i in range(n) :
        if isinstance(qs.sites[i], fermion) :
            new_sites.append(qubit(new_qs))
        elif isinstance(qs.sites[i], qubit) :
            new_sites.append(qubit(new_qs))
        elif isinstance(qs.sites[i], boson) :
            new_sites.append(boson(new_qs))

    for evo_ind in range(len(qs.evos)) :
        (h, t) = qs.evos[evo_ind]
        new_h = 0
        for (prod, c) in h.ham :
            new_prod = c
            prev_Z = 1
            for i in range(n) :
                if isinstance(qs.sites[i], fermion) :
                    for op in prod[i] :
                        if op == 'a' :
                            new_prod *= prev_Z * (new_sites[i].X + 1j * new_sites[i].Y) / 2
                        elif op == 'c' :
                            new_prod *= prev_Z * (new_sites[i].X - 1j * new_sites[i].Y) / 2
                    prev_Z *= new_sites[i].Z
                elif isinstance(qs.sites[i], qubit) :
                    for op in prod[i] :
                        if op == 'X' :
                            new_prod *= new_sites[i].X
                        elif op == 'Y' :
                            new_prod *= new_sites[i].Y
                        elif op == 'Z' :
                            new_prod *= new_sites[i].Z
                elif isinstance(qs.sites[i], boson) :
                    for op in prod[i] :
                        if op == 'a' :
                            new_prod *= new_sites[i].a
                        elif op == 'c' :
                            new_prod *= new_sites[i].c
            new_h += new_prod
        new_qs.add_evolution(new_h, t)

    return new_qs, new_sites
