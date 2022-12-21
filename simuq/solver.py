"""
The solver of SimuQ.

The solver contains three steps:
(1) find an alignment of target sites and machine sites.
(2) synthesize the target Hamiltonian with instructions
without considering mutually exclusiveness.
(3) trotterize the resulting boxes of instructions.

The detailed descriptions are left in the implementations.

TODO: optimize synthesize procedure for cases without
global variables.
"""

import logging
import math

import networkx as nx
import numpy as np
import scipy.optimize as opt

from simuq.hamiltonian import TIHamiltonian

logging.basicConfig(level=logging.CRITICAL)

logger = logging.getLogger()
# Functions that locate a local variable.
def locate_switch_evo(mach, evo_index):
    return mach.num_gvars + evo_index * (mach.num_inss + mach.num_lvars)


def switch_term(mach, evo_index, ins_index, mc):
    return lambda x: x[locate_switch_evo(mach, evo_index) + ins_index] * mc.exp_eval(
        x[: mach.num_gvars],
        x[
            locate_switch_evo(mach, evo_index)
            + mach.num_inss : locate_switch_evo(mach, evo_index + 1)
        ],
    )


def switch_fun(mach, evo_index, ins_index):
    return lambda x: x[locate_switch_evo(mach, evo_index) + ins_index]


def locate_nonswitch_evo(mach, evo_index):
    return mach.num_gvars + evo_index * mach.num_lvars


def non_switch_term(mach, evo_index, ins_index, mc):
    return lambda x: mc.exp_eval(
        x[: mach.num_gvars],
        x[
            locate_nonswitch_evo(mach, evo_index) : locate_nonswitch_evo(
                mach, evo_index + 1
            )
        ],
    )


# Equation builder and solver when an alignment is supplied.
def solve_aligned(ali, qs, mach, tol=1e-3):
    logger.info(ali)

    """
    Build equations. The equation variables are the global variables
    and copies of local variables (# evolution of copies).

    Each product term of the target Hamiltonian is used to build an
    equation, as the constant side. The other side contains all the
    possible contributor of this product in the machine's instructions.

    Since a instruction is either selected or not, we represent it
    by a switch variable for each instruction. Then it turns to a mixed
    equation solving problem ({0, 1} variables and continuous varibales).

    Currently, this is solved by first treating the switch variables as
    continuous variables, truncating them and solving the equations 
    again with switch variables set to 0 or 1.
    """
    nvar = mach.num_gvars + len(qs.evos) * (mach.num_inss + mach.num_lvars)

    def match_prod(tprod, mprod):
        return tprod == mprod
        """
        untouched = [True for j in range(mach.num_sites)]
        for k in range(qs.num_sites) :
            if tprod[k] != mprod[ali[k]] :
                return False
            untouched[ali[k]] = False
        for k in range(mach.num_sites) :
            if untouched[k] and mprod[k] != '' :
                return False
        return True
        """

    def to_mprod(tprod):
        ret = ["" for j in range(mach.num_sites)]
        for k in range(qs.num_sites):
            ret[ali[k]] = tprod[k]
        return ret

    def is_id(tprod):
        ret = True
        for i in range(mach.num_sites):
            if not tprod[i] == "":
                ret = False
                break
        return ret

    def build_eqs(ins_locator, switch_locator=None):
        eqs = []
        for evo_index in range(len(qs.evos)):
            (h, t) = qs.evos[evo_index]
            mark = [[0 for ins in line.inss] for line in mach.lines]
            targ_terms = [(to_mprod(ham), c) for (ham, c) in h.ham]
            ind = 0
            while ind < len(targ_terms):
                tprod, tc = targ_terms[ind]
                if is_id(tprod):
                    ind += 1
                    continue
                eq = (lambda c: lambda x: -c)(tc)
                for i in range(len(mach.lines)):
                    line = mach.lines[i]
                    for j in range(len(line.inss)):
                        ins = line.inss[j]
                        for (mprod, mc) in ins.h.ham:
                            if match_prod(tprod, mprod):
                                eq = (lambda eq_, f_: lambda x: eq_(x) + f_(x))(
                                    eq, ins_locator(mach, evo_index, ins.index, mc)
                                )
                                mark[i][j] = 1
                                # Check if other terms of machine instruction exists in targer Ham terms.
                                for (mprod_prime, _) in ins.h.ham:
                                    exists_in_targ_terms = False
                                    for (tprod_prime, _) in targ_terms:
                                        if match_prod(tprod_prime, mprod_prime):
                                            exists_in_targ_terms = True
                                            break
                                    if not exists_in_targ_terms:
                                        targ_terms.append((mprod_prime, 0))
                                break
                eqs.append((lambda eq_: lambda x: t * eq_(x))(eq))
                ind += 1
            if switch_locator != None:
                for i in range(len(mach.lines)):
                    line = mach.lines[i]
                    for j in range(len(line.inss)):
                        ins = line.inss[j]
                        if mark[i][j] == 0 and not ins.is_sys_ham:
                            eqs.append(switch_locator(mach, evo_index, ins.index))
            if mach.with_sys_ham and mark[-1][0] == 0:
                line = mach.lines[-1]
                ins = line.inss[0]
                for (mprod, mc) in ins.h.ham:
                    eqs.append(
                        (lambda eq_: lambda x: t * eq_(x))(
                            ins_locator(mach, evo_index, ins.index, mc)
                        )
                    )
        return eqs

    eqs = build_eqs(switch_term, switch_fun)
    
    logger.info('Number of vars and equations', nvar, len(eqs))
    
    f = (
        lambda eqs_: lambda x: [(lambda i_: eqs_[i_](x))(i) for i in range(len(eqs_))]
    )(eqs)
    lbs = [mach.gvars[j].lower_bound for j in range(mach.num_gvars)]
    ubs = [mach.gvars[j].upper_bound for j in range(mach.num_gvars)]
    init = [mach.gvars[j].init_value for j in range(mach.num_gvars)]
    for i in range(len(qs.evos)):
        if mach.with_sys_ham:
            lbs += (
                [0 for j in range(mach.num_inss - 1)]
                + [1]
                + [mach.lvars[j].lower_bound for j in range(mach.num_lvars)]
            )
            ubs += (
                [1 for j in range(mach.num_inss - 1)]
                + [1 + tol]
                + [mach.lvars[j].upper_bound for j in range(mach.num_lvars)]
            )
            init += (
                [0.5 for j in range(mach.num_inss - 1)]
                + [1 + tol / 2]
                + [mach.lvars[j].init_value for j in range(mach.num_lvars)]
            )
        else:
            lbs += [0 for j in range(mach.num_inss)] + [
                mach.lvars[j].lower_bound for j in range(mach.num_lvars)
            ]
            ubs += [1 for j in range(mach.num_inss)] + [
                mach.lvars[j].upper_bound for j in range(mach.num_lvars)
            ]
            init += [0.5 for j in range(mach.num_inss)] + [
                mach.lvars[j].init_value for j in range(mach.num_lvars)
            ]

    logger.info("Init values: ", init)
    # logger.info("Upper bound: ", ubs)
    # logger.info("Lower bound: ", lbs)
    sol_detail = opt.least_squares(f, init, bounds=(lbs, ubs))
    sol = sol_detail.x

    logger.info(np.linalg.norm(f(sol)))
    # logger.info(sol)
    if np.linalg.norm(f(sol)) > tol:
        return False

    # Solve it again, with initial value set as the previous solution.

    lbs = [mach.gvars[j].lower_bound for j in range(mach.num_gvars)]
    ubs = [mach.gvars[j].upper_bound for j in range(mach.num_gvars)]
    for i in range(len(qs.evos)):
        lbs += [mach.lvars[j].lower_bound for j in range(mach.num_lvars)]
        ubs += [mach.lvars[j].upper_bound for j in range(mach.num_lvars)]

    nvar = mach.num_gvars + len(qs.evos) * mach.num_lvars
    initvars = sol[: mach.num_gvars].tolist()
    switch = [
        [[0 for ins in line.inss] for line in mach.lines]
        for evo_index in range(len(qs.evos))
    ]
    for evo_index in range(len(qs.evos)):
        for i in range(len(mach.lines)):
            line = mach.lines[i]
            for j in range(len(line.inss)):
                ins = line.inss[j]
                if abs(sol[locate_switch_evo(mach, evo_index) + ins.index]) < 1e-3:
                    switch[evo_index][i][j] = 0
                else:
                    switch[evo_index][i][j] = 1
        initvars += sol[
            locate_switch_evo(mach, evo_index)
            + mach.num_inss : locate_switch_evo(mach, evo_index + 1)
        ].tolist()

    eqs = build_eqs(non_switch_term)
    f = (
        lambda eqs_: lambda x: [(lambda i_: eqs_[i_](x))(i) for i in range(len(eqs_))]
    )(eqs)
    var_lb = -np.inf
    var_ub = np.inf
    lbs = [var_lb for i in range(mach.num_gvars)]
    ubs = [var_ub for i in range(mach.num_gvars)]
    for i in range(len(qs.evos)):
        lbs += [var_lb for j in range(mach.num_lvars)]
        ubs += [var_ub for j in range(mach.num_lvars)]

    global gsol
    global gswitch
    global alignment
    alignment = ali
    if len(initvars) == 0:
        gsol = np.array([])
        gswitch = switch
        if np.linalg.norm(f([])) > tol:
            return False
        return True
    sol_detail = opt.least_squares(f, initvars, bounds=(lbs, ubs), )
    sol = sol_detail.x

    gsol = sol
    gswitch = switch

    logger.info(np.linalg.norm(f(sol)))
    if np.linalg.norm(f(sol)) > tol:
        return False

    logger.info(switch)
    logger.info(sol)
    return True

def solve_aligned_wrapper(ali, qs, mach, tol) :
    if mach.num_gvars != 0 or mach.with_sys_ham :
        return solve_aligned(ali, qs, mach, tol)
    else :
        switch = []
        sol = np.zeros(mach.num_lvars * len(qs.evos))
        evos = qs.evos
        global gswitch
        global gsol
        for evo_index in range(len(evos)) :
            qs.evos = [evos[evo_index]]
            if not solve_aligned(ali, qs, mach, tol) :
                return False
            switch.append(gswitch[0])
            sol[
                locate_nonswitch_evo(mach, evo_index) :
                locate_nonswitch_evo(mach, evo_index + 1)
            ] = gsol
        qs.evos = evos
        gswitch = switch
        gsol = sol
        return True


"""
The search of a valid alignemnt.

Optimization used here: if the current alignment dooms a
non-zero product term in the target Hamiltonian, then break.
"""


def align(i, ali, qs, mach, tol):
    def is_id(tprod):
        ret = True
        for i in range(qs.num_sites):
            if not tprod[i] == "":
                ret = False
                break
        return ret

    if i == qs.num_sites:
        if solve_aligned_wrapper(ali, qs, mach, tol):
            return True
        return False
    for x in range(mach.num_sites):
        available = True
        ali[i] = x
        for j in range(i):
            if ali[j] == x:
                available = False
                break
        if available == False:
            continue
        for (h, t) in qs.evos:
            for (tprod, tc) in h.ham:
                if is_id(tprod):
                    continue
                found = False
                for line in mach.lines:
                    for ins in line.inss:
                        for (mprod, mc) in ins.h.ham:
                            prod_partial_match = True
                            for k in range(i + 1):
                                if tprod[k] != mprod[ali[k]]:
                                    prod_partial_match = False
                            if prod_partial_match:
                                found = True
                                break
                        if found:
                            break
                    if found:
                        break
                if not found:
                    available = False
                    break
            if not available:
                break
        if available:
            if align(i + 1, ali, qs, mach, tol):
                return True
    return False


def find_sol(qs, mach, ali=[], tol=1e-3):
    if ali == []:
        return align(0, [0 for i in range(qs.num_sites)], qs, mach, tol)
    else:
        return solve_aligned_wrapper(ali, qs, mach, tol)


# The generation of abstract schedule
# Trotterize the solution provided by the second step.
def generate_as(qs, mach, trotter_step=4, solver_tol=1e-1):
    mach.instantiate_sys_ham()
    mach.extend_instruction_sites()
    if find_sol(qs, mach, tol=solver_tol):
        sol = gsol
        switch = gswitch

        sol_gvars = sol[: mach.num_gvars]
        boxes = []
        edges = []
        ending_boxes = []

        # Deal with machines with system Hamiltonian first.
        if mach.with_sys_ham:
            for evo_index in range(len(qs.evos)):
                (h, t) = qs.evos[evo_index]
                sol_lvars = sol[
                    locate_nonswitch_evo(mach, evo_index) : locate_nonswitch_evo(
                        mach, evo_index + 1
                    )
                ].tolist()
                box = ([], t)
                for i in range(len(mach.lines)):
                    line = mach.lines[i]
                    for j in range(len(line.inss)):
                        ins = line.inss[j]
                        if switch[evo_index][i][j] == 1 and not ins.is_sys_ham:
                            ins_lvars = []
                            for l in range(len(ins.vars_index)):
                                ins_lvars.append(sol_lvars[ins.vars_index[l]])
                            box[0].append(
                                (
                                    (i, j),
                                    ins,
                                    ins.exp_eval(sol_gvars, sol_lvars),
                                    ins_lvars,
                                )
                            )
                if box != []:
                    boxes.append(box)

            logger.info(alignment)
            logger.info(sol_gvars)
            for box in boxes:
                logger.info(
                    (
                        [((i, j), ins_lvars) for ((i, j), ins, h, ins_lvars) in box[0]],
                        box[1],
                    )
                )
            edges = [(i, i + 1) for i in range(len(qs.evos) - 1)]
            for idx in range(len(qs.evos) - 1):
                logger.info((idx, idx + 1))

            return (alignment, sol_gvars, boxes, edges)

        # Then deal with other cases
        for evo_index in range(len(qs.evos)):
            (h, t) = qs.evos[evo_index]
            next_ending_boxes = []
            coloring = [i for i in range(mach.num_sites)]
            sol_lvars = sol[
                locate_nonswitch_evo(mach, evo_index) : locate_nonswitch_evo(
                    mach, evo_index + 1
                )
            ].tolist()
            # Detach by touched sites
            """
            Note that we assume that a derived instruction only 
            affects those sites it touches (it commutes with any
            instruction not touching these sites). We can separate 
            the instructions in the box by the sites they touch.

            Then we use a connectivity coloring to separate them first.
            """
            first_touch = [[-1 for i in range(len(line.inss))] for line in mach.lines]
            for i in range(len(mach.lines)):
                line = mach.lines[i]
                for j in range(len(line.inss)):
                    ins = line.inss[j]
                    if switch[evo_index][i][j] == 1:
                        t_sites = (ins.h.exp_eval(sol_gvars, sol_lvars)).touched_sites()
                        for k in range(mach.num_sites):
                            if t_sites[k] == 1:
                                if first_touch[i][j] == -1:
                                    first_touch[i][j] = k
                                else:
                                    if coloring[k] != coloring[first_touch[i][j]]:
                                        c = coloring[k]
                                        for l in range(mach.num_sites):
                                            if coloring[l] == c:
                                                coloring[l] = coloring[
                                                    first_touch[i][j]
                                                ]
            color_part = []
            for k in range(mach.num_sites):
                ins_set = []
                for i in range(len(mach.lines)):
                    line = mach.lines[i]
                    for j in range(len(line.inss)):
                        ins = line.inss[j]
                        if (
                            switch[evo_index][i][j] == 1
                            and coloring[first_touch[i][j]] == k
                        ):
                            ins_lvars = []
                            for l in range(len(ins.vars_index)):
                                ins_lvars.append(sol_lvars[ins.vars_index[l]])
                            ins_set.append(
                                (
                                    (i, j),
                                    ins,
                                    ins.exp_eval(sol_gvars, sol_lvars),
                                    ins_lvars,
                                )
                            )
                if ins_set != []:
                    color_part.append(ins_set)
            # Detach if commute with all others
            """
            If an instruction's Hamiltonian commutes with every
            other ones in the box, then we can set itself as a
            stand-alone box, detaching from the original box.
            """
            logger.info(color_part)
            for i in range(len(color_part)):
                ins_set = color_part[i]
                j = 0
                while j < len(ins_set):
                    all_commute = True
                    for k in range(len(ins_set)):
                        if j == k:
                            continue
                        # logger.info(j, k)
                        # logger.info(ins_set[j][2].ham, ins_set[k][2].ham)
                        if not TIHamiltonian.commutativity_test(
                            ins_set[j][2],
                            ins_set[k][2],
                            ins_set[j][1].prop == "derived"
                            or ins_set[k][1].prop == "derived",
                        ):
                            all_commute = False
                            break
                    if all_commute:
                        color_part.append([ins_set[j]])
                        del ins_set[j]
                        j -= 1
                    j += 1
            # Trotterization
            """
            The mutually exclusiveness can be modeled as a graph
            where vertices are instructions.

            An edge is added either both vertices share the same
            signal line, or they have shared sites and one of them
            is an derivied instruction.

            We use a coloring algorithm to minimize the number of
            colors, such that no edge has two vertices in the same
            color. With these coloring we can further divide a box
            into multiple parallel boxes, repeated multiply times
            (due to trotterization).
            """
            for ins_set in color_part:
                local_ending_boxes = ending_boxes
                n = len(ins_set)
                G = nx.Graph()
                G.add_nodes_from(range(n))
                for i in range(n):
                    for j in range(n):
                        if i == j:
                            continue
                        if (ins_set[i][0][0] == ins_set[j][0][0]) or (
                            (
                                ins_set[i][1].prop == "derived"
                                or ins_set[j][1].prop == "derived"
                            )
                            and (
                                not TIHamiltonian.commutativity_test(
                                    ins_set[i][2], ins_set[j][2], True
                                )
                            )
                        ):
                            G.add_edge(i, j)
                col = nx.coloring.greedy_color(G, strategy="largest_first")
                nodes_of_color = [[] for i in set(col.values())]
                for i in range(n):
                    nodes_of_color[col[i]].append(i)
                # Check if the partitions are pair-wise commute
                sumh = []
                for i in range(len(nodes_of_color)):
                    h = 0
                    for j in range(len(nodes_of_color[i])):
                        h = h + ins_set[nodes_of_color[i][j]][2]
                    sumh.append(h)
                all_commute = True
                for i in range(len(nodes_of_color)):
                    for j in range(i):
                        if TIHamiltonian.commutativity_test(sumh[i], sumh[j]) == False:
                            all_commute = False
                            break
                    if not all_commute:
                        break
                if all_commute:
                    steps = 1
                else:
                    steps = trotter_step
                for i in range(steps):
                    next_local_ending_boxes = []
                    for k in range(len(nodes_of_color)):
                        list_ins = nodes_of_color[k]
                        box_label = len(boxes)
                        box_ins = []
                        for j in range(len(list_ins)):
                            box_ins.append(ins_set[list_ins[j]])
                        boxes.append((box_ins, t / steps, sumh[k]))
                        for label in local_ending_boxes:
                            edges.append((label, box_label))

                        #'''
                        # 1st order Trotter
                        if k > 0 :
                            edges.append((box_label - 1, box_label))
                        #'''
                        
                        next_local_ending_boxes.append(box_label)
                    local_ending_boxes = next_local_ending_boxes
                next_ending_boxes += next_local_ending_boxes

            ending_boxes = next_ending_boxes

        # Delete commutative edges
        """
        Clean up the result.
        Generate the minimal graph containing all relations.
        """
        G = nx.DiGraph()
        G.add_nodes_from(range(len(boxes)))
        s = 0
        while s < len(edges):
            edge = edges[s]
            if nx.has_path(G, edge[0], edge[1]):
                del edges[s]
                continue
            G.add_edge(*edge)
            s += 1
        s = 0
        while s < len(edges):
            edge = edges[s]
            G.remove_edge(*edge)
            if nx.has_path(G, edge[0], edge[1]):
                del edges[s]
                continue
            G.add_edge(*edge)
            h1 = boxes[edge[0]][2]
            h2 = boxes[edge[1]][2]
            if TIHamiltonian.commutativity_test(h1, h2):
                del edges[s]
                G.remove_edge(*edge)
                for i in range(len(edges)):
                    if edges[i][1] == edge[0]:
                        new_edge = (edges[i][0], edge[1])
                        if not nx.has_path(G, new_edge[0], new_edge[1]):
                            edges.append(new_edge)
                            G.add_edge(*new_edge)
                s -= 1
            s += 1

        logger.info(alignment)
        logger.info(sol_gvars)
        output_boxes = []
        for box in boxes:
            output_boxes.append(
                (
                    [((i, j), ins_lvars) for ((i, j), ins, h, ins_lvars) in box[0]],
                    box[1],
                )
            )
            logger.info(output_boxes[-1])
        for edge in edges:
            logger.info(edge)
        return (alignment, sol_gvars, output_boxes, edges)
    else:
        logger.info("No solution is found.")
        raise Exception("No solution is found.")
