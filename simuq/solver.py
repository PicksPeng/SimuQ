import scipy.optimize as opt
import numpy as np
import networkx as nx
from simuq.hamiltonian import TIHamiltonian
import math

def locate_switch_evo(mach, evo_index) :
    return mach.num_gvars + evo_index * (mach.num_inss + mach.num_lvars)

def switch_term(mach, evo_index, ins_index, mc) :
    return (lambda x : x[locate_switch_evo(mach, evo_index) + ins_index] * mc.exp_eval(x[:mach.num_gvars], x[locate_switch_evo(mach, evo_index) + mach.num_inss : locate_switch_evo(mach, evo_index + 1)]))

def switch_fun(mach, evo_index, ins_index) :
    return (lambda x : x[locate_switch_evo(mach, evo_index) + ins_index])

def locate_nonswitch_evo(mach, evo_index) :
    return mach.num_gvars + evo_index * mach.num_lvars

def non_switch_term(mach, evo_index, ins_index, mc) :
    return (lambda x : mc.exp_eval(x[:mach.num_gvars], x[locate_nonswitch_evo(mach, evo_index) : locate_nonswitch_evo(mach, evo_index + 1)]))

def solve_aligned(ali, qs, mach) :
    print(ali)
    
    eqs = []
    nvar = mach.num_gvars + len(qs.evos) * (mach.num_inss + mach.num_lvars)
    for evo_index in range(len(qs.evos)) :
        (h, t) = qs.evos[evo_index]
        mark = [[0 for ins in line.inss] for line in mach.lines]
        for (tprod, tc) in h.ham :
            eq = (lambda c : lambda x : -c)(tc)
            for i in range(len(mach.lines)) :
                line = mach.lines[i]
                for j in range(len(line.inss)) :
                    ins = line.inss[j]
                    for (mprod, mc) in ins.h.ham :
                        prod_match = True
                        for k in range(qs.num_sites) :
                            if tprod[k] != mprod[ali[k]] :
                                prod_match = False
                                break
                        if prod_match :
                            eq = (lambda eq_, f_ : lambda x : eq_(x) + f_(x))(eq, switch_term(mach, evo_index, ins.index, mc))
                            mark[i][j] = 1
                            break
            eqs.append(eq)
        for i in range(len(mach.lines)) :
            line = mach.lines[i]
            for j in range(len(line.inss)) :
                ins = line.inss[j]
                if mark[i][j] == 0 :
                    eqs.append(switch_fun(mach, evo_index, ins.index))

    #if len(eqs) < nvar :
    #    eqs = eqs + [(lambda x : 0) for i in range(nvar - len(eqs))]
    f = (lambda eqs_ : lambda x : [(lambda i_ : eqs_[i_](x))(i) for i in range(len(eqs_))])(eqs)
    var_lb = -np.inf
    var_ub = np.inf
    lbs = [var_lb for i in range(mach.num_gvars)]
    ubs = [var_ub for i in range(mach.num_gvars)]
    init = [0 for i in range(mach.num_gvars)]
    for i in range(len(qs.evos)) :
        lbs += [0 for j in range(mach.num_inss)] + [var_lb for j in range(mach.num_lvars)]
        ubs += [1 for j in range(mach.num_inss)] + [var_ub for j in range(mach.num_lvars)]
        init += [0.5 for j in range(mach.num_inss)] + [0 for j in range(mach.num_lvars)]

    sol_detail = opt.least_squares(f, init, bounds = (lbs, ubs))
    sol = sol_detail.x

    print(np.linalg.norm(f(sol)))
    if np.linalg.norm(f(sol)) > 1e-3 :
        return False

    lbs = [var_lb for i in range(mach.num_gvars)]
    ubs = [var_ub for i in range(mach.num_gvars)]
    for i in range(len(qs.evos)) :
        lbs += [var_lb for j in range(mach.num_lvars)]
        ubs += [var_ub for j in range(mach.num_lvars)]

    eqs = []
    nvar = mach.num_gvars + len(qs.evos) * mach.num_lvars
    initvars = sol[:mach.num_gvars].tolist()
    switch = [[[0 for ins in line.inss] for line in mach.lines] for evo_index in range(len(qs.evos))]
    for evo_index in range(len(qs.evos)) :
        for i in range(len(mach.lines)) :
            line = mach.lines[i]
            for j in range(len(line.inss)) :
                ins = line.inss[j]
                if abs(sol[locate_switch_evo(mach, evo_index) + ins.index]) < 1e-3 :
                    switch[evo_index][i][j] = 0
                else :
                    switch[evo_index][i][j] = 1
        initvars += sol[locate_switch_evo(mach, evo_index) + mach.num_inss : locate_switch_evo(mach, evo_index + 1)].tolist()
    for evo_index in range(len(qs.evos)) :
        (h, t) = qs.evos[evo_index]
        for (tprod, tc) in h.ham :
            eq = (lambda c : lambda x : -c)(tc)
            for i in range(len(mach.lines)) :
                line = mach.lines[i]
                for j in range(len(line.inss)) :
                    ins = line.inss[j]
                    if switch[evo_index][i][j] == 1 :
                        for (mprod, mc) in ins.h.ham :
                            prod_match = True
                            for k in range(qs.num_sites) :
                                if tprod[k] != mprod[ali[k]] :
                                    prod_match = False
                                    break
                            if prod_match :
                                eq = (lambda eq_, f_ : lambda x : eq_(x) + f_(x))(eq, non_switch_term(mach, evo_index, ins.index, mc))
                                break
            eqs.append(eq)

    f = (lambda eqs_ : lambda x : [(lambda i_ : eqs_[i_](x))(i) for i in range(len(eqs_))])(eqs)
    var_lb = -np.inf
    var_ub = np.inf
    lbs = [var_lb for i in range(mach.num_gvars)]
    ubs = [var_ub for i in range(mach.num_gvars)]
    for i in range(len(qs.evos)) :
        lbs += [var_lb for j in range(mach.num_lvars)]
        ubs += [var_ub for j in range(mach.num_lvars)]

    global gsol
    global gswitch
    global alignment
    alignment = ali
    if len(initvars) == 0 :
        gsol = np.array([])
        gswitch = switch
        if np.linalg.norm(f([])) > 1e-5 :
            return False
        return True
    sol_detail = opt.least_squares(f, initvars, bounds = (lbs, ubs))
    sol = sol_detail.x

    gsol = sol
    gswitch = switch
    
    print(np.linalg.norm(f(sol)))
    if np.linalg.norm(f(sol)) > 1e-5 :
        return False

    print(switch)
    print(sol)
    return True


def align(i, ali, qs, mach) :
    if i == qs.num_sites :
        if solve_aligned(ali, qs, mach) :
            return True
        return False
    for x in range(mach.num_sites) :
        available = True
        ali[i] = x
        for j in range(i) :
            if ali[j] == x :
                available = False
                break
        if available == False :
            continue
        for (h, t) in qs.evos :
            for (tprod, tc) in h.ham :
                found = False
                for line in mach.lines :
                    for ins in line.inss :
                        for (mprod, mc) in ins.h.ham :
                            prod_partial_match = True
                            for k in range(i + 1) :
                                if tprod[k] != mprod[ali[k]] :
                                    prod_partial_match = False
                            if prod_partial_match :
                                found = True
                                break
                        if found :
                            break
                    if found :
                        break
                if not found :
                    available = False
                    break
            if not available :
                break
        if available :
            if align(i + 1, ali, qs, mach) :
                return True
    return False


def find_sol(qs, mach, ali = []) :
    if ali == [] :
        return align(0, [0 for i in range(qs.num_sites)], qs, mach)
    else :
        return solve_aligned(ali, qs, mach)


def generate_as(qs, mach, trotter_step = 4) :
    if find_sol(qs, mach) :
        sol = gsol
        switch = gswitch
        sol_gvars = sol[:mach.num_gvars]
        boxes = []
        edges = []
        ending_boxes = []
        for evo_index in range(len(qs.evos)) :
            (h, t) = qs.evos[evo_index]
            next_ending_boxes = []
            coloring = [i for i in range(mach.num_sites)]
            sol_lvars = sol[locate_nonswitch_evo(mach, evo_index) : locate_nonswitch_evo(mach, evo_index + 1)].tolist()
            # Detach by touched sites
            first_touch = [[-1 for i in range(len(line.inss))] for line in mach.lines]
            for i in range(len(mach.lines)) :
                line = mach.lines[i]
                for j in range(len(line.inss)) :
                    ins = line.inss[j]
                    if switch[evo_index][i][j] == 1 :
                        t_sites = (ins.h.exp_eval(sol_gvars, sol_lvars)).touched_sites()
                        for k in range(mach.num_sites) :
                            if t_sites[k] == 1 :
                                if first_touch[i][j] == -1 :
                                    first_touch[i][j] = k
                                else :
                                    if coloring[k] != coloring[first_touch[i][j]] :
                                        c = coloring[k]
                                        for l in range(mach.num_sites) :
                                            if coloring[l] == c :
                                                coloring[l] = coloring[first_touch[i][j]]
            color_part = []
            for k in range(mach.num_sites) :
                ins_set = []
                for i in range(len(mach.lines)) :
                    line = mach.lines[i]
                    for j in range(len(line.inss)) :
                        ins = line.inss[j]
                        if switch[evo_index][i][j] == 1  and  coloring[first_touch[i][j]] == k :
                            ins_lvars = []
                            for l in range(len(ins.vars_index)) :
                                ins_lvars.append(sol_lvars[ins.vars_index[l]])
                            ins_set.append(((i, j), ins, ins.exp_eval(sol_gvars, sol_lvars), ins_lvars))
                if ins_set != [] :
                    color_part.append(ins_set)
            # Detach if commute with all others
            for i in range(len(color_part)) :
                ins_set = color_part[i]
                j = 0
                while j < len(ins_set) :
                    all_commute = True
                    for k in range(len(ins_set)) :
                        if j == k :
                            continue
                        if not TIHamiltonian.commutativity_test(ins_set[j][2], ins_set[k][2], \
                                                                ins_set[j][1].prop == 'derived' or ins_set[k][1].prop == 'derived') :
                            all_commute = False
                            break
                    if all_commute :
                        color_part.append([ins_set[j]])
                        del ins_set[j]
                        j -= 1
                    j += 1
            # Trotterization
            for ins_set in color_part :
                local_ending_boxes = ending_boxes
                n = len(ins_set)
                G = nx.Graph()
                G.add_nodes_from(range(n))
                for i in range(n) :
                    for j in range(n) :
                        if i == j :
                            continue
                        if (ins_set[i][0][0] == ins_set[j][0][0]) \
                           or ((ins_set[i][1].prop == 'derived' or ins_set[j][1].prop == 'derived') \
                               and (not TIHamiltonian.commutativity_test(ins_set[i][2], ins_set[j][2], True))) :
                            G.add_edge(i, j)
                col = nx.coloring.greedy_color(G, strategy = 'largest_first')
                nodes_of_color = [[] for i in set(col.values())]
                for i in range(n) :
                    nodes_of_color[col[i]].append(i)
                # Check if the partitions are pair-wise commute
                sumh = []
                for i in range(len(nodes_of_color)) :
                    h = TIHamiltonian.empty(mach.num_sites)
                    for j in range(len(nodes_of_color[i])) :
                        h = h + ins_set[nodes_of_color[i][j]][2]
                    sumh.append(h)
                all_commute = True
                for i in range(len(nodes_of_color)) :
                    for j in range(i) :
                        if TIHamiltonian.commutativity_test(sumh[i], sumh[j]) == False :
                            all_commute = False
                            break
                    if not all_commute :
                        break
                if all_commute :
                    steps = 1
                else :
                    steps = trotter_step
                for i in range(steps) :
                    next_local_ending_boxes = []
                    for k in range(len(nodes_of_color)) :
                        list_ins = nodes_of_color[k]
                        box_label = len(boxes)
                        box_ins = []
                        for j in range(len(list_ins)) :
                            box_ins.append(ins_set[list_ins[j]])
                        boxes.append((box_ins, t / steps, sumh[k]))
                        for label in local_ending_boxes :
                            edges.append((label, box_label))
                        next_local_ending_boxes.append(box_label)
                    local_ending_boxes = next_local_ending_boxes
                next_ending_boxes += next_local_ending_boxes

            ending_boxes = next_ending_boxes

        # Delete commutative edges
        G = nx.DiGraph()
        G.add_nodes_from(range(len(boxes)))
        s = 0
        while s < len(edges) :
            edge = edges[s]
            if nx.has_path(G, edge[0], edge[1]) :
                del edges[s]
                continue
            G.add_edge(*edge)
            s += 1
        s = 0
        while s < len(edges) :
            edge = edges[s]
            G.remove_edge(*edge)
            if nx.has_path(G, edge[0], edge[1]) :
                del edges[s]
                continue
            G.add_edge(*edge)
            h1 = boxes[edge[0]][2]
            h2 = boxes[edge[1]][2]
            if TIHamiltonian.commutativity_test(h1, h2) :
                del edges[s]
                G.remove_edge(*edge)
                for i in range(len(edges)) :
                    if edges[i][1] == edge[0] :
                        new_edge = (edges[i][0], edge[1])
                        if not nx.has_path(G, new_edge[0], new_edge[1]) :
                            edges.append(new_edge)
                            G.add_edge(*new_edge)
                s -= 1
            s += 1


        print()
        print(alignment)
        print(sol_gvars)
        for box in boxes :
            print(([((i, j), ins_lvars) for ((i, j), ins, h, ins_lvars) in box[0]], box[1]))
        for edge in edges :
            print(edge)
    else :
        print('No solution is found.')


