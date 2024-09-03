import numpy as np
import math

def generate_as_ionq(qs, mach, aais, trotter_args, precision = 1e-5, verbose=0):
    n = mach.num_sites
    boxes = []
    edges = []

    ins_ind_2site = np.zeros((n, n), dtype=np.int32)
    cnt = n
    for i in range(n):
        for j in range(i + 1, n):
           ins_ind_2site[i, j] = cnt
           cnt += 1

    order = trotter_args["order"]
    steps = trotter_args["num"] * (1 << (order - 1))
    sequential = trotter_args["sequential"]

    P = {"X":0, "Y":1, "Z":2}
    ending_boxes = []
    for evo_index in range(len(qs.evos)):
        uni_coeff = np.zeros((n, 3))
        bi_coeff = np.zeros((n, n, 3, 3))
        h, t = qs.evos[evo_index]
        for (ham, c) in h.ham:
            if len(ham) == 2:
                pl = ham.to_list()
                (q0, p0), (q1, p1) = pl[0], pl[1]
                bi_coeff[q0, q1, P[p0], P[p1]] += c
            elif len(ham) == 1:
                [(q0, p0)] = ham.to_list()
                uni_coeff[q0, P[p0]] += c

        local_ending_boxes = ending_boxes

        if aais == "heisenberg":
            for q0 in range(n):
                for q1 in range(q0 + 1, n):
                    for k1 in range(3):
                        for k2 in range(3):
                            if k1 != k2 and np.abs(bi_coeff[q0, q1, k1, k2]) > precision:
                                raise Exception("The quantum system cannot be compiled to 'heisenberg' AAIS. There exist terms with mixing product Hamiltonian (i.e., 'q0.X * q1.Y'). Try 'two_pauli' AAIS.")

        for i in range(steps):
            new_boxes = []

            def add_box(box):
                new_boxes.append(len(boxes))
                boxes.append(box)

            for q0 in range(n):
                if np.sum(np.abs(uni_coeff[q0, 0:2])) > precision:
                    amp = np.linalg.norm(uni_coeff[q0, 0:2])
                    phase = math.atan2(uni_coeff[q0, 1], uni_coeff[q1, 0])
                    box = ([((q0, 0), [amp, phase])], t / steps)
                    add_box(box)
                if np.abs(uni_coeff[q0, 2]) > precision:
                    amp = uni_coeff[q0, 2]
                    box = ([((q0, 1), [amp])], t / steps)
                    add_box(box)
            
            for q0 in range(n):
                for q1 in range(q0 + 1, n):
                    if aais == "two_pauli" or "2pauli":
                        if np.sum(np.abs(bi_coeff[q0, q1])) > precision:
                            box = ([((ins_ind_2site[q0, q1], 0), bi_coeff[q0, q1].reshape(9))], t / steps)
                            add_box(box)
                    elif aais == "heisenberg":
                        for k in range(3):
                            if np.abs(bi_coeff[q0, q1, k, k]) > precision:
                                box = ([((ins_ind_2site[q0, q1], k), [bi_coeff[q0, q1, k, k]])], t / steps)
                                add_box(box)

            if sequential:
                if order == 1 or (order == 2 and i % 2 == 0):
                    for b in local_ending_boxes:
                        edges.append((b, new_boxes[0]))
                    for j in range(1, len(new_boxes)):
                        edges.append((new_boxes[j - 1], new_boxes[j]))
                    local_ending_boxes = [new_boxes[-1]]
                elif order == 2 and i % 2 == 1:
                    for b in local_ending_boxes:
                        edges.append((b, new_boxes[-1]))
                    for j in range(len(new_boxes) - 1, 0, -1):
                        edges.append((new_boxes[j], new_boxes[j - 1]))
                    local_ending_boxes = [new_boxes[0]]
            else:
                for j in new_boxes:
                    for b in local_ending_boxes:
                        edges.append((b, j))
                local_ending_boxes = new_boxes

        ending_boxes = local_ending_boxes
    
    return boxes, edges