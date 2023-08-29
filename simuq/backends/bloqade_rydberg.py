def gen_bloqade_code_constant(pos, clocks, pulse):
    n = len(pos)
    code = "using Bloqade\n"
    pos_1d = list(map(lambda x: (x,), pos))
    code += f"atoms = AtomList({ pos_1d }) \n"
    code_delta = "delta = ["
    code_omega = "omega = ["
    code_phi = "phi = ["
    for i in range(n):
        code_delta += f"piecewise_constant(clocks = { clocks }, values = { pulse[0][i] }), "
        code_omega += f"piecewise_constant(clocks = { clocks }, values = { pulse[1][i] }), "
        code_phi += f"piecewise_constant(clocks = { clocks }, values = { pulse[2][i] }), "
    code_omega += "]\n"
    code_phi += "]\n"
    code_delta += "]\n"
    code += code_omega + code_phi + code_delta
    code += "h = rydberg_h(atoms; Δ = delta, Ω = omega, ϕ = phi)"
    return code


def gen_bloqade_code_linear(pos, clocks, pulse):
    n = len(pos)
    code = "using Bloqade\n"
    pos_1d = list(map(lambda x: (x,), pos))
    if len(pulse[0][0]) == 1:
        for i in range(3):
            for j in range(n):
                pulse[i][j].append(pulse[i][j][-1])
    else:
        for i in range(3):
            for j in range(n):
                k = (pulse[i][j][-1] - pulse[i][j][-2]) / (clocks[-2] - clocks[-3])
                pulse[i][j].append(pulse[i][j][-1] + k * (clocks[-1] - clocks[-2]))
    code += f"atoms = AtomList({ pos_1d }) \n"
    code_delta = "delta = ["
    code_omega = "omega = ["
    code_phi = "phi = ["
    for i in range(n):
        code_delta += f"piecewise_linear(clocks = { clocks }, values = { pulse[0][i] }), "
        code_omega += f"piecewise_linear(clocks = { clocks }, values = { pulse[1][i] }), "
        code_phi += f"piecewise_linear(clocks = { clocks }, values = { pulse[2][i] }), "
    code_omega += "]\n"
    code_phi += "]\n"
    code_delta += "]\n"
    code += code_omega + code_phi + code_delta
    code += "h = rydberg_h(atoms; Δ = delta, Ω = omega, ϕ = phi)"
    return code


def gen_clocks(times):
    clocks = [0]
    for t in times:
        clocks.append(clocks[-1] + t)
    return clocks


def clean_as(alignment, sol_gvars, boxes):
    pos = [0, *sol_gvars]
    n = len(pos)
    m = len(boxes)
    times = []
    pulse = [[[0 for i in range(m)] for j in range(n)] for k in range(3)]
    for evo_idx in range(m):
        box = boxes[evo_idx]
        times.append(box[1])
        for (i, j), ins, h, ins_lvars in box[0]:
            if i < n:
                pulse[0][i][evo_idx] = ins_lvars[0]
            else:
                pulse[1][i - n][evo_idx] = ins_lvars[0]
                pulse[2][i - n][evo_idx] = ins_lvars[1]

    return (pos, gen_clocks(times), pulse)


def transpile(alignment, sol_gvars, boxes, edges, inter_order=0):
    if inter_order == 0:
        code = gen_bloqade_code_constant(*clean_as(alignment, sol_gvars, boxes))
    elif inter_order == 1:
        code = gen_bloqade_code_linear(*clean_as(alignment, sol_gvars, boxes))
    # with open('transpiled.jl', 'w') as f :
    #    print(code, file = f)
    return code
