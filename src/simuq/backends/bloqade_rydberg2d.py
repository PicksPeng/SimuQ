def gen_bloqade_code(pos2d, clocks, pulse):
    n = len(pos2d)
    code = "using Bloqade\n"
    code += f"atoms = AtomList({ pos2d }) \n"
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


def gen_clocks(times):
    clocks = [0.0]
    for t in times:
        clocks.append(clocks[-1] + t)
    return clocks


def clean_as(alignment, sol_gvars, boxes):
    # pos = [0., 0., 0.] + sol_gvars
    pos = [0.0, 0.0, *sol_gvars]
    pos2d = []
    for i in range(len(pos) // 2):
        pos2d.append((pos[2 * i], pos[2 * i + 1]))
    n = len(pos2d)
    m = len(boxes)
    times = []
    pulse = [[[0.0 for i in range(m)] for j in range(n)] for k in range(3)]
    for evo_idx in range(m):
        box = boxes[evo_idx]
        times.append(box[1])
        for (i, j), ins, h, ins_lvars in box[0]:
            if i < n:
                pulse[0][i][evo_idx] = ins_lvars[0]
            else:
                pulse[1][i - n][evo_idx] = ins_lvars[0]
                pulse[2][i - n][evo_idx] = ins_lvars[1]

    return (pos2d, gen_clocks(times), pulse)


def transpile(alignment, sol_gvars, boxes, edges):
    code = gen_bloqade_code(*clean_as(alignment, sol_gvars, boxes))
    with open("transpiled.jl", "w") as f:
        print(code, file=f)
    return code
