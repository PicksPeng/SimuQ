# Generate a QHD Hamiltonian for non-convex optimization
# The formulation comes from arXiv: 2303.01471
# Optimizing f(x, y) on [0, 1]^2
# H(α, β) = α (-1/2 ∇^2) + β f(x, y)
# ∇^2 = q^2 Σ_j (X_{x,j} + X_{y,j})
# x = 1/q Σ_j n_{x,j}
# y = 1/q Σ_j n_{y,j}

from simuq.environment import Qubit
from simuq.qsystem import QSystem


def GenQS(q, T, alpha, beta, f=lambda x, y: -x * x + y * y + x * y):
    qs = QSystem()
    x, y = [Qubit(qs) for i in range(q)], [Qubit(qs) for i in range(q)]
    H, xoper, yoper = 0, 0, 0
    for i in range(q):
        H += alpha * (-0.5) * q**2 * (x[i].X + y[i].X)
    for i in range(q):
        xoper += 1.0 / q * (x[i].I - x[i].Z) / 2
        yoper += 1.0 / q * (y[i].I - y[i].Z) / 2
    for i in range(q):
        H += beta * f(xoper, yoper)
    qs.add_evolution(H, T)
    return qs
