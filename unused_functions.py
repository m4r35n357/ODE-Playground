def t_sqrt(r, u, k):
    if k == 0:
        return sqrt(u[0])
    else:
        rt = 0.0
        if k % 2 == 1:
            for j in range(1, (k - 1) // 2 + 1):
                rt += r[j] * r[k - j]
            return (u[k] - 2.0 * rt) / (2.0 * u[0])
        else:
            for j in range(1, (k - 2) // 2 + 1):
                rt += r[j] * r[k - j]
            return (u[k] - 2.0 * rt - r[k // 2]**2) / (2.0 * u[0])


def t_ln(l, u, k):
    if k == 0:
        return exp(u[0])
    else:
        return (u[k] - ddot(u, l, k)) / u[0]


def ad_ln(u):
    n = len(u)
    ln_jet = jet_0(n)
    for k in range(n):
        ln_jet[k] = t_ln(ln_jet, u, k)
    return ln_jet


def ad_sqrt(u):
    n = len(u)
    sqrt_jet = jet_0(n)
    for k in range(n):
        sqrt_jet[k] = t_sqrt(sqrt_jet, u, k)
    return sqrt_jet



