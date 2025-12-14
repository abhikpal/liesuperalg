# -*-mode: python; -*-
# This file is part of the liesuperalg package for SageMath
# Copyright (C) 2023-2025 Abhik Pal
#
# liesuperalg is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# liesuperalg is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# liesuperalg. If not, see <https://www.gnu.org/licenses/>.

from sage.combinat.partition import PartitionsInBox
from sage.misc.latex import LatexExpr

# format string for partitions
def add(p, q):
    return [pi + qi for (pi, qi) in zip(p, q)]

def det (mu):
    return [mu + 1 for mi in mu]

def det_inv(mu):
    return [mi - 1 for mi in mu]

def dual(mu):
    return [-mi for mi in reversed(mu)]

def det_factor(weight):
    """factors the given gl_n weight into a power of determinant and a partition"""
    k = min(min(weight), 0)
    return k, [wi - k for wi in weight]

def common_det_factor(*weights):
    k = min([min(w) for w in weights])
    return tuple([k] + [[wi - k for wi in w] for w in weights])

def is_weakly_decreasing(seq):
    """Determine if the given sequence is weakly decreasing.
    """
    return all([seq(i) - seq(i + 1) >= 0 for i in range(len(seq) - 1)])

def zero_pad(part, n):
    return part + [0] * (n - len(part))

def B(N, h, w):
    for p in PartitionsInBox(h, w):
        if p.size() == N:
            yield p


# --- PRINTING METHODS ---

def pfmt(part, omit_zero=False, zero='.'):
    if omit_zero:
        return "[" + "  ".join(["{:>2}".format(pi) for pi in part if pi != 0]) + " ]"
    return "[" + "  ".join(["{:>2}".format(zero if pi == 0 else pi) for pi in part]) + " ]" # 4 * len(part)

def pp_sum(weight, tex=True):
    sum_str = ''
    for i in range(1, weight.m + weight.n + 1):
        if weight[i] == 0:
            continue

        if i > weight.m:
            sum_str += "-" if weight[i] > 0 else "+"
        else:
            sum_str += "+" if weight[i] > 0 else "-"

        sum_str += str(abs(weight[i])) if abs(weight[i]) != 1 else ""

        if tex:
            sum_str += "\\delta_{" if i > weight.m else "\\epsilon_{"
        else:
            sum_str += "d_{" if i > weight.m else "e_{"

        sum_str += (str(i - weight.m) if i > weight.m else str(i) ) + "}"

    if tex:
        return LatexExpr(sum_str)
    return sum_str

def pp_diagram(diagram, low=None, high=None, convention='sz', tex=True):
    if not low:
        low = min(diagram.keys()) - 1

    if not high:
        high = max(diagram.keys()) + 1

    idxs = []

    if convention == 'sz':
        tex_keys = {
            'B': "\\times", 
            'L': "\\succ",
            'R': "\\prec",
            '.': "\\ ",
        }

        str_keys = {
            'B': "x",
            'L': ">",
            'R': "<",
            '.': " ",
        }
    elif convention == 'bs':
        tex_keys = {
            'B': "\\vee", 
            'L': "\\times",
            'R': "\\circ",
            '.': "\\wedge",
        }

        str_keys = {
            'B': "v",
            'L': "x",
            'R': "o",
            '.': "^",
        }
    else:
        raise NotImplementedError

    if tex:
        for idx in range(low, high + 1):
            symb = diagram.get(idx, ".")
            idxs.append("\\overset{" + tex_keys[symb] + "}{" + str(idx) + "}")
        return LatexExpr("-".join(idxs))

    symbs = []
    for idx in range(low, high + 1):
        symb = str_keys[diagram.get(idx, '.')]
        sfmt = "{:<" + str(len(str(idx))) + "}"
        symbs.append(sfmt.format(symb))

    ints = " - ".join([str(idx) for idx in range(low, high + 1)])
    dstr = "   ".join(symbs) + "\n" + ints
    return dstr

def pp_kd(kd, fmt=None):
    """Pretty-print Kac-Dynkin labels"""
    l, m, r = kd

    if fmt is None:
        d = max([len(str(si)) for si in l + [m] + r])
        fmt = "{:" + str(d) + "d}"

    lstr = ", ".join([fmt.format(li) for li in l])
    rstr = ", ".join([fmt.format(ri) for ri in r])
    return "[" + lstr + "; " + fmt.format(m) + "; " + rstr + "]"

def pp_grid(data, r, fmt=str, default="", labels=True):
    max_len = max([len(fmt(dt)) for _, dt in data.items()] + [len(str(r)), 2])
    if labels:
        print("".rjust(max_len) + " t:".rjust(max_len))
        print("s:".rjust(max_len) + " " + \
              " ".join([str(ti).rjust(max_len) for ti in range(1, r + 1)]))
    for s in range(1, r + 1):
        if labels:
            print(str(s).rjust(max_len), end=' ')
        print(" ".join([fmt(data.get((s, t), default)).rjust(max_len)\
                        for t in range(1, r + 1)]))

def pp_nqc(nqc, r, oneline=False, **kwargs):
    if oneline:
        for s in range(1, r + 1):
            for t in range(s, r + 1):
                print(nqc[s, t], end="")
            print("/", end="")
        print()
        return None
    pp_grid(nqc, r, **kwargs)
