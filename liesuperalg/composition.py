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

"""Composition Factors of Kac Modules
The compution here is based on Def # 3.7, Eq (3.12), and the main
Theorem 1.1 in [su2006:composition].
"""

from .weight import rho
from .weight import Weight

def lowering_op(theta, w_lambda):
    """Lowering operator ([su2006:composition] Eq (3.12))

    """
    if len(theta) != w_lambda.adeg:
        raise ValueError

    wrho = w_lambda.rho()
    data = wrho._data
    for ts, mns in zip(theta, w_lambda.atypical_roots):
        ms, ns = mns

        k = wrho[ms] - 1

        while ts > 0:
            x = w_lambda._diagram.get(k, None)
            # print(f"checking {k} in {w_lambda._S}")
            if x is None:
                ts = ts - 1
                data[ms - 1] = k
                data[ns - 1] = k
            k = k - 1

    L = sorted(data[:w_lambda.m], reverse=True)
    R = sorted(data[w_lambda.m:])
    return Weight(L, R) - rho(w_lambda.m, w_lambda.n)

def kac_composition_factors(w_lambda):
    """Iterator for composition factors corresponding to the Kac module of `w_lambda`.

    """
    def Theta(r):
        """Generate all tuples (t1, ..., tr) where
        0 <= t_s <=s.
        """
        if r == 1:
            yield (0,)
            yield (1,)

        for s in range(r + 1):
            for ttuple in Theta(r - 1):
                yield ttuple + (s,)

    def C1(ttuple, nqc, r):
        """condition C1 ([su2006:composition] Def 3.7)"""
        for s, ts in enumerate(ttuple, start=1):
            if 1 <= ts and ts < s:
                if nqc[s - ts, s] != 'c':
                    continue
                return False
        return True

    def C2(ttuple, nqc, r):
        """condition C2 ([su2006:composition] Def 3.7)"""
        for s, ts in enumerate(ttuple, start=1):
            if 1 < ts and ts <= s:
                if nqc[s + 1 - ts, s] != 'n':
                    continue
                return False
        return True

    def C3(ttuple, nqc, r, verbose=False):
        """condition C3 ([su2006:composition] Def 3.7)"""
        def _c3_i_b_find_pp(s, p):
            # find p' satisfying criteria in C3.i.b
            for pp in range(p + 1, s):
                if ttuple[pp - 1] > pp - p and nqc[p, pp] == 'q':
                    return pp
            return None

        for s, ts in enumerate(ttuple, start=1):
            if ts < 2:
                continue

            for p in range(s + 1 - ts, s):
                tp = ttuple[p - 1]

                # print(f"ts = {ts}, s = {s}, p = {p}")

                if nqc[p, s] == 'c':
                    c3_i_a = (1 <= tp) and (tp <= ts - s + p)
                    if c3_i_a:
                        continue

                    c3_i_b1 = tp == 0
                    c3_i_b2 = (_c3_i_b_find_pp(s, p) is not None)
                    c3_i_b = c3_i_b1 and c3_i_b2
                    if c3_i_b:
                        continue

                if nqc[p, s] == 'q':
                    c3_ii_a = tp == 0
                    if c3_ii_a:
                        continue

                    c3_ii_b1 = s + 1 - ts < p
                    c3_ii_b2 = 1 <= tp
                    c3_ii_b3 = tp < ts - s + p
                    c3_ii_b = c3_ii_b1 and c3_ii_b2 and c3_ii_b3
                    if c3_ii_b:
                        continue

                if nqc[p, s] == 'n':
                    c3_iii = (0 <= tp) and (tp < ts - s + p)
                    if c3_iii:
                        continue

                return False
        return True

    def Theta_lambda(nqc, r):
        """Generate all tuples in Theta^lambda ([su2006:composition], Def 3.7)

        :param nqc: the nqc-table of a w_lambda
        :type nqc: dict

        :param r: degree of atypicality of a w_lambda
        :type r: int
        """
        valid_tuple = lambda tt: C1(tt, nqc, r) and C2(tt, nqc, r) and C3(tt, nqc, r)
        for ttuple in filter(valid_tuple, Theta(r)):
            yield ttuple

    # Main Theorem 1.1 in [su2006:composition]
    for theta in Theta_lambda(w_lambda.nqc, w_lambda.adeg):
        yield lowering_op(theta, w_lambda)

