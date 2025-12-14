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

from .utils import add
from .utils import det_factor

class Weight:
    # lambda
    #   = lambda_{0, 1} epsilon_1 + ... + lamdba_{0, m} epsilon_m - ...
    #         - lambda_{1, 1} delta_1 - ... - mu_{1, n} delta_n                             # sum notation
    #   = (lambda_{0 , 1}, ..., lambda_{0, m} | lambda_{1, 1}, ..., lambda_{1, n})          # default notation
    #   = Weight([lambda_{0, 1}, ..., lambda_{0, m}], [lambda_{1, 1}, ..., lambda_{1, n}])  # CONSTRUCTOR
    #
    def __init__(self, L, R, setup=False):
        self.L = L
        self.R = R
        self.m = len(L)
        self.n = len(R)
        self._data = L + R

        self._is_diag_setup = False
        if setup:
            self._diag_setup()

        self._atyp_roots = None
        self._lens = None
        self._kd = None
        self._nqc = None
        self._cr = None
        self._scr = None

    def _diag_setup(self):
        srho = self.rho()
        self._SL = set(srho.L)
        self._SR = set(srho.R)
        self._S = self._SL.union(self._SR)
        self._SB = self._SL.intersection(self._SR)

        diag = dict()
        for i in self._SB:
            diag[i] = 'B'
        for i in self._SR.difference(self._SB):
            diag[i] = 'R'
        for i in self._SL.difference(self._SB):
            diag[i] = 'L'
        self._diagram = diag

        self._is_diag_setup = True

    def _needs_diagram(f):
        """Decorator for functions that require a diagram to be setup"""
        def _f(self, *args, **kwargs):
            if not self._is_diag_setup:
                self._diag_setup()

            return f(self, *args, **kwargs)
        return _f

    @property
    def coeff_L(self):
        return self.L

    @property
    def coeff_R(self):
        return [-ri for ri in self.R]

    def det_factor(self):
        """If self = lambda, then factor
        V(lambda) = det(V0)^k0 * S_{part0}V0 * det(V1)^k1 * S_{part1}V1

        where ki is an integer and part0, part1 are partitions.

        """
        k0, part0 = det_factor(self.coeff_L)
        k1, part1 = det_factor(self.coeff_R)
        return ([k0, k1], [part0, part1])

    @property
    def kac_dynkin(self):
        """Return Kac-Dynkin labels"""

        # Definition of Kac-Dynkin labels from (3.5) in [hkv1992:composition]
        if self._kd == None:
            left = [self[i] - self[i + 1] for i in range(1, self.m)]

            # based on our sign convention, the weight
            #
            #     \Lambda = \sum_i \mu_i \epsilon_i + \sum_a \nu_a \delta_a
            #
            # in [hkv1992:composition] in par. before eq (3.5) becomes
            #
            #     \Lambda = (\mu_1, ..., \mu_m | -\nu_1, ..., -\nu_n).
            #
            # Hence, we need to change the sign in the middle term of
            # the Kac-Dynkin label
            mid = self[0, self.m] - self[1, 1]
            right = [self[i + 1] - self[i] for i in range(self.m + 1, self.m + self.n)]
            self._kd = (left, mid, right)
        return self._kd

    def rho(self, *args, **kwargs):
        """Return the rho translate"""
        return self + rho(self.m, self.n, *args, **kwargs)

    @_needs_diagram
    def is_dominant(self):
        """Checks if the weight is dominant."""
        # A weight is dominant iff number of < + x = m and > + x = n.
        # See Gruson-Serganova Sec 6 after the definition of a weight
        # diagram. [https://arxiv.org/pdf/0906.0918.pdf]
        # return (len(self._SL) == self.m) and (len(self._SR) == self.n)
        dsc_L = [li >= lj for li, lj in zip(self.L[:-1], self.L[1:])]
        asc_R = [ri <= rj for ri, rj in zip(self.R[:-1], self.R[1:])]
        return self.L[-1] >= 0 and self.R[0] >= 0 and all(dsc_L) and all(asc_R)

    def is_positive(self):
        """Checks if the weight is positive"""
        # defn due to Serganova "Kazhdan-Lusztig Polynomials and
        # Character Formula" (1996). See Remark 1.2
        coeffs = self.coeff_L + self.coeff_R
        all_sums_non_neg = all(sum(coeffs[:i]) >= 0 for i in range(self.m + self.n - 1))
        zero_sum = sum(coeffs) == 0
        return zero_sum and all_sums_non_neg

    @property
    @_needs_diagram
    def adeg(self):
        """Degree of atypicality """
        # [sz2012:generalised] after (5.3)
        return len(self._SB)

    @property
    @_needs_diagram
    def diagram(self):
        return self._diagram

    @property
    def atypical_roots(self):
        if self._atyp_roots is None:
            srho = self.rho()
            self._atyp_roots = []
            for i in range(self.m, 0, -1):
                for j in range(1, self.n + 1):
                    if srho[i] == srho[self.d(j)]:
                        self._atyp_roots.append((i, self.d(j)))
        return self._atyp_roots

    @property
    def atypicality_matrix(self):
        srho = self.rho()
        A = []
        for ai in srho.L:
            A.append([])
            for aj in srho.R:
                A[-1].append(ai - aj)
        return A

    @property
    @_needs_diagram
    def lens(self):
        if self._lens is None:
            sb = sorted(self._SB)
            self._lens = dict()
            for t in range(1, self.adeg + 1):
                for s in range(1, t + 1):
                    interval = set(range(sb[s - 1], sb[t - 1] + 1))
                    self._lens[(s, t)] = len(interval.difference(self._S))
        return self._lens

    # To account for Def 3.6 in [sz2007:character].
    dist = lens

    @property
    @_needs_diagram
    def nqc(self):
        if self._nqc is None:
            self._nqc = dict()
            for st, l in self.lens.items():
                s, t = st
                if l > t - s:
                    self._nqc[s, t] = 'n'
                elif l == t - s:
                    self._nqc[s, t] = 'q'
                else:
                    self._nqc[s,  t] = 'c'

        return self._nqc

    @property
    def cr(self):
        """c-relation between atypical roots"""
        if self._cr is None:
            # Def P:3.6 and eq p:(3.23) in [sz2007:character]
            self._cr = dict()
            for s in range(1, self.adeg + 1):
                self._cr[s, s] = True
                for t in range(s + 1, self.adeg + 1):
                    self._cr[s, t] = self.nqc[s, t] == 'c'
        return self._cr

    @property
    def scr(self):
        """strong c-relation between atypical roots"""
        if self._scr is None:
            self._scr = dict()
            for s in range(1, self.adeg + 1):
                # Def P:3.6 and eq p:(3.23) in [sz2007:character]
                self._scr[s, s] = self.cr[s, s]
                for t in range(s + 1, self.adeg + 1):
                    self._scr[s, t] = self.cr[s, t] and self.cr[s, t - 1]
        return self._scr

    @property
    def height(self):
        h = []
        for s, gamma in enumerate(self.atypical_roots, start=1):
            ms, dot_ns = gamma
            ns = dot_ns - self.m
            h.append(self[ms] - ns + s)
        return h

    @property
    def atyp(self):
        atyp = []
        srho = self.rho()
        for gamma in self.atypical_roots:
            ms, _ = gamma
            atyp.append(srho[ms])
        self._atyp = atyp
        return self._atyp

    @property
    def typ(self):
        srho = self.rho()
        mss = [gamma[0] for gamma in self.atypical_roots]
        ns_dots = [gamma[1] for gamma in self.atypical_roots]

        typ_L = [srho[0, s] for s in range(1, self.m + 1) if not (s in mss)]
        typ_R = [srho[1, s] for s in range(1, self.n + 1) if not (self.d(s) in ns_dots)]

        return [typ_L, typ_R]

    def d(self, i):
        """Internal method to get dot-indexing"""
        return i + self.m

    def dot(self):
        """Expose the internal dot-indexing method.
        """
        return lambda i: self.d(i)

    def respects_scr(self, sigma):
        """Check whether sigma respects the order of all strongly
        c-related atypical roots. That is,

            sigma_inv(s) < sigm_inv(t)

        whenever gamma_s and gamma_t are strongly c-related. See Def
        P:3.16 in [sz2007:character].

        """
        r = len(sigma.tuple())
        if r is not self.adeg:
            raise ValueError(f"Expected sigma in Sym_{self.adeg} got Sym_{r}")

        sigma_inv = sigma.inverse()
        is_scr_idx = lambda st: self.scr[st[0], st[1]] and (st[0] != st[1])
        return all([sigma_inv(st[0]) < sigma_inv(st[1]) \
                    for st in filter(is_scr_idx, self.scr.keys())])

    def atyp_dot_action(self, sigma):
        r = len(sigma.tuple())
        if r is not self.adeg:
            raise ValueError(f"Expected sigma in Sym_{self.adeg} got Sym_{r}")

        srho = self.rho()
        permuted_L = srho.L
        permuted_R = srho.R
        sigma_inv = sigma.inverse()

        for mns, perm_root in zip(self.atypical_roots, sigma_inv(self.atyp)):
            ms, ns_dot = mns
            ns = ns_dot - self.m
            permuted_L[ms - 1] = perm_root
            permuted_R[ns - 1] = perm_root

        return Weight(permuted_L,  permuted_R) - rho(self.m, self.n)

    def __getitem__(self, key):
        if type(key) == type((0, 1)):
            parity, idx = key
            return self._data[parity * self.m + idx - 1]

        # type casting to int to convert from sage.rings.integer.Integer
        return self._data[int(key) - 1]

    def __add__(self, other):
        compatible = (self.n == other.n) and (self.m == other.m)
        if not compatible:
            raise TypeError("unsupported operand (+) for weights"
                            f"in gl({self.m} | {self.n}) and gl({other.m} | {other.n})")

        return Weight(add(self.L, other.L), add(self.R, other.R))

    def __sub__(self, other):
        return self + (- other)

    def __neg__(self):
        return Weight([-li for li in self.L], [-ri for ri in self.R])

    def __eq__(self, other):
        return self._data == other._data

    def __le__(self, other):
        """Partial order on weights as defined in [sz2007:character] eq (3.10)
        """
        adeg_eq = self.adeg == other.adeg
        typ_eq = self.typ == other.typ
        atyp_leq = all([ai <= bi for (ai, bi) in zip(self.atyp, other.atyp)])
        return adeg_eq and typ_eq and atyp_leq

    def __lt__(self, other):
        return (self <= other) and (self != other)

    def pp(self, fmt):
        l_str = ", ".join([fmt.format(li) for li in self.L])
        r_str = ", ".join([fmt.format(ri) for ri in self.R])
        return "(" + l_str + " | " + r_str + ")"

    def _latex_(self):
        tstr = ""
        l_str = ", ".join([str(li) for li in self.L])
        r_str = ", ".join([str(ri) for ri in self.R])
        tstr = "\\left(" + l_str + " \\mid " + r_str + "\\right)"
        return tstr

    def __str__(self):
        strs = [str(li) for li in self.L] + [str(ri) for ri in self.R]
        sz = max(2, max([len(s) for s in strs]))
        fmt = "{:>" + str(sz) + "}"
        return self.pp(fmt)

    def __repr__(self):
        return f"gl({self.m} | {self.n}) weight {str(self)}"

    def __hash__(self):
        return hash((tuple(self.L), tuple(self.R)))

# --- WEIGHTS: MISC --------------------------------------------------
def rho(m, n, *args, **kwargs):
    return Weight(list(range(m, 0, -1)), list(range(1, n + 1, 1)), *args, **kwargs)

def one(m, n, *args, **kwargs):
    return Weight([1] * m, [1] * n)



