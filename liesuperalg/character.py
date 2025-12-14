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
"""Character Formula for Irreducible Representations

The main algorithm is based on [sz2007:character].
"""

from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.permutation import Permutation
from sage.groups.perm_gps.permgroup_named import SymmetricGroup

from .utils import is_weakly_decreasing
from .weight import Weight

def length(w_lambda, w_mu):
    """Length between two weights"""
    # eq (3.15) in [sz2007:character]
    return sum(w_lambda.height) - sum(w_mu.height)

def S(w_lambda, w_mu):
    """Iterator for the set S^{lamdba, mu}.

    See Def P:3.16 [sz2007:character]
    """
    for sigma in SymmetricGroup(w_lambda.adeg):
        if not w_lambda.respects_scr(sigma):
            continue
        if not (w_mu <= w_lambda.atyp_dot_action(sigma)):
            continue
        yield sigma

def gen_KL(w_lambda, w_mu, q):
    """Calculate the generalized Kazhdan-Lusztig polynomial
    K_{lambda, mu}(q).

    """
    S_terms = 0
    for sigma in S(w_lambda, w_mu):
        bruhat_len = Permutation(sigma).number_of_inversions()
        S_terms += q**(- 2 * bruhat_len)
    return (q ** length(w_lambda, w_mu)) * S_terms

def mult_kac_in_irrd(w_lambda, w_mu):
    """Returns the multiplicity of ch[K(mu)] in ch[L(lambda)]

    In the notation of [sz2007:character] this multiplicity is equal
    to m(weight)_mu = |S^{weight, mu}| (eq P:(4.1))
    """
    # if lambda and mu aren't in the same block, return 0.
    if w_lambda.typ != w_mu.typ:
        return 0

    m = sum([1 for s in S(w_lambda, w_mu)])
    return (-1) ** length(w_lambda, w_mu) * m

def phi(w_lambda):
    """The isomoprhism between the block containing w_lambda and the
    unique block in #w_lambda.

    """
    h = w_lambda.height
    return Weight(list(reversed(h)), h)

def phi_inv(w_phi_lambda, typ):
    """Invert the isomorphism and recover the weight corresponding to
    phi(lambda) in the block parameterized by typ.

    """
    if not (2 * w_phi_lambda.adeg == len(w_phi_lambda.L) + len(w_phi_lambda.R)):
        raise ValueError(f"{w_phi_lambda} not in unique block of adeg = {w_phi_lambda.adeg}")
    h = w_phi_lambda.R
    atyp = height_to_atyp(h, typ)
    return typ_atyp_to_weight(typ, atyp)

def irrd_char_summands(w_lambda):
    """Return the summands ch(K(mu)) together with their
    multiplicities appearing in the character of L(lambda)

    """
    r = w_lambda.adeg
    typ = w_lambda.typ

    w_phi_lambda = phi(w_lambda)
    A = w_phi_lambda.L

    N = 0
    # iterate thorugh all weakely decreasing integer vectors B bounded
    # above by A
    while N >= 0:
        for delta in IntegerVectors(N, r):
            B = [ai - delta_i for ai, delta_i in zip(A, delta)]
            if not is_weakly_decreasing(B):
                continue
            w_phi_mu = Weight(B, list(reversed(B)))
            w_mu = phi_inv(w_phi_mu, typ)
            yield mult_kac_in_irrd(w_phi_lambda, w_phi_mu), w_mu
                
        N = N + 1
