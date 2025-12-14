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

from sage.combinat.partition import Partition
from sage.combinat.partition import Partitions
from sage.combinat.partition import PartitionsInBox
import sage.libs.lrcalc.lrcalc as lrcalc

from .character import mult_kac_in_irrd
from .convertors import parts_to_weight
from .utils import det_factor
from .utils import common_det_factor
from .utils import zero_pad
from .utils import dual
from .weight import Weight

def _dfmt(seq):
    d = len(str(len(seq)))
    fmt = '{:>' + str(d) + 'd}'
    return d, fmt

# --- COMPARISION OF WEIGHTS -------------------------------------------
def _is_gt(w_lambda, w_mu):
    """Check if lambda - mu is a positive linear combination of
    postive roots.

    :returns:
        (-1) if w_lambda < w_mu
        0    if w_lambda = w_mu
        (+1) if w_lambda > w_mu
        None if w_lambda and w_mu are uncomparable

    """
    if w_lambda == w_mu:
        return 0
    
    diff = w_lambda - w_mu
    if diff.is_positive():
        return 1
    
    if (- diff).is_positive():
        return (-1)
    
    return None

def max_weight(weights, comp=_is_gt):
    for w in weights:
        num_smaller = 0
        for v in weights:
            c = comp(w, v)
            if c is None:
                continue
            num_smaller += 1 if c > 0 else 0
        if num_smaller == len(weights) - 1:
            return w
    return None

def _is_lt(a, b):
    gt = _is_gt(a, b)
    if gt is None:
        return None
    return -gt

def min_weight(M):
    return max_weight(M, comp=_is_lt)

def print_comp_matrix(weights, comp=_is_gt, show_weights=False):
    d, fmt = _dfmt(weights)
    symb =  {
         1: ' ' * (d - 1) + '+',
         0: ' ' * (d - 1) + '_' ,
        -1: ' ' * (d - 1) + "-",
        None: ' ' * d,
    }

    for i, w_lambda in enumerate(weights, start=1):
        for w_mu in weights:
            g = comp(w_lambda, w_mu)
            if g == 0:
                print(fmt.format(i), end=' ')
            else:
                print(symb[g], end=' ')    
        if show_weights:
            print(f"\t{w_lambda}")
        else:
            print()
    print()
# --- END COMPARISION OF WEIGHTS ---------------------------------------

# def kacs_containing(w_lambda):
#     """Iterate over all (c, alpha) such that the Kac module K(alpha)
#     has a nonzero lambda weight space with multiplicity c.

#     """
#     k0, lambda_0 = det_factor(w_lambda.coeff_L)
#     k1, lambda_1 = det_factor(w_lambda.coeff_R)

#     for beta in PartitionsInBox(w_lambda.m, w_lambda.n):
#         beta_conj = zero_pad(beta.conjugate(), w_lambda.n)
#         b, beta_conj_dual = det_factor(dual(beta_conj))
#         for alpha_0, c0 in lrcalc.mult(lambda_0, beta).items():
#             if len(alpha_0) > w_lambda.m:
#                 continue
#             for alpha_1, c1 in lrcalc.mult(lambda_1, beta_conj_dual).items():
#                 if len(alpha_1) > w_lambda.n:
#                     continue
#                 a0 = k0
#                 a1 = k1 + b
#                 w_alpha_L = [ai + a0 for ai in zero_pad(alpha_0, w_lambda.m)]
#                 w_alpha_R = [-(ai + a1) for ai in zero_pad(alpha_1, w_lambda.n)]
#                 yield (c0 * c1, Weight(w_alpha_L, w_alpha_R))

def mult_v_in_kac(w_mu, w_alpha):
    """Return the multiplicity of V(mu) in K(alpha)
    """
    eq_dim = (w_mu.m == w_alpha.m) and (w_mu.n == w_alpha.n)
    if not eq_dim:
        return 0

    m, n = w_mu.m, w_mu.n

    # Factor out a common power of the det.
    k0, part_mu_0, part_alpha_0 = common_det_factor(w_mu.coeff_L, w_alpha.coeff_L)
    k1, part_mu_1, part_alpha_1 = common_det_factor(w_mu.coeff_R, w_alpha.coeff_R)


    mult = 0
    
    b = sum(part_mu_1) - sum(part_alpha_1)
    eta_0 = Partition([n + mi for mi in part_mu_0])
    eta_1 = Partition(part_mu_1)
    for beta in Partitions(b, max_length=m, max_part=n, outer=eta_1.conjugate()):
        beta_c = [n - bi for bi in reversed(zero_pad(beta, m))]
        c0 = lrcalc.lrcoef_unsafe(eta_0, part_alpha_0, beta_c)
        
        beta_t = beta.conjugate()
        c1 = lrcalc.lrcoef_unsafe(eta_1, part_alpha_1, beta_t)
        mult += c0 * c1

    return mult

def contains_as_kac_weight(w_mu):
    """Iterate over all alpha such that w_mu is a weight of K(alpha)
    """
    k_mu, part_mu = w_mu.det_factor()
    
    eta_0 = [pi + w_mu.n for pi in part_mu[0]]
    eta_1 = part_mu[1]

    for beta in PartitionsInBox(w_mu.m, w_mu.n):
        beta_c = [w_mu.n - bi for bi in reversed(zero_pad(beta, w_mu.m))]
        beta_t = zero_pad(beta.conjugate(), w_mu.n)

        parts_alpha_0 = lrcalc.skew(eta_0, beta_c)
        for part_alpha_0 in parts_alpha_0.keys():
            if len(part_alpha_0) > w_mu.m:
                continue
                
            parts_alpha_1 = lrcalc.skew(eta_1, beta_t)
            for part_alpha_1 in parts_alpha_1.keys():
                if len(part_alpha_1) > w_mu.n:
                    continue

                alpha_0 = [ai + k_mu[0] for ai in zero_pad(part_alpha_0, w_mu.m)]
                alpha_1 = [ai + k_mu[1] for ai in zero_pad(part_alpha_1, w_mu.n)]
                yield parts_to_weight(alpha_0, alpha_1)

def all_nonzero_kacs(weights, cache=None, verbose=False):
    kac_mults = dict()
    all_kacs = []
    _, fmt = _dfmt(weights)
    for i, w_lambda in enumerate(weights, start=1):
        if verbose:
            print(fmt.format(i) + ": " + "*" * i +  "_" * (len(weights) - i), end="\t")
        
        num_kacs = 0
        for w_alpha in contains_as_kac_weight(w_lambda):
            num_kacs += 1
            for w_mu in weights:
                same_block = w_mu.typ == w_alpha.typ
                if not same_block:
                    continue

                m = mult_kac_in_irrd(w_mu, w_alpha)
                if abs(m) == 0:
                    continue

                _add = False
                kac_mults[w_mu, w_alpha] = m
                if w_alpha not in all_kacs:
                    all_kacs.append(w_alpha)
                    _add = True

                if verbose:
                    print("+" if _add else ".", end='')
        if verbose:
            print(f"\t[{num_kacs}]")

    irrd_mults = dict()
    for w_alpha in all_kacs:
        for w_mu in weights:
            irrd_mults[w_alpha, w_mu] = mult_v_in_kac(w_mu, w_alpha)

    return all_kacs, kac_mults, irrd_mults

def print_mult_matrix(data, row, col):
    d, fmt = _dfmt(col)
    for i in range(len(col)):
        print(fmt.format(i + 1), end=' ')
    print("\n" + '_' * (d + 1) * len(col))

    for i, ri in enumerate(row, start=1):
        for j, cj in enumerate(col, start=1):
            m = data.get((ri, cj), 0)
            print(' ' * (d - 1) + '.' if m == 0 else fmt.format(m), end=' ')
        print(f"| {i}")
    print()

def decompose(module, verbose=False, cache=None):
    """Given weight spaces of a gl(m|n) module M determine a
    decomposition of M in terms of irreducibles L(lambda) in 
    the Grothendieck group. That is given a finite sum
        [M] = sum d_mu [V_mu]
    we determine the coefficients
        [M] = sum e_lambda [L_lambda]
    
    :param module: dict { mu : d_mu }
    :returns:       dict { lambda : e_lambda }
    """
    if verbose:
        print("Given:")
        for i, w_lambda in enumerate(module.keys(), start=1):
            mult = module[w_lambda]
            print(f"{i:4d}: {w_lambda} {'' if mult == 1 else ' x ' + str(mult)}")

        print("\nOrdering:")
        print_comp_matrix(module.keys())

    _data = all_nonzero_kacs(module.keys(), verbose=verbose, cache=cache)
    alphas, kac_in_irrd, weights_in_kac = _data

    if verbose:
        print("\nSupported Kac modules:")
        for i, w_alpha in enumerate(alphas, start=1):
            print(f"{i:4d}: {w_alpha}")

        print("\nMultiplicity in ch(L_lambda) [row] of ch(K_alpha) [col]:")
        print_mult_matrix(kac_in_irrd, module.keys(), alphas)

        print("\nMultiplicity V_mu [row] in K_alpha [col]:")
        print_mult_matrix(weights_in_kac, module.keys(), alphas)

    def _decompose(weights):
        """
        :param weights: dict { weight : mult }
        """
        all_zero = all(m == 0 for w, m in weights.items())
        if len(weights) == 0 or all_zero:
            return dict()

        # if all multplicities are negative, sth has horrible gone wrong(?)
        all_neg = all(m <= 0 for w, m in weights.items())
        if all_neg:
            raise ValueError("all multiplicites are negative!?")

        # determine the larget weight
        w_lambda = max_weight(weights.keys())
        if w_lambda is None:
            w_lambda = next(iter(weights))

        # decompose the largest weight in terms of Kacs
        decomp_lambda = dict()
        for p in filter(lambda p: p[0] == w_lambda, kac_in_irrd.keys()):
            _, w_alpha = p
            b_lambda_alpha = kac_in_irrd[w_lambda, w_alpha]
            for q in filter(lambda q: q[0] == w_alpha, weights_in_kac.keys()):
                _, w_mu = q

                c_alpha_mu = weights_in_kac[w_alpha, w_mu]
                decomp_lambda[w_mu] = decomp_lambda.get(w_mu, 0) + b_lambda_alpha * c_alpha_mu

        # determine the coeffs of [M] - [L_lambda]
        quot_weights = dict()
        for w_mu, d_mu in weights.items():
            r_mu = d_mu - decomp_lambda.get(w_mu, 0)
            if r_mu != 0:
                quot_weights[w_mu] = r_mu

        if verbose:
            print(f"... + {w_lambda}")
            for w_mu, d_mu in quot_weights.items():
                print(f"\t {d_mu:+3d}  {w_mu}")
        rest = _decompose(quot_weights)

        # combine_mults
        mults = { w_lambda : 1 }
        for w_eta, e_eta in rest.items():
            mults[w_eta] = mults.get(w_eta, 0) + e_eta
        return mults

    return _decompose(module)
