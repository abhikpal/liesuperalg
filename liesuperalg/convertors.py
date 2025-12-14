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

from .weight import Weight
from .weight import rho

def parts_to_weight(part_L, part_R):
    """Convert two partitions to weight"""
    return Weight(part_L, [-ri for ri in part_R])

def typ_atyp_to_diagram(typ, atyp):
    """Returns a diagram corresponding to a typ and atyp tuples.

    :returns: diagram
    :rtype: dict
    """
    diag = dict()
    for li in typ[0]:
        diag[li] = 'L'
    for ri in typ[1]:
        diag[ri] = 'R'
    for ai in atyp:
        diag[ai] = 'B'
    return diag

def diagram_to_weight(diagram):
    # the weight returned is dominant!
    _L = []
    _R = []
    for num, label in diagram.items():
        if label in ['B', 'L']:
            _L.append(num)
        if label in ['B', 'R']:
            _R.append(num)
    L_rho = list(sorted(_L, reverse=True))
    R_rho = list(sorted(_R))
    weight_rho = Weight(L_rho, R_rho)
    return weight_rho - rho(weight_rho.m, weight_rho.n)

def typ_atyp_to_weight(typ, atyp):
    # the weight returned is dominant!
    return diagram_to_weight(typ_atyp_to_diagram(typ, atyp))

def height_to_atyp(h, typ):
    """
    Given a height vector and typical tuple detemrine the atypical tuple.
    """
    a = [hs + s for s, hs in enumerate(h, start=1)]
    for xi in sorted(typ[0] + typ[1]):
        for t, at in enumerate(a):
            if at >= xi:
                a[t] = at + 1
    return a

def str_to_weight(w_str):
    str_to_int_list = lambda s: [int(i) for i in s.split(',')]
    parts = w_str[1:-1].split('|')
    L = str_to_int_list(parts[0])
    R = str_to_int_list(parts[1])
    return Weight(L, R)
