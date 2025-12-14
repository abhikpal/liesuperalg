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

import math

from sage.plot.graphics import Graphics
from sage.plot.text import text
from sage.plot.arc import arc
from sage.plot.plot import list_plot

from .weight import Weight

def cup_diagram(weight, color='blue', **plot_kwargs):
    """Return a plot of the cap diagram of the given weight. 
    """
    diag = Graphics()
    diag.set_aspect_ratio(1)

    # text
    for x in weight.diagram.keys():
        diag += text(f"${x}$", (x, 0.1), vertical_alignment='bottom', color='black')

    cap_starts = sorted([p for p, tag in weight.diagram.items() if tag == 'B'], reverse=True)
    spots_taken = list(weight.diagram.keys())
    for p in cap_starts:
        q = p
        while q in spots_taken:
            q = q + 1
        spots_taken.append(q)
        r = (q - p) / 2
        diag += arc((p + r, 0), r1=r, r2=0.7 * r, sector=(math.pi, 2 * math.pi), axes=False, ymax=1, color=color)

    # markers
    kwargs = {
        'size': 32,
        'color':color
    }
    markers = {
        'N': '^',
        'B': 'v',
        'L': 'x',
        'R': 'o',
    }
    for tag, marker in markers.items():
        points = [(xi, 0) for xi in spots_taken if weight.diagram.get(xi, 'N') == tag]
        diag += list_plot(points, marker=marker, **kwargs, **plot_kwargs)

    diag += list_plot([(min(spots_taken) - 1, 0), (max(spots_taken) + 1, 0)], marker='.', color=color)
    return diag

