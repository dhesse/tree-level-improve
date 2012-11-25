# -*- coding: utf-8 -*-
r"""
:mod:`dirac` -- Minimalistic Dirac algebra.
==================================================

Dirac algebra obeying
:math:`\{\gamma_\mu, \gamma_\nu\} = 2\delta_{\mu\nu}`.
"""

import numpy as np
import sys

mat = np.matrix

#:
#: Dirac matrices :math:`\gamma_{0,1,2,3}`.
gamma = (
    mat([[0,0,1,0],
        [0,0,0,1],
        [1,0,0,0],
        [0,1,0,0]]),
    mat([[0,0,0,-1J],
        [0,0,-1J,0],
        [0,1J,0,0],
        [1J,0,0,0]]),
    mat([[0,0,0,-1],
        [0,0,1,0],
        [0,1,0,0],
        [-1,0,0,0]]),
    mat([[0,0,-1J,0],
        [0,0,0,1J],
        [1J,0,0,0],
        [0,-1J,0,0]])
    )

#: Fifth Dirac matrix :math:`\gamma_5`
g5 = gamma[1] * gamma[2] * gamma[3] * gamma[0]
#: :math:`\bar \gamma = \frac 1 {\sqrt 3} \sum_{i = 1}^3 \gamma_i`
gbar = (gamma[1] + gamma[2] + gamma[3])/np.sqrt(3)
#: unit matrix
one = mat(np.identity(4))
#: :math:`P_+ = \frac 1 2 (1 + \gamma_0)`
Pp = 0.5 * (one + gamma[0])
#: :math:`P_- = \frac 1 2 (1 - \gamma_0)`
Pm = 0.5 * (one - gamma[0])
