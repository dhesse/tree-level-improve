import numpy as np
import sys

mat = np.matrix
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

g5 = gamma[1] * gamma[2] * gamma[3] * gamma[0]

gbar = (gamma[1] + gamma[2] + gamma[3])/np.sqrt(3)

one = mat(np.identity(4))

Pp = 0.5 * (one + gamma[0])
Pm = 0.5 * (one - gamma[0])
