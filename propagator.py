import numpy as np
from dirac import Pp, Pm, gamma

# Fermion propagator helper functions
def hat(p):
    return 2. * np.sin(0.5 * p);

def ring(p):
  return np.sin(p);


def funsq(f, p):
    assert isinstance(p, np.ndarray)
    return np.sum(f(p)*f(p))

def A(p, m):
    result = 1. + m
    result += 0.5 * funsq(hat, p)
    return result

def w(p, m):
    arg = A(p, m)
    arg += 1. / arg * (1. + funsq(ring, p))
    return np.arccosh(0.5 * arg)


def Mp(p, m):
  return A(p, m) - np.exp(w(p, m))

def Mm(p, m):
  return A(p, m) - np.exp(-w(p, m))

def R(p, T, m):
    arg = w(p, m) * T
    return Mm(p, m) * np.exp(arg) - Mp(p, m) * np.exp(-arg)

def oneoN(p, T, m):
    #print A(p,m), np.sinh(w(p,m)), R(p,T,m), w(p,m)
    return 2. * A(p, m) * np.sinh(w(p, m)) * R(p, T, m)

def sqrtN(p, T, m):
    #print np.sqrt(1. / oneoN(p, T, m))
    return np.sqrt(1. / oneoN(p, T, m))

def f1(p, t, T, m):
    arg = w(p, m) * t
    if np.linalg.norm(p) < 1e-7 and m == 0:
        return t
        #print sqrtN(p, T, m) * (np.exp(arg) - np.exp(-arg)),p[0],t,T,m,"f1"
    return sqrtN(p, T, m) * (np.exp(arg) - np.exp(-arg))


def f2(p, t, T, m):
    arg = w(p, m) * t
    if np.linalg.norm(p) < 1e-7 and m == 0:
        return 1
        #print sqrtN(p, T, m) * (Mm(p, m) * np.exp(arg) - Mp(p, m) * np.exp(-arg)), p[0],t,T,m,"f2"
    return sqrtN(p, T, m) * (Mm(p, m) * np.exp(arg) - Mp(p, m) * np.exp(-arg))

def Gp(p, x, y, T, m):
  if x > y:
    return f1(p, (T - x), T, m) * f2(p, y, T, m)
  else:
    return f2(p, x, T, m) * f1(p, (T - y), T, m)


def Gm(p, x, y, T, m):
  if x > y:
    return f2(p, (T - x), T, m) * f1(p, y, T, m)
  else:
    return f1(p, x, T, m) * f2(p, (T - y), T, m)


def Hp(p, x, y, T, m):
  if x >= y:
    return f2(p, (T - x), T, m) * f2(p, y, T, m)
  else:
    return Mm(p, m) * Mp(p, m) * f1(p, x, T, m) * f1(p, (T - y), T, m)


def Hm(p, x, y, T, m):
  if x > y:
      return Mm(p, m) * Mp(p, m) * f1(p, (T - x), T, m) * f1(p, y, T, m);
  else:
    return f2(p, x, T, m) * f2(p, (T - y), T, m)


def Sf(p, x, y, T, m = 0.):
    result = Pp * Hp(p, x, y, T, m) + Pm * Hm(p, x, y, T, m) + 1J
    tmp = (Pp * Gp(p, x, y, T, m) + Pm * Gm(p, x, y, T, m))
    for k in range(1,4):
        result -= 1J* gamma[k] * ring(p[k-1]) * tmp
    return result;

