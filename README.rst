==================================
 Tree level improvement in the SF
==================================

These scripts calculate tree level improvement of observables as
described by Divitiis et al. (**ALPHA** collaboration) in
``[hep-lat/9411017]``. The observables in observables.py are named as
in ``[arXiv:1001.4783]`` (Blossier et al. for the **ALPHA**
collaboration) and the references therein.

**DO NOT FORGET TO CITE THE SOURCES ABOVE IN PUBLICATIONS OF ANY FORM**


Files
=====

``calculate_tly.py``
    the main script

``observables.py``
    store your observables here

``tli.py``
    contains the actual tli function

``propagator.py``
    tree level SF quark propagator

``dirac.py``
    Dirac algebra

Usage
=====

Further observables should be defined in observables.py.
Their tli can be calculated by adding them to the arguments dict
contained int calculate_tli.py. Type::

  python tli.py # to cross check with internal notes by M. della Morte
  python calculate_tli.py # to calculate tli of observables

Unit tests
==========

There are unit tests included in the ``test`` directory.