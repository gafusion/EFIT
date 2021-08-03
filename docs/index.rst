.. efit documentation master file, created by
   sphinx-quickstart on Fri Jul 30 11:16:13 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to EFIT's documentation!
================================


EFIT Equilibrium Fitting code
-----------------------------

EFIT (Equilibrium Fitting) is a computer code developed to translate measurements from plasma diagnostics into useful information like plasma geometry, stored energy, and current profiles. The measuremnts are obtained from diagnostics such as external magnetic probes, external poloidal flux loops, and the Motional Stark Effect (MSE), which measures the direction of the magnetic field lines inside the plasma. The Grad-Shafranov equilibrium equation, which describes the force balance in a plasma, is solved using the available measurements as constraints on the toroidal current density. Since the current also depends on the solution of the equation, the poloidal flux function, this is a nonlinear optimization problem. The equilibrium constraint allows the two-dimensional current density to be represented by two one-dimensional stream functions (functions only of flux), which significantly reduces the complexity of the problem.



Documentation contents
----------------------
.. toctree::
   :maxdepth: 3
   
   subroutines
   modules
   namelist

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
