

EFIT Equilibrium Fitting code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EFIT (Equilibrium Fitting) is a computer code developed to translate measurements from plasma diagnostics into useful information like plasma geometry, stored energy, and current profiles. The measuremnts are obtained from diagnostics such as external magnetic probes, external poloidal flux loops, and the Motional Stark Effect (MSE), which measures the direction of the magnetic field lines inside the plasma. The Grad-Shafranov equilibrium equation, which describes the force balance in a plasma, is solved using the available measurements as constraints on the toroidal current density. Since the current also depends on the solution of the equation, the poloidal flux function, this is a nonlinear optimization problem. The equilibrium constraint allows the two-dimensional current density to be represented by two one-dimensional stream functions (functions only of flux), which significantly reduces the complexity of the problem.

Key references are: :cite:`lao85,lao05`.


User documentation
~~~~~~~~~~~~~~~~~~
.. toctree::
   :maxdepth: 2
   
   install
   quick-start
   namelist
   dictionary
   license

Developer documentation
~~~~~~~~~~~~~~~~~~~~~~~~
.. toctree::
   :maxdepth: 2
   
   developer
   docker
   subroutines
   modules
   
EFIT-AI documentation
~~~~~~~~~~~~~~~~~~~~~~~~

Here are topics that are not formally related to EFIT itself, but rather are
part of the larger EFIT-AI project.  

.. toctree::
   :maxdepth: 2
   
   database

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Bibliography
------------
.. bibliography:: /efit.bib
