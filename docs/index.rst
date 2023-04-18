EFIT Equilibrium Fitting code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EFIT (Equilibrium Fitting) is a computer code developed to translate measurements from plasma diagnostics into useful information like plasma geometry, stored energy, and current profiles. The measuremnts are obtained from diagnostics such as external magnetic probes, external poloidal flux loops, and the Motional Stark Effect (MSE), which measures the direction of the magnetic field lines inside the plasma. The Grad-Shafranov equilibrium equation, which describes the force balance in a plasma, is solved using the available measurements as constraints on the toroidal current density. Since the current also depends on the solution of the equation, the poloidal flux function, this is a nonlinear optimization problem. The equilibrium constraint allows the two-dimensional current density to be represented by two one-dimensional stream functions (functions only of flux), which significantly reduces the complexity of the problem.

The version of the code described here is machine independent, portable, and has better performance than the original GA version.  It is also being designed with the goal of connecting to other components of the `EFIT-AI project <https://fusion.gat.com/conference/e/EFIT-AI>`__.

Key references are: :cite:`lao85,lao05,lao22`.


EFIT-AI documentation
~~~~~~~~~~~~~~~~~~~~~~~~

EFIT is part of the larger `EFIT-AI project <https://fusion.gat.com/conference/e/EFIT-AI>`__.

Documentation on other project can be found at:

  + `EFIT-AI database <https://efit-ai.gitlab.io/efitai_database>`__.
  + `UNBAFFELD <https://efit-ai.gitlab.io/unbaffeld>`__.
  + `EFIT-MORNN <https://efit-ai.gitlab.io/efit_mornn>`__.


User documentation
~~~~~~~~~~~~~~~~~~
.. toctree::
   :maxdepth: 2
   
   install
   quick-start
   namelist
   files
   dictionary
   license

Developer documentation
~~~~~~~~~~~~~~~~~~~~~~~~
.. toctree::
   :maxdepth: 2
   
   developer
   gen_subroutines
   gen_functions
   gen_modules
   
EFIT-AI documentation
~~~~~~~~~~~~~~~~~~~~~~~~

EFIT is part of the `EFIT-AI project <https://fusion.gat.com/conference/e/EFIT-AI>`__.

Related documentation from this project:

  + `EFIT-AI database <https://efit-ai.gitlab.io/efitai_database>`__.
  + `unbaffeld <https://efit-ai.gitlab.io/unbaffeld>`__.

.. Indices and tables

Search
==================

* :ref:`search`
* :ref:`genindex`

.. TODO: The following pages still do not render properly online
.. * :ref:`modindex`


Bibliography
------------
.. bibliography:: /efit.bib
