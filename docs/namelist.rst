Input Namelists
===============

.. _namelist:

SETUP
-----

SETUP is the optional namelist that can be read from the efit.input 
that can be used to set table and machine directories and variable
size limits.

.. csv-table:: SETUP
   :file: tables/setup.csv
   :widths: 15,15,70
   :header-rows: 1

OPTIN
-----

OPTIN is the optional namelist that can be read from efit.input 
instead of using the interactive command line inputs.

.. csv-table:: OPTIN
   :file: tables/optin.csv
   :widths: 15,15,70
   :header-rows: 1

IN0
---

IN0 is an optional namelist that can be read from efit.input 
and supercedes the values set in the IN1 namelist.  It only
includes scalar control parameters.  These values are ignored
when running in snap mode (no reason to use this).

.. csv-table:: IN0
   :file: tables/in0.csv
   :widths: 15,15,70
   :header-rows: 1

IN1 General
-----------

IN1 is the main namelist specified in input K, R or X file.
See input namelists for other options.

.. csv-table:: IN1 general
   :file: tables/in1.csv
   :widths: 15,15,70
   :header-rows: 1

IN1 - BASIS
-----------

Plasma toroidal current is modeled by a number of filamentary toroidal current
elements whose `(R,Z)` locations are defined by a set of points on a
(cross-sectional) grid.  Available grid sizes are::

    33 x 33
    33 x 65
    65 x 65
    129 x 129
    257 x 257
    513 x 513
    1025 x 1025
    2049 x 2049

The choice of grid size is determined by the command line parameters given.
For example, `efit 129` uses a `129 x 129` grid.
The amount of current in each filament ("at" each grid point) is represented 
by 

.. math::
    JT(R,\psi)=R[P'(\psi)+\mu _{0}FF'(\psi)/(4\pi^{2}R^{2})]

where :math:`\psi`  is the value of normalized flux at each grid point, 
and :math:`P'(\psi)` and :math:`FF'(\psi)` are represented by 
linear combinations of a small number of basis functions.

The variables below are used to specify the representation of the functions :math:`P'` and
:math:`FF'` and any desired constraints on those representations. All variables are input
via the IN1 or EFITIN namelists.

.. csv-table:: IN1 basis functions
   :file: tables/in1_basis.csv
   :widths: 15,15,70
   :header-rows: 1

IN1 - Fitting
-------------

Many variables used in fitting mode represent data which is extracted from the shot data
base and written into a K-file containing the IN1 and INWANT namelists by running EFIT 
with input mode 5. Values for these variables are typically never entered by the user.

.. csv-table:: IN1 fitting
   :file: tables/in1_fitting.csv
   :widths: 15,15,70
   :header-rows: 1

IN1 - Equilibrium
-----------------

.. csv-table:: IN1 equilibirum
   :file: tables/in1_equilibrium.csv
   :widths: 15,15,70
   :header-rows: 1

INWANT
------

INWANT is specified in input file (K or boundary inputs) for advanced options.

.. csv-table:: INWANT
   :file: tables/inwant.csv
   :widths: 15,15,70
   :header-rows: 1
   
INS
---

INS is specified in input K file for MSE data.

.. csv-table:: INS
   :file: tables/ins.csv
   :widths: 15,15,70
   :header-rows: 1

INVT 
----

INVT is specified in input file (K or boundary inputs) for toroidal rotation.
To turn on toroidal rotation, must set KPRFIT=3 and ICURRT=5 in namelist
IN1. If no pressure data, set NPRESS=0.

.. csv-table:: INVT
   :file: tables/invt.csv
   :widths: 15,15,70
   :header-rows: 1

INK 
---

INK is specified in input K or boundary file for vertical stablization.

.. csv-table:: INK
   :file: tables/ink.csv
   :widths: 15,15,70
   :header-rows: 1

INMS
----

(TODO: need to add/describe namelist)

IN_MSELS
--------

(TODO: need to add/describe namelist)

INA
---

(TODO: need to add/describe namelist)

INLIBIM
-------

(TODO: need to add/describe namelist)

INECE
-----

INECE is specified in input K file for using ece temperature constraints.

.. csv-table:: INECE
   :file: tables/inece.csv
   :widths: 15,15,70
   :header-rows: 1


INER
----

(TODO: need to add/describe namelist)

INSXR
-----

INSXR contains plotting options for the Soft X-Ray diagnostic.
(TODO: need to add/describe namelist)

EDGEP
-----

(TODO: need to add/describe namelist)

EDAT
----

(TODO: need to add/describe namelist)

PROFILE_EXT
-----------

(TODO: need to add/describe namelist)

MACHINEIN (EFUND)
-----------------

MACHINEIN is the namelist specifying machine specific variable sizes for efund and efit
in the mhdin.dat file.

.. csv-table:: MACHINEIN
   :file: tables/machinein.csv
   :widths: 15,15,70

INCHECK
-------

INCHECK is the namelist specifying machine specific limits for efit to use when
IERCHK is set (in the IN0, IN1, or EFITIN namelist).

.. csv-table:: INCHECK
   :file: tables/incheck.csv
   :widths: 15,15,70

EFITIN
------

EFITIN is the main namelist specified by a snap file.

.. csv-table:: EFITIN
   :file: tables/efitin.csv
   :widths: 15,15,70
   :header-rows: 1

EFITINK
-------

EFITINK is for vertical stabilization, the same as INK, but specified in the snap file. It can be included in the input file in file mode.

IN3 (EFUND)
-----------

IN3 is used to read diagnostic parameters from an mhdin.dat file.

.. csv-table:: IN3
   :file: tables/in3.csv
   :widths: 15,15,70
   :header-rows: 1

IN4
---

IN4 is used to read an alternate set of pointnames
from a file.
(TODO: need to add/describe namelist)

IN5 (EFUND)
-----------

IN5 is used to read efund specific parameters from an mhdin.dat file.

.. csv-table:: IN5
   :file: tables/in5.csv
   :widths: 15,15,70
   :header-rows: 1

OUT1
----

OUT1 has part of the results that are written to g-files.  Most variables are repeated from the IN1 or other namelists, but several are updated with computed values.  Unlike the main part of a g-file, this namelist is not designed to be useful externally.

.. csv-table:: OUT1
   :file: tables/out1.csv
   :widths: 15,15,70
   :header-rows: 1


Hardcoded
---------

Unfortunately some variables are not a part of any namelists and can only be manipulated from
within the source code...
Example: see IBOUND

.. csv-table:: hardcoded
   :file: tables/hardcoded.csv
   :widths: 15,15,70
   :header-rows: 1
