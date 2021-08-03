EFIT input namelist
================================

IN1 General
---------------------------------------

IN1 is the main namelist specified in input K, R or X file.
Variables which have no Defaults and thus require inputs are printed in bold. 
See input namelists for other options.

.. csv-table:: IN1 general
   :file: tables/in1.csv
   :widths: 20, 100
   :header-rows: 1

IN1 BASIS
-----------------------------------------

Plasma toroidal current is modeled by a number of filamentary toroidal current elements 
whose (r,z) locations are defined by a set of points on a (cross-sectional) grid. 
Available grid sizes are 33 x 33, 33 x 65, 65 x 65, 129 x 129, 257 x 257, 513 x 513. 
(The choice of grid size is determined by which version of EFIT you run or the command line parameters given.
For example, efitdu3365 uses a 33 x 65 grid, or efitd90 129 129 uses a 129 x 129 grid.)
The amount of current in each filament ("at" each grid point) is represented 
by{\displaystyle JT(R,\psi )=R[P'(\psi )+\mu _{0}FF'(\psi )/(4\pi ^{2}R^{2})]} 
where {\displaystyle \psi } is the value of normalized flux at each grid point, 
and {\displaystyle P'(\psi )} and {\displaystyle FF'(\psi )} are represented by 
linear combinations of a small number of basis functions.

The variables below are used to specify the representation of the functions P' and FF' and 
any desired constraints on those representations. All variables are input via the IN1 or 
EFITIN namelists. Variables which have no defaults (and thus require inputs??) are printed 
in bold.

.. csv-table:: IN1 basis functions
   :file: tables/in1_basis.csv
   :widths: 20, 100
   :header-rows: 1

IN1 Fitting
----------------------------------------------------

Many variables used in fitting mode represent data which is extracted from the shot data
base and written into a K-file containing the IN1 and INWANT namelists by running EFIT 
with input mode 5. Values for these variables are typically never entered by the user.


.. csv-table:: IN1 fitting
   :file: tables/in1_fitting.csv
   :widths: 20, 100
   :header-rows: 1

IN1 - Equilibrium
--------------------------------------------------------

Variables which have no Defaults (and thus require inputs??) are printed in bold.

.. csv-table:: IN1 equilibirum
   :file: tables/in1_equilibrium.csv
   :widths: 20, 100
   :header-rows: 1

INWANT
------------------------------------------

INWANT is specified in input file (K or boundary inputs) for advanced options. Variables which have no defaults (and thus require inputs??) are printed in bold.

.. csv-table:: INWANT
   :file: tables/inwant.csv
   :widths: 20, 100
   :header-rows: 1


INVT 
----------------------------------------

INVT is specified in input file (K or boundary inputs) for toroidal rotation. Variables which have no defaults (and thus require inputs??) are printed in bold.

To turn on toroidal rotation, must set KPRFIT=3 and ICURRT=5 in namelist IN1. If no pressure data, set NPRESS=0.


.. csv-table:: INVT
   :file: tables/invt.csv
   :widths: 20, 100
   :header-rows: 1