EFIT subroutines
================================

The logic flow chart::

	EFIT: Main driver
		GETSETS: Initialize run, input snap file, Green's tables.
		LOOP: Loop over timeslices KTIME.
			DATA: Read input file or fetch data, set up fitting weight.
			INICUR: Initialize plasma current distribution.
			FIT: Equilibrium and fitting iterations.
				LOOP: Loop over MXITER, current profile.
					GREEN: Update respond matrix.
					PRESUR: Update pressure profile data.
					MATRIX: Update current and pressure profile.
					LOOP: Loop over equilibrium iteration NXITER.
						CURRNT: Compute plasma current.
						FCURRT: Get external shaping coil currents.
						STEPS: Find magnetic axis and trace boundary.
						FINDAX: Set up bicubic spline.
						BOUND: Get boundary.
						FINDAX: Get magnetic axis.
						RESIDU: Check equilibrium convergence.
			SHAPE: Compute plasma parameters, graphics after done.
			PLTOUT: Generate graphics.
			PRTOUT: Print out.
			WEQDSK: Write G EQDSK.
			SHIPIT: Write A EQDSK.


modules-efit.F90
----------------

.. doxygenfunction:: fch5init
.. doxygenfunction:: fch5resetvars
.. doxygenfunction:: errctrl_setstate
.. doxygenfunction:: errctrl_msg

set_extvars.F90
---------------

.. doxygenfunction:: set_extvars

read_namelist.F90
-----------------

.. doxygenfunction:: read_optin
.. doxygenfunction:: read_eparm
.. doxygenfunction:: read_dirs_shot
.. doxygenfunction:: read_dirs_shot_omas

get_opt_input.F90
-----------------

.. doxygenfunction:: get_opt_input

set_eparm.f90
-------------

.. doxygenfunction:: set_eparm_defaults
.. doxygenfunction:: set_eparm_dependents

tables.F90
----------

.. doxygenfunction:: table_name_ch
.. doxygenfunction:: set_table_dir
.. doxygenfunction:: read_tables

getsets.F90
-----------

.. doxygenfunction:: getsets

set_defaults.f90
----------------

.. doxygenfunction:: set_defaults

set_basis.F90
-------------

.. doxygenfunction:: set_basis_params
.. doxygenfunction:: setff
.. doxygenfunction:: setfp
.. doxygenfunction:: setfpp
.. doxygenfunction:: setpp
.. doxygenfunction:: setppp
.. doxygenfunction:: setpr
.. doxygenfunction:: setpw
.. doxygenfunction:: setpwp
.. doxygenfunction:: setpwpp
.. doxygenfunction:: seter
.. doxygenfunction:: seterp

ppbasisfunc.f90
---------------

.. doxygenfunction:: ppcnst
.. doxygenfunction:: ppstore

ffbasisfunc.f90
---------------

.. doxygenfunction:: ffcnst
.. doxygenfunction:: ffstore

wwbasisfunc.f90
---------------

.. doxygenfunction:: wwcnst
.. doxygenfunction:: wwstore

eebasisfunc.f90
---------------

.. doxygenfunction:: eecnst
.. doxygenfunction:: eestore

spline.f90
----------

.. doxygenfunction:: seva2d
.. doxygenfunction:: sets2d
.. doxygenfunction:: spl2bc
.. doxygenfunction:: spl2pp
.. doxygenfunction:: eknot
.. doxygenfunction:: spli2d
.. doxygenfunction:: bsplvb
.. doxygenfunction:: banslv
.. doxygenfunction:: banfac
.. doxygenfunction:: interv
.. doxygenfunction:: linv1f

zpline.f90
----------

.. doxygenfunction:: zpline
.. doxygenfunction:: spleen
.. doxygenfunction:: splaan

get_exp_data.f90
----------------

.. doxygenfunction:: getlim
.. doxygenfunction:: getsxr

constraints.F90
---------------

.. doxygenfunction:: get_constraints
.. doxygenfunction:: avdata
.. doxygenfunction:: amdata
.. doxygenfunction:: apdata
.. doxygenfunction:: gettanh
.. doxygenfunction:: avdiam
.. doxygenfunction:: zmooth
.. doxygenfunction:: smoothit
.. doxygenfunction:: smoothit2
.. doxygenfunction:: zplines
.. doxygenfunction:: magsigma
.. doxygenfunction:: get_constraints_mpi

getdat.F90
----------

.. doxygenfunction:: getdat

diamagnetic.f90
---------------

.. doxygenfunction:: getdia
.. doxygenfunction:: dlcomp
.. doxygenfunction:: lowpass
.. doxygenfunction:: interp

stark.F90
---------

.. doxygenfunction:: getstark
.. doxygenfunction:: getstark_mpi
.. doxygenfunction:: setstark
.. doxygenfunction:: fixstark

msels.f90
---------

.. doxygenfunction:: getmsels
.. doxygenfunction:: msels_data
.. doxygenfunction:: msels_hist

data_input.f90
--------------

.. doxygenfunction:: data_input

set_filenames.f90
-----------------

.. doxygenfunction:: setfnmeq
.. doxygenfunction:: setfnmd
.. doxygenfunction:: setfnmt
.. doxygenfunction:: setfnmpl
.. doxygenfunction:: setfnmq

auto_knot.F90
------------

.. doxygenfunction:: autoknot
.. doxygenfunction:: restore_autoknotvals
.. doxygenfunction:: store_autoknotvals

inicur.F90
----------

.. doxygenfunction:: inicur

fit.F90
-------

.. doxygenfunction:: fit
.. doxygenfunction:: residu

ece.F90
-------

.. doxygenfunction:: setece
.. doxygenfunction:: geteceb
.. doxygenfunction:: getecer
.. doxygenfunction:: gettir

green.F90
---------

.. doxygenfunction:: green

matrix.f90
----------

.. doxygenfunction:: matrix

current.f90
-----------

.. doxygenfunction:: currnt

external_current.f90
--------------------

.. doxygenfunction:: fcurrt
.. doxygenfunction:: vescur

pflux.90
--------

.. doxygenfunction:: pflux

buneman.f90
-----------

.. doxygenfunction:: buneto
.. doxygenfunction:: rzpois

cyclic.F90
----------

.. doxygenfunction:: cyclic_reduction
.. doxygenfunction:: pflux_cycred
.. doxygenfunction:: vsma_
.. doxygenfunction:: ef_vvmul
.. doxygenfunction:: ef_tridiag2
.. doxygenfunction:: ef_tridiag1
.. doxygenfunction:: ef_vadd_shrt
.. doxygenfunction:: ef_vmul_const_shrt

update_parameters.90
--------------------

.. doxygenfunction:: update_params
.. doxygenfunction:: weight
.. doxygenfunction:: chisqr

boundary.f90
------------

.. doxygenfunction:: bound
.. doxygenfunction:: cellb
.. doxygenfunction:: chkcrn
.. doxygenfunction:: cntour
.. doxygenfunction:: extrap
.. doxygenfunction:: findax
.. doxygenfunction:: fqlin
.. doxygenfunction:: maxpsi
.. doxygenfunction:: minmax
.. doxygenfunction:: order
.. doxygenfunction:: packps
.. doxygenfunction:: qfit
.. doxygenfunction:: surfac
.. doxygenfunction:: zlim

pressure.F90
------------

.. doxygenfunction:: presur
.. doxygenfunction:: presurw

get_kinetic_data.f90
--------------------

.. doxygenfunction:: getne
.. doxygenfunction:: getbeam
.. doxygenfunction:: gette
.. doxygenfunction:: gettion

shapesurf.f90
-------------

.. doxygenfunction:: shapesurf
.. doxygenfunction:: dslant

beta_li.F90
-----------

.. doxygenfunction:: betali
.. doxygenfunction:: betsli

utils.F90
---------

.. doxygenfunction:: fluxav
.. doxygenfunction:: splitc
.. doxygenfunction:: tsorder
.. doxygenfunction:: fitpp
.. doxygenfunction:: fitfp
.. doxygenfunction:: lenco2

mat_solve.f90
--------------

.. doxygenfunction:: decomp
.. doxygenfunction:: solve
.. doxygenfunction:: sdecm

weq.f90
-------

.. doxygenfunction:: shipit
.. doxygenfunction:: weqdsk
.. doxygenfunction:: timdot

wmeasure.F90
------------

.. doxygenfunction:: wmeasure

write_K.F90
-----------

.. doxygenfunction:: write_K
.. doxygenfunction:: write_K2

wtime.F90
---------

.. doxygenfunction:: wtime

prtout.F90
----------

.. doxygenfunction:: prtout
.. doxygenfunction:: prtoutheader
