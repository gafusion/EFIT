Database
========

This page gives instructions for running accessing the EFIT-AI equilibrium database which 
was generated for the EFIT-AI project, a database of equilibrium data have saved for machine learning
studies. This database consists of X magnetic equilibrium reconstructions, X MSE reconstructions
from the from the 2019 run campaign, and  X kinetic EFIT reconstructions.


Access 
------

The DIII-D EFIT-AI equilibrium database is available to download through `GA sharepoint
<https://fusionga.sharepoint.com/sites/EFIT-AIProject/SitePages/Home.aspx>`_. If you do not 
have a GA sharepoint account and would like access, email ____ for an account.



Data Format
-----------

The EFIT-AI data uses `OMAS <https://gafusion.github.io/omas/schema/schema_equilibrium.html>`_.
to save the data into the IMAS data format as hdf5 files. The IMAS format is a very complete 
record of discharge and EFIT runs. This format can easily converts to IMAS if desired.

If omas and the OMFIT geqdsk libraries are installed, the data can be converted into the 
standard geqdsk format using from_omas()::

	from omfit_classes.omfit_eqdsk import OMFITgeqdsk, OMFITsrc
	from omas import *
	ods = load_omas_h5(filename)
	eq = OMFITgeqdsk('geqdsk').from_omas(ods,time_index=$TIME_INDEX)

Quick-start
-----------

For quick examining the h5 files, the h5py library can be used::

	import glob 
	import h5py
	my_database= "/fusion/projects/theory/mcclenaghanj/magnetic_db_2019/"
	myfiles = glob.glob(my_database+'*.h5')
	myh5 = h5py.File(myfile[0],"r")
	
	
Examining tree structure::

    print(myh5['equilibrium'].keys())
    print(myh5['equilibrium']['time_slice'].keys())
    print(myh5['equilibrium']['time_slice']['0'].keys())
    print(myh5['equilibrium']['time_slice']['0']['global_quantities'].keys())
    
    
which returns::
    
    <KeysViewHDF5 ['code', 'ids_properties', 'time', 'time_slice', 'vacuum_toroidal_field']>
    <KeysViewHDF5 ['0', '1', '10', '100', '101', '102', '103']>
    <KeysViewHDF5 ['boundary', 'constraints', 'global_quantities', 'profiles_1d', 'profiles_2d', 'time']>
    <KeysViewHDF5 ['beta_normal', 'beta_tor', 'ip', 'li_3', 'magnetic_axis', 'psi_axis', 'psi_boundary', 'q_95', 'q_axis', 'q_min']>

	
Plotting a time-trace for a single shot::	
	
    import matplotlib as plt
    import numpy as np
    ip = []
    times = [] 
    betan = []
    for itime in myh5['equilibrium']['time_slice']:
        times.append(np.array(myh5['equilibrium']['time_slice'][itime]['time']['data'])
        ip.append(np.array(myh5['equilibrium']['time_slice'][itime]['global_quantities']['ip']))
        betan.append(np.array(myh5['equilibrium']['time_slice'][itime]['global_quantities']['beta_normal']))
    	
    plot(times, ip*1e-6, '.')
    plot(times, betan*1e-6, '.')
    	
The EFIT-AI team developed a script to extract the necessary data from ...


