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

Building Input (feature) and Output spaces for EFIT-MOR
-------------------------------------------------------

The EFIT-AI team developed a script to extract the necessary data from the OMAS database (DB). 
The said script, called NN_io_fromMagDB.py, contains a series of modules that operate on either
each time slice or discharge file from the DB (i.e. the .h5 files contained therein). 

For any discharge contained in the 2019 EFIT01 or EFIT02 DB, the IO vector for the EFIT-MOR 
can be built via the following Pyhton commands::

     from NN_io_fromMagDB import extract_inputs, extract_outputs
     hin = h5py.File(input_file, mode="r")
     inputs, times = extract_inputs(input_file, '01')
     EFIToutputs = extract_outputs(input_file)

where the flag 01 in extract_inputs() indicates that data is being pulled from EFIT01 (magnetics 
only). The user must switch this flag to 02 to pull the MSE data along with the magnetics from
EFIT02. This way of doing things will be phased out eventually. 

EFIToutputs a list that contains 4 3D arrays: the normalized poloidal flux, Btor, 
JtorDS ( curl curl of psi) and JtorFB = R*p' + FF'/(mu0*R) on the RZ grid as a function of the 
time slices contained within the discharge.

The inputs as a function of time and one of the 4 outputs as a function of the RZ grid can be
plotted via::

    plot_data(inputs, EFIToutputs, times, R.T, Z.T)

where the RZ grid coords can be obtained with::

    R, Z  = get_grid(time_slice)

after invoking from NN_io_fromMagDB import get_grid(). 
The time-slice for the 2D contour plot of the EFIToutput (be it psi or Btor or Jtor etc) is 
chosen randomly over the available time slices encoded in the 1D vector times. 

Time slices with BetaN>10.0 are automatically discarded as bad data inside
the functions extract_inputs() amd extract_outputs()

Building the IO vectors for a series of discharges
++++++++++++++++++++++++++++++++++++++++++++++++++

An additional parallel Python script, save_efitmor_inputs.py, has been created to loop over a 
series of discharge files and monotonically build an input and output vector for the ML training. 
Because a small subset of the discharges are either a vacuum or a test shot, the user must invoke
the function find_goodShots(discharge_files), which checks the list of .h5 files, discharge_files, for
both vacuum and test shots, and only keeps non-vacuum and non-test shots. 

.. Warning::
     The user has to be sourcing the text file GoodShots.txt or GoodShots2019.txt 
     for the function find_goodShots() to sift through the discharges correctly.  
    
The code snippet to prune out the:: 

    # determine the "good" shots
    num_good_shots, goodID = findGoodShots(discharge_files)

The main routine is as follows::
    
    # // process the data 
    numthreads = mp.cpu_count()
    pool = mp.Pool(processes=numthreads)
    results = pool.map( pull_data, [infile for infile in discharge_files[goodID]])
    pool.close()
    pool.join()

    # parse the results into the separate arrays to create the IO of ML
    for ii in np.arange(num_good_shots):
        if ii == 0:
            inputs = results[ii][0]
            times = results[ii][1]
            psi = results[ii][2]
            Jtor = results[ii][3]
            JtorDS = results[ii][4]
        else:
            inputs = np.append(inputs, results[ii][0], axis=1)  
            times = np.append(times, results[ii][1])  
            psi = np.append(psi, results[ii][2], axis=2)  
            Jtor = np.append(Jtor, results[ii][3], axis=2)  
            JtorDS = np.append(JtorDS, results[ii][4], axis=2)  
 
