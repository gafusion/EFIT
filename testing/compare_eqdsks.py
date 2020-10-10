from omfit.classes.omfit_eqdsk import OMFITgeqdsk, OMFITaeqdsk, OMFITmeqdsk
from omfit.classes.omfit_dir import OMFITdir
import os
import sys
import numpy
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")

def compare_afiles(time):
    aname = 'a' + str(shot).zfill(6) +'.'+ str(time).zfill(5)
    items = ['ali3', 'bcentr','betan', 'cpasma', 'doutl', 'sibdry','simagx','wplasm', 'terror']
    print ('shot: ', shot, 'time ', time)
    for item in items:
        perror = abs((dirs[0][aname][item] - dirs[1][aname][item])/dirs[0][aname][item])
        print(item,'(parall): ', dirs[0][aname][item], dirs[1][aname][item], perror)
        serror = abs((dirs[2][aname][item] - dirs[3][aname][item])/dirs[3][aname][item])
        print(item,'(serial): ', dirs[2][aname][item], dirs[3][aname][item], serror)
    print ('\n')

def compare_gfiles(time):
    fig = plt.figure()
    gname = 'g' + str(shot).zfill(6) +'.'+ str(time).zfill(5)
    for mydir in dirs:
        mydir[gname].plot()
    fig.savefig(plotdir +gname+'.pdf')

    nw = len(dirs[0][gname]['PSIRZ'])
    R = numpy.linspace(dirs[0][gname]['RLEFT'], dirs[0][gname]['RDIM'], nw)              
    Z = numpy.linspace(dirs[0][gname]['ZMID']-dirs[0][gname]['ZDIM'], dirs[0][gname]['ZMID']-dirs[0][gname]['ZDIM'], nw)
    RR, ZZ = numpy.meshgrid(R, Z)
    dpsi1 = numpy.abs((dirs[1][gname]['PSIRZ'] - dirs[0][gname]['PSIRZ'])/dirs[0][gname]['PSIRZ'])
    dpsi2 = numpy.abs((dirs[3][gname]['PSIRZ'] - dirs[2][gname]['PSIRZ'])/dirs[2][gname]['PSIRZ'])

    fig, axs = plt.subplots(1,2)
    fig.suptitle('Error is PSIRZ')
    axs[0].contourf(RR, ZZ, dpsi1)
    axs[1].contourf(RR, ZZ, dpsi2)

    axs[0].set_xlabel('R')
    axs[0].set_ylabel('Z')
    axs[1].set_xlabel('R')
                       
    fig.savefig(plotdir +'DPSI_' +gname +'.pdf')

def compare_mfiles(time):
    fig = plt.figure()
    mname = 'm' + str(shot).zfill(6) +'.'+ str(time).zfill(5)
    for mydir in dirs:
        mydir[mname].plot()
    fig.savefig(plotdir +mname+'.pdf')
    
    fig = plt.figure()
    for mydir in dirs:
        plt.semilogy(mydir[mname]['cerror']['data'][0,:])

    plt.ylabel('chi^2')
    plt.xlabel('iteration')
    fig.savefig(plotdir +"CHI2_" + mname+'.pdf')
        

print("Compare two gEQDSKs")
shot_dir = sys.argv[1]
shot = shot_dir.split('_')[-1]
dirs = []
dirs.append(OMFITdir('./{}/public/parallel'.format(shot_dir)))
dirs.append(OMFITdir('./{}/new/parallel'.format(shot_dir))) 
dirs.append(OMFITdir('./{}/public/serial'.format(shot_dir))) 
dirs.append(OMFITdir('./{}/new/serial'.format(shot_dir)))  

all_times = []
for mydir in dirs:
    times = []
    for item in mydir:
        if item[0] == 'g':
            mydir[item] = OMFITgeqdsk(mydir[item].filename)
            times.append(item.split('.')[-1])
        if item[0] == 'a':
            mydir[item] = OMFITaeqdsk(mydir[item].filename)
        if item[0] == 'm':
            mydir[item] = OMFITmeqdsk(mydir[item].filename)
    all_times.append(times)


plotdir = shot_dir+'/plots/'
if os.path.exists(plotdir):
   os.system("rm -r "+plotdir)
os.mkdir(plotdir)


for time in all_times[0]:
    compare_gfiles(time)
    compare_mfiles(time)
    compare_afiles(time)
