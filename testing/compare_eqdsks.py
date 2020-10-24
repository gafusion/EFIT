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
    serror_all = 0.
    perror_all = 0.
    for item in items:
        perror = abs((dirs[0][aname][item] - dirs[1][aname][item])/dirs[0][aname][item])
        print(item,'(parall): ', dirs[0][aname][item], dirs[1][aname][item], perror)
        serror = abs((dirs[2][aname][item] - dirs[3][aname][item])/dirs[3][aname][item])
        print(item,'(serial): ', dirs[2][aname][item], dirs[3][aname][item], serror)
        serror_all += serror
        perror_all += perror
    serror /= len(items)
    perror /= len(items)
    print ('\n')
    return serror, perror

def compare_gfiles(time):
    fig = plt.figure()
    gname = 'g' + str(shot).zfill(6) +'.'+ str(time).zfill(5)
    for mydir in dirs:
        mydir[gname].plot()
    fig.savefig(plotdir +gname+'.pdf')

    nw = len(dirs[0][gname]['PSIRZ'])
    nh = len(dirs[0][gname]['PSIRZ'][0,:])
    R = numpy.linspace(dirs[0][gname]['RLEFT'], dirs[0][gname]['RDIM'], nw)              
    Z = numpy.linspace(dirs[0][gname]['ZMID']-dirs[0][gname]['ZDIM'], dirs[0][gname]['ZMID']-dirs[0][gname]['ZDIM'], nw)
    RR, ZZ = numpy.meshgrid(R, Z)
    dpsi1 = numpy.abs((dirs[1][gname]['PSIRZ'] - dirs[0][gname]['PSIRZ'])/dirs[0][gname]['PSIRZ'])
    dpsi2 = numpy.abs((dirs[3][gname]['PSIRZ'] - dirs[2][gname]['PSIRZ'])/dirs[2][gname]['PSIRZ'])
    print(numpy.max(dpsi1),numpy.max(dpsi2))
    fig, axs = plt.subplots(1,2)
    fig.suptitle('Error is PSIRZ')
    axs[0].contourf(RR, ZZ, dpsi1)
    axs[1].contourf(RR, ZZ, dpsi2)

    axs[0].set_xlabel('R[m]')
    axs[0].set_ylabel('Z[m]')
    axs[1].set_xlabel('R[m]')
                       
    fig.savefig(plotdir +'DPSI_' +gname +'.pdf')
    perror = numpy.sum(dpsi1)/nw/nh
    serror = numpy.sum(dpsi2)/nw/nh
    return serror, perror

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
        

print("Compare gEQDSKs")
shot_dir = sys.argv[1]
shot = shot_dir.split('/')[0]
dirs = []
print('./{}/public/parallel'.format(shot_dir))
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

serror_all = 0.
perror_all = 0.
for time in all_times[0]:
    serror, perror = compare_gfiles(time)
    serror_all += serror
    perror_all += perror

serror_all /= len(all_times[0])
perror_all /= len(all_times[0])
print ("Average total difference in PSI (serial) = ", serror_all)
print ("Average total difference in PSI (parallel) = ", perror_all)

serror_all = 0.
perror_all = 0.
for time in all_times[0]:
    serror, perror = compare_afiles(time)
    serror_all += serror
    perror_all += perror

serror_all /= len(all_times[0])
perror_all /= len(all_times[0])
print ("Average total difference in aEQDSK (serial) = ", serror_all)
print ("Average total difference in aEQDSK (parallel) = ", perror_all)

for time in all_times[0]:
    compare_mfiles(time)
