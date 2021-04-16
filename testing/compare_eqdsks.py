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
    serror_imax = (0.,'','')
    perror_imax = (0.,'','')
    for item in items:
        perror = (abs((dirs[0][aname][item] - dirs[1][aname][item])/dirs[0][aname][item]), item, aname)
        print(item,'(parall): ', dirs[0][aname][item], dirs[1][aname][item], perror[0])
        serror = (abs((dirs[2][aname][item] - dirs[3][aname][item])/dirs[3][aname][item]), item, aname)
        print(item,'(serial): ', dirs[2][aname][item], dirs[3][aname][item], serror[0])
        if serror[0] > serror_imax[0]:
            serror_imax = serror
        if perror[0] > perror_imax[0]:
            perror_imax = perror

    print ('\n')
    return serror_max, perror_imax

def compare_gfiles(time):
    fig = plt.figure()
    gname = 'g' + str(shot).zfill(6) +'.'+ str(time).zfill(5)
    for mydir in dirs:
        mydir[gname].plot()
    fig.savefig(plotdir +gname+'.pdf')

    nw = len(dirs[0][gname]['PSIRZ'])
    nh = len(dirs[0][gname]['PSIRZ'][0,:])
    R = numpy.linspace(dirs[0][gname]['RLEFT'], dirs[0][gname]['RLEFT']+dirs[0][gname]['RDIM'], nw)              
    Z = numpy.linspace(dirs[0][gname]['ZMID']-0.5*dirs[0][gname]['ZDIM'], 
                       dirs[0][gname]['ZMID']+0.5*dirs[0][gname]['ZDIM'], nw)
    RR, ZZ = numpy.meshgrid(R, Z)
    dpsi1 = numpy.abs((dirs[1][gname]['PSIRZ'] - dirs[0][gname]['PSIRZ'])/dirs[0][gname]['PSIRZ'])
    dpsi2 = numpy.abs((dirs[3][gname]['PSIRZ'] - dirs[2][gname]['PSIRZ'])/dirs[2][gname]['PSIRZ'])

    perror = (numpy.amax(dpsi1),numpy.unravel_index(dpsi1.argmax(), dpsi1.shape), gname)
    serror = (numpy.amax(dpsi1),numpy.unravel_index(dpsi1.argmax(), dpsi1.shape), gname)

    fig, axs = plt.subplots(1, 2, figsize=(7, 5.))
    fig.suptitle('Difference in PSIRZ (public vs. new)')
    
    cntr0 = axs[0].contourf(RR, ZZ, dpsi1)
    cntr1 = axs[1].contourf(RR, ZZ, dpsi2)
    
    axs[0].axis('equal')
    axs[1].axis('equal')

    axs[0].plot(dirs[0][gname]['RBBBS'], dirs[0][gname]['ZBBBS'])
    axs[0].plot(dirs[1][gname]['RBBBS'], dirs[1][gname]['ZBBBS'])
    axs[1].plot(dirs[2][gname]['RBBBS'], dirs[2][gname]['ZBBBS'])
    axs[1].plot(dirs[3][gname]['RBBBS'], dirs[3][gname]['ZBBBS'])

    axs[0].plot(dirs[0][gname]['RLIM'], dirs[2][gname]['ZLIM'],color='w')
    axs[1].plot(dirs[0][gname]['RLIM'], dirs[3][gname]['ZLIM'],color='w')

    fig.colorbar(cntr0, ax=axs[0])
    fig.colorbar(cntr1, ax=axs[1])
    
    axs[0].set_xlabel('R[m]')
    axs[0].set_ylabel('Z[m]')
    axs[1].set_xlabel('R[m]')
    axs[0].set_title('parallel')
    axs[1].set_title('serial')
                       
    fig.savefig(plotdir+'DPSI_'+gname +'.pdf')
    # Return maximum error, indicies of maximum, and gfile name

    return serror, perror

def compare_mfiles(time):
    fig = plt.figure()
    mname = 'm' + str(shot).zfill(6) +'.'+ str(time).zfill(5)
    for mydir in dirs:
        mydir[mname].plot()
    fig.savefig(plotdir +mname+'.pdf')
    
    fig = plt.figure()
    labels = ['public parallel', 'new parallel','public serial','new serial']
    for i, mydir in enumerate(dirs):
        plt.semilogy(mydir[mname]['cerror']['data'][0,:], label=labels[i])

    plt.ylabel('chi^2')
    plt.xlabel('iteration')
    plt.legend(loc='upper right')
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

serror_max = (0., (0, 0), 'all gEQDSKs')
perror_max = (0., (0, 0), 'all gEQDSKs')
for time in all_times[0]:
    serror, perror = compare_gfiles(time)
    if serror[0] > serror_max[0]:
        serror_max = serror
    if perror[0] > perror_max[0]:
        perror_max = perror
print (serror)
print ("Maximum difference in gEQDSK PSI (serial) = {} at i={}, j={} for : {}".format(
    serror_max[0], serror_max[1][0]+1, serror_max[1][1]+1, serror_max[2]))
print ("Maximum difference in gEQDSK PSI (parallel) = {} at i={}, j={} for : {}".format(
    perror_max[0], perror_max[1][0]+1, perror_max[1][1]+1, perror_max[2]))


serror_max = (0., 'all times', 'all aEQKSKs')
perror_max = (0., 'all times', 'all aEQKSKs')
for time in all_times[0]:
    serror, perror = compare_afiles(time)
    if serror[0] > serror_max[0]:
        serror_max = serror
    if perror[0] > perror_max[0]:
        perror_max = perror


print ("Maximum difference in aEQDSK (serial) = {} found in {} for {}".format(
    serror_max[0], serror_max[1],  serror_max[2]))
print ("Maximum difference in aEQDSK (parallel) = {} found in {} for {}".format(
    perror_max[0], perror_max[1],  perror_max[2]))


for time in all_times[0]:
    compare_mfiles(time)
