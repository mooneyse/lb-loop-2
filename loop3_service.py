# loop3_service   Neal Jackson 1.10.2018
# Service routines for loop3.
# Response to complaints that this looks too much like AIPS will involve adding
#    APARM arrays as arguments.

import numpy as np,os,sys,glob,time,losoto.h5parm,scipy,pickle
from scipy import interpolate
import pyrap.tables as pt
import matplotlib; from matplotlib import pyplot as plt
import h5py

# Given two calibration structures H1 and H2, and antennas to interpolate, replace the
# phase calibration of H1 for each antenna with an interpolated version of H2.

def clcal (H1,H2,ant_interp=None,dozero=False):
    h1,h2 = h5py.File(H1,'r+'),h5py.File(H2)
    n1,n2 = h1.get('sol000/phase000'),h2.get('sol000/phase000')
    t1,t2 = np.array(n1['time']),np.array(n2['time'])
    v1,v2 = np.array(n1['val']),np.array(n2['val'])
    a1 = np.array(h1.get('sol000/phase000/ant'))
    ant_interp = a1 if ant_interp==None else ant_interp
    for i in range(len(a1)):
        if a1[i] not in ant_interp:
            continue
        for iz in range(v1.shape[1]):
            for ipol in range(v1.shape[3]):
                if dozero:
                    z = np.zeros_like (v1[0])
                else:
                    z = scipy.interpolate.griddata(t2,np.unwrap(v2[:,iz,i,ipol]),t1,method='linear')
                    while z.max()>np.pi:
                        np.putmask(z,z>np.pi,z-2.*np.pi)
                    while z.min()<-np.pi:
                        np.putmask(z,z<-np.pi,z+2.*np.pi)
                h1['sol000/phase000/val'][:,iz,i,ipol] = z
    h1.close(); h2.close()
    h1 = h5py.File(H1,'r+')
    n1 = h1.get('sol000/phase000')
    v1 = np.array(n1['val'])
    h1.close()

def calib (vis,incol='DATA',outcol='DATA',solint=180,solmode='P',\
           model=None,outms='.',outcal=None):
    outcal = vis+'_cal' if outcal==None else outcal
    mgain = 'sourcedb=%s\n'%model if model else 'usemodelcolumn=true\n'
    caltype = 'phaseonly' if solmode=='P' else 'diagonal'
    f=open('calib.parset','w')
    f.write('msin=%s\n'%vis)
    f.write('msin.datacolumn=%s\n'%incol)
    f.write('msout=%s\n'%outms)
    f.write('msout.datacolumn=%s\n'%outcol)
    f.write('steps=[gaincal]\n')
    f.write('gaincal.'+mgain)
    f.write('gaincal.caltype=%s\n'%caltype)
    f.write('gaincal.solint=%i\n'%(solint))
    f.write('gaincal.usebeammodel=False\n')
    f.write('gaincal.parmdb=%s\n'%outcal)
    f.write('gaincal.applysolution=%s\n'%('False' if incol==outcol else 'True'))
    f.close()
    time_start = time.time()
    # Bug fix here: NDPPP leaves the .h5 files unclosed. So we have to start a separate python
    # session to run the NDPPP on calib.parset, which closes the .h5 files on exit.
    fo=open('calib.py','w')
    fo.write ('import os\nos.system(\'NDPPP calib.parset\')\n')
    fo.close()
    os.system('python calib.py')
    #    os.system('NDPPP calib.parset')
    time_end = time.time()
    print 'NDPPP took %d s' % int(time_end-time_start)


# Make the coherence parameter. This relies on the difference in the phase
# solutions in XX and YY remaining constant if the solutions are coherent.
# Also need to return an incoherent answer (2.0) if there are too many NaN
# solutions (here >10%)
def coherence_metric (htab='1327_test.ms_cal.h5',solset='sol000',soltab='phase000'):
    NANFRAC, INCOH = 0.1, 2.0
    # do this in another python instance otherwise the hparm does not close properly
    # write out the dictionaries/arrays, then load again from pickle. Horrible and
    # should be replaced once anyone can be found who knows how to close h5 files.
    fo = open('tmp_cohmetric.py','w')
    fo.write ('import h5parm,os,pickle\n')
    fo.write ('tab = h5parm.openSoltab(\'%s\',solsetName=\'%s\',soltabName=\'%s\')\n' % \
                          (htab,solset,soltab)   )
    fo.write ('v, vm = tab.getValues()[0], tab.getValues()[1]\n')
    fo.write ('pickle.dump(v,open(\'v.pkl\',\'wb\'))\n')
    fo.write ('pickle.dump(vm,open(\'vm.pkl\',\'wb\'))\n')
    fo.close()
    os.system ('python tmp_cohmetric.py')
    v = pickle.load (open('v.pkl','rb'))
    vm = pickle.load (open('vm.pkl','rb'))
    ant,freq,pol,time = vm['ant'],vm['freq'],vm['pol'],vm['time']
    coh = np.array([])
    for i in range(len(ant)):   # assumes two polarizations XX YY
    #     changed this (njj) - note that np.unwrap gives an array full of NaN
    #     if even the first element of the input array is NaN
    #        diff = np.unwrap(v[:,0,i,0]-v[:,0,i,1])
        diff = v[:,0,i,0]-v[:,0,i,1]
        if float(len(diff[np.isnan(diff)]))>NANFRAC*float(len(diff)):
            coh = np.append(coh,INCOH)
        else:
            diff = np.unwrap(diff[~np.isnan(diff)])
            coh = np.append(coh,np.nanmean(np.gradient(abs(diff))**2))
    return coh


def snplt (htab='1327_test.ms_cal.h5',solset='sol000',soltab='phase000',antenna=None,nplot=6,
           outpng=None):
    outpng = outpng if outpng else htab
    fo = open('tmp_snplt.py','w')
    fo.write ('import h5parm,os,pickle\n')
    fo.write ('tab = h5parm.openSoltab(\'%s\',solsetName=\'%s\',soltabName=\'%s\')\n' % \
                          (htab,solset,soltab)   )
    fo.write ('v, vm = tab.getValues()[0], tab.getValues()[1]\n')
    fo.write ('pickle.dump(v,open(\'v.pkl\',\'wb\'))\n')
    fo.write ('pickle.dump(vm,open(\'vm.pkl\',\'wb\'))\n')
    fo.close()
    os.system ('python tmp_snplt.py')
    v = pickle.load (open('v.pkl','rb'))
    vm = pickle.load (open('vm.pkl','rb'))
    #    tab = h5parm.openSoltab(htab,solsetName=solset,soltabName=soltab)
    #    v, vm = tab.getValues()[0], tab.getValues()[1]
    ant,freq,pol,time = vm['ant'],vm['freq'],vm['pol'],vm['time']
    time = 24.*(time/86400. - int(time[0])/86400)
    iplot = 0
    antenna = antenna if antenna else ant
    while iplot<len(antenna):
        a = antenna[iplot]
        aidx = np.argwhere(ant==a)[0][0]
        sys.stdout.write(a+' ')
        for ipol in range(v.shape[3]):
            if not (iplot+1)%nplot:
                plt.subplot(nplot,1,1+iplot%nplot)
            else:
                plt.subplot(nplot,1,1+iplot%nplot,xticks=[])
            plt.plot(time,np.rad2deg(v[:,0,aidx,ipol]),'+')
            plt.ylim(-180.,180.);plt.xlim(time[0],time[-1])
            plt.text(time[0],180.0-12.*nplot,a)
            plt.subplots_adjust(wspace=0,hspace=0)
        iplot+=1
        if not iplot%nplot:
            os.system('rm '+outpng+'_%d.png'%(iplot//nplot -1))
            print '-> '+outpng+'_%d.png'%(iplot//nplot -1)
            plt.savefig(outpng+'_%d.png'%(iplot//nplot -1),bbox_inches='tight')
            plt.clf()
    if iplot%nplot:
        os.system('rm '+outpng+'_%d.png'%(iplot//nplot))
        print '-> '+outpng+'_%d.png'%(iplot//nplot)
        plt.savefig(outpng+'_%d.png'%(iplot//nplot),bbox_inches='tight')


# Because I don't like writing enormous command lines in code. Also only have to change once if the
# wsclean arguments change - or indeed if we use a different imager.
def imagr (vis,threads=0,mem=100,doupdatemodel=True,tempdir='',dosaveweights=False,doprimary=False,\
           robust=0,domfsweight=False,gausstaper=0.0,tukeytaper=0.0,dostoreweights=False,outname='wsclean',\
           imsize=1024,cellsize='0.1asec',dopredict=False,niter=10000,pol='I',datacolumn='',autothreshold=0.,\
           dolocalrms=False,gain=0.1,mgain=1.0,domultiscale=False,dojoinchannels=False,channelsout=1,fitsmask='',\
           baselineaveraging=0.0,maxuvwm=0.0,minuvwm=0.0,maxuvl=0.0,minuvl=0.0,dostopnegative=False,automask=0.,\
           dosavesourcelist=False,weightingrankfilter=0.0,weightingrankfiltersize=0.0):
    cmd = 'wsclean '
    cmd += ('' if not threads else '-j '+str(threads)+' ')
    cmd += ('' if not mem==100 else '-mem '+str(mem)+' ')
    cmd += ('' if doupdatemodel else '-no-update-model-required ')
    cmd += tempdir+' '
    cmd += ('' if not dosaveweights else '-save-weights ')
    cmd += ('' if not doprimary else '-apply-primary-beam ')
    if robust >=5:
        cmd += '-weight natural '
    elif robust <=-5:
        cmd += '-weight uniform '
    else:
        cmd += '-weight briggs %f '%robust
    cmd += ('' if not domfsweight else '-mfs-weighting ')
    cmd += ('' if gausstaper==0.0 else '-taper-gaussian %f '%gausstaper)
    cmd += ('' if tukeytaper==0.0 else '-taper-tukey %f '%tukeytaper)
    cmd += ('' if not dostoreweights else '-store-imaging-weights ')
    cmd += '-name '+outname+' '
    cmd += '-size '+str(imsize)+' '+str(imsize)+' '
    cmd += '-scale '+str(cellsize)+' '
    cmd += ('' if not dopredict else '-predict ')
    cmd += ('-niter '+str(niter)+' ')
    cmd += ('' if pol=='I' else '-pol '+pol+' ')
    cmd += ('' if datacolumn=='' else '-datacolumn %s '%datacolumn)
    cmd += ('' if autothreshold==0. else '-auto-threshold %f '%autothreshold)
    cmd += ('' if not dolocalrms else '-local-rms ')
    cmd += ('' if not domultiscale else '-multiscale ')
    cmd += ('' if not dojoinchannels else '-join-channels ')
    cmd += ('-channels-out %d '%channelsout if channelsout!=1 else '')
    cmd += ('' if gain==0.1 else '-gain %f '%gain)
    cmd += ('' if mgain==1.0 else '-mgain %f '%mgain)
    cmd += ('' if fitsmask=='' else '-fits-mask %s '%fitsmask)
    cmd += ('' if baselineaveraging==0.0 else '-baseline-averaging %f '%baselineaveraging)
    cmd += ('' if maxuvl==0.0 else '-maxuv-l %f '%maxuvl)
    cmd += ('' if minuvl==0.0 else '-minuv-l %f '%minuvl)
    cmd += ('' if maxuvwm==0.0 else '-maxuvw-m %f '%maxuvwm)
    cmd += ('' if minuvwm==0.0 else '-minuvw-m %f '%minuvwm)
    cmd += ('' if not dostopnegative else '-stop-negative ')
    cmd += ('' if automask==0. else '-auto-mask %f '%automask)
    cmd += ('' if not dosavesourcelist else '-save-source-list ')
    cmd += ('' if weightingrankfilter==0.0 else '-weighting-rank-filter %f '%weightingrankfilter)
    cmd += ('' if weightingrankfiltersize==0.0 else '-weighting-rank-filter-size %f '%weightingrankfiltersize)
    cmd += vis+ '>>wsclean_chunterings'
    print 'Executing: '+cmd
    os.system (cmd)


def montage_plot(vis):
    pass
