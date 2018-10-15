#!/usr/bin/env python

import numpy as np,os,sys,glob,astropy,argparse,time
import pyrap.tables as pt
import loop3_service
CCRIT=1.6

def imaging(vis,niters,threshold):
    loop3_service.imagr (vis,cellsize='0.05asec',imsize=4096,maxuvl=800000,gain=0.1,\
                  mgain=0.85,dostopnegative=True,niter=niters,autothreshold=threshold,\
                  weightingrankfiltersize=256,weightingrankfilter=3,domultiscale=True,\
                  automask=7.0,outname='test',dolocalrms=True)

# Needs an initial model. May be provided as:
#    model=None        Make an image and use that. (Unlikely to be a good idea)
#    model='MODEL'     Look in the MODEL_DATA column
#    model=[filename]  LOFAR sourcedb format e.g. converted from FIRST, LoTSS, EHT imager....
#                      [Not working yet, maybe require calling routine to do this?]
# zinc is the integration time of input file. Can this be determined automatically?

def selfcal(vis,model='MODEL',outcal_root='',zinc=8.,max_sol=3600.0,init_sol=30.0,\
            incol='DATA',outcol='DATA'):
    if not model:
        imaging(vis,1000,10)
    # need a predict step to deal with sourcedb here if necessary
    ptant = pt.table(vis+'/ANTENNA')
    antenna_list = np.array([],dtype='S')
    for i in ptant.select('NAME'):
        antenna_list = np.append(antenna_list,i.values())
    sol_int_range = np.arange(np.ceil(np.log(max_sol/init_sol)/np.log(3.)))
    sol_int_range = np.ceil(init_sol*3.**sol_int_range/zinc)
    nant, nsol = len(antenna_list), len(sol_int_range)
    coh = CCRIT*np.ones((nsol,nant))
    for i in range(nsol):
        solint = sol_int_range[i]
        outcal_root = outcal_root if len(outcal_root) else vis
        outcal = outcal_root+'_c%d.h5'%i
        print 'Beginning pass with solint %.1f sec' % (solint*zinc)
        loop3_service.calib (vis, solint=solint, outcal=outcal, incol=incol, outcol=outcol)
        loop3_service.snplt (htab=outcal,outpng=outcal)
        coh[i] = loop3_service.coherence_metric (outcal)
        print 'Coherences: '
        for j in range(nant):
            sys.stdout.write('%s:%f '%(antenna_list[j],coh[i,j]))
        if len(coh[i][coh[i]>=CCRIT])==0:
            break
    # For each antenna in the antenna list, find the selfcal table with the shortest        
    # solution interval that contains coherent solutions. If there is no such table, report
    # -1 in order to signal that they should all be set to zero.
    ncoh = np.ones(nant,dtype=int)*-1
    allcoh = np.ones(nant,dtype=float)*CCRIT
    for i in range(nant):
        try:
            ncoh[i] = np.min(np.ravel(np.argwhere(coh[:,i]<CCRIT)))
            allcoh[i] = coh[:,i][np.argwhere(coh[:,i]==ncoh[i])[0][0]]
        except:
            pass
    # For each selfcal table containing the shortest solution interval with coherence on some
    # antennas, replace the entries in the first selfcal table with the interpolated values 
    # from that antenna
    for i in range(1,coh.shape[0]):
        iant = antenna_list[ncoh==i]
        if len(iant):
            loop3_service.clcal (outcal_root+'_c0.h5',outcal_root+'_c%d.h5'%i,ant_interp=iant)
    # For each antenna without any coherence at all, zero the phase solutions for that antenna
    iant = antenna_list[ncoh==-1]
    if len(iant):
        loop3_service.clcal (outcal_root+'_c0.h5',outcal_root+'_c0.h5',ant_interp=iant,dozero=True)
    return allcoh

# following is based on Frits's algorithm with measure_statistic

def measure_statistic (filename):
    img = fits.open(filename)[0].data.squeeze()
    return abs (img.max()/img.min())

def hybridloops (vis,nloops=5):
    import bdsf
    # make an initial model using wsclean and write into the model column (doupdatemodel True)
    # if we have another image from somewhere else (or sourcedb) better to use that
    loop3_service.imagr (vis,threads=6,domultiscale=True,outname=vis+'_im00',channelsout=8,\
           dojoinchannels=True,robust=0.)
    prevstat = 0.0
    for i in range(nloops):
        coh = selfcal(vis,model='MODEL',incol='DATA',outcol='CORRECTED_DATA',outcal_root=vis+'_%02d'%i)
        ####### need a line here to bomb out if no coherence
        img=bdsf.process_image('%s_im%02d-MFS-image.fits'%(vis,i),atrous_do=True,thresh_isl=5.0)
        img.export_image(img_type='island_mask')
        # rename to something shorter
        os.system('mv %s_im%02d-MFS-image.pybdsm_island_mask.fits %s_im%02d-mask.fits'%(vis,i,vis,i))
        # on iteration 0, overwrites the original image
        loop3_service.imagr(vis,cellsize='0.1asec',domultiscale=True,\
                  outname=vis+'_im%02d'%i,dojoinchannels=True,channelsout=8,robust=0,\
                  fitsmask=vis+'_im%02d-mask.fits'%i,dolocalrms=True)
        thisstat = measure_statistic(vis+'_im%02d-MFS-image.fits'%i)
        # exit loop if clean finishing
        if thisstat-prevstat<0.01:
            break
        else:   # otherwise predict image back for next selfcal loop
            loop3_service.imagr(vis,dopredict=True,automask=5.0,autothreshold=2.5,dolocalrms=True,\
                  robust=0,outname=vis+'_im%02d-MFS-image.fits'%i)
 
    
    
