#!/usr/bin/env python

# MIT License

# Copyright (c) 2015, 2017 Marie Lemoine-Busserolle

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

################################################################################
#                Import some useful Python utilities/modules                   #
################################################################################

import os, sys, glob, log, ConfigParser, gnirsHeaders, gnirs_unc
from pyraf import iraf
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
# NOTE:  Import matplotlib to avoid (in cron) 'tkinter.TclError: no display name and no $DISPLAY' error 
# (jmiller,mpohlen; September 22, 2015)
import matplotlib
matplotlib.use('Agg')

#---------------------------------------------------------------------------------------------------------------------#

def start(configfile):
    """
    Looks for offsets between orders and then, combines orders.
    
    Note:  
      1. No order scaling currently done on full-slit or stepwise extractions.
      2. offsets=manual will only allow user to scale orders in regular extraction (no full-slit or stepwise)
    """
    logger = log.getLogger('gnirsCombineOrdersXD.start')

    ###########################################################################
    ##                                                                       ##
    ##                  BEGIN - GENERAL COMBINE 2D SETUP                     ##
    ##                                                                       ##
    ###########################################################################


    # Store current working directory for later use.
    path = os.getcwd()

    logger.info('#################################################')
    logger.info('#                                               #')
    logger.info('#       Start Combining GNIRS XD Orders         #')
    logger.info('#                                               #')
    logger.info('#################################################\n')

# Set up/prepare IRAF.
    iraf.gemini()
    iraf.gemtools()
    iraf.gnirs()
    iraf.onedspec()
    iraf.imutil()

    # Reset to default parameters the used IRAF tasks.
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs
    ,iraf.imutil)

    # From http://bishop.astro.pomona.edu/Penprase/webdocuments/iraf/beg/beg-image.html:
    # Before doing anything involving image display the environment variable stdimage must be set to the correct frame 
    # buffer size for the display servers (as described in the dev$graphcap file under the section "STDIMAGE devices") 
    # or to the correct image display device. The task GDEVICES is helpful for determining this information for the 
    # display servers.
    iraf.set(stdimage='imt1024')

    # Prepare the IRAF package for GNIRS.
    # NSHEADERS lists the header parameters used by the various tasks in the GNIRS package (excluding headers values 
    # which have values fixed by IRAF or FITS conventions).
    iraf.nsheaders("gnirs",logfile=logger.root.handlers[0].baseFilename)

    # Set clobber to 'yes' for the script. This still does not make the gemini tasks overwritewrite files, so: YOU WILL 
    # LIKELY HAVE TO REMOVE FILES IF YOU RE_RUN THE SCRIPT.
    us_clobber=iraf.envget("clobber")
    iraf.reset(clobber='yes')
    
    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)
    
    # Read general config
    manualMode = config.getboolean('defaults','manualMode')
    overwrite = config.getboolean('defaults','overwrite')
    
    # config required for flux calibration
    calculateSpectrumSNR = config.getboolean('gnirsPipeline', 'calculateSpectrumSNR')
    extractFullSlit = config.getboolean('extractSpectra1D','extractFullSlit')
    extractStepwise = config.getboolean('extractSpectra1D','extractStepwise')
    extractStepSize = config.getfloat('extractSpectra1D','extractStepSize')

    combinedsrc = config.get('runtimeFilenames', 'combinedsrc')
    combinedsky = config.get('runtimeFilenames', 'combinedsky')
    extractRegularPrefix = config.get('runtimeFilenames', 'extractRegularPrefix')
    extractFullSlitPrefix = config.get('runtimeFilenames', 'extractFullSlitPrefix')
    extractStepwisePrefix = config.get('runtimeFilenames', 'extractStepwisePrefix')
    hLinePrefix = config.get('runtimeFilenames', 'hLinePrefix')
    dividedTelContinuumPrefix = config.get('runtimeFilenames', 'dividedTelContinuumPrefix')
    telluricPrefix = config.get('runtimeFilenames', 'telluricPrefix')
    bb_unscaled = config.get('runtimeFilenames', 'bb_unscaled')
    bb_scaled = config.get('runtimeFilenames','bb_scaled')
    fluxCalibPrefix = config.get('runtimeFilenames', 'fluxCalibPrefix')
    orderOffsetLog = config.get('runtimeFilenames', 'orderOffsetLog')
    orderScaledPrefix = config.get('runtimeFilenames', 'orderScaledPrefix')
    orderResampledSrc = config.get('runtimeFilenames', 'orderResampledSrc')
    orderResampledSky = config.get('runtimeFilenames', 'orderResampledSky')
    
    # combineOrdersXD specific config
    redshift = config.get('combineOrdersXD','redshift')  # import this as string for check later
    shiftToRestframe = config.getboolean('combineOrdersXD','shiftToRestframe')
    offsetCorrectionMethod = config.get('combineOrdersXD','offsetCorrectionMethod')
    orderScalingRegions = config.items('combineOrdersXD','orderScalingRegions')
    orderResampling = config.getboolean('combineOrdersXD','orderResampling')

    ###########################################################################
    ##                                                                       ##
    ##                 COMPLETE - GENERAL COMBINE 2D SETUP                   ##
    ##                                                                       ##
    ###########################################################################

    for scipath in config.options('ScienceDirectories'):
        
        if not config.getboolean("ScienceDirectories", scipath):  # only process directories that are marked True
            logger.info('Skipping flux calibaration in %s', scipath)
            continue

        ###########################################################################
        ##                                                                       ##
        ##                  BEGIN - OBSERVATION SPECIFIC SETUP                   ##
        ##                                                                       ##
        ###########################################################################
        
        scipath += '/Intermediate'
        logger.info("Moving to science directory: %s\n", scipath)
        os.chdir(scipath)
        iraf.chdir(scipath)

        finalpath = '/Final'

        # Check if required sci and std source spectra available in their respective paths
        logger.info("Checking if required science spectra available in %s", scipath)

        sci_fluxCalibrated = []
        sci_fluxCalibrated.append(sorted(glob.glob(scipath + '/' + fluxCalibPrefix + dividedTelContinuumPrefix + \
            telluricPrefix + extractRegularPrefix + nofits(combinedsrc) + '_order*_MEF.fits')))
        if len(sci_fluxCalibrated) > 0:
            logger.info("Required flux calibrated science source spectra available.")
            sci_header_info = gnirsHeaders.info(sci_telluricCorrected[0])
            sciName = sci_header_info[os.path.basename(sci_fluxCalibrated[0])]['OBJECT']
        else:
            logger.warning("Required flux calibrated science source spectra not available.")
            logger.warning("Please run gnirsFluxCalibrate.py to create the telluric corected spectra or provide them")
            logger.warning("manually in %s.", scipath)
            logger.warning("Exiting script.\n")
            raise SystemExit

        if calculateSpectrumSNR:
            # Check if extracted sky images available for science
            sky_extracted_spectra = []
            sky_extracted_spectra.append(sorted(glob.glob(scipath + '/' + extractRegularPrefix + nofits(combinedsky) + \
                '.fits')))
            if len(sky_extracted_spectra) > 0:
                logger.info("Required extracted science sky spectra available.")
            else:
                logger.warning("Parameter 'calculateSpectrumSNR' is 'True', but required extracted science sky")
                logger.warning("spectra not available in %s.", scipath)
                logger.warning("Setting the 'calculateSpectrumSNR' parameter for the current set of observations to")
                logger.warning("'False'.\n")
                calculateSpectrumSNR = False

        logger.info("Required science spectra check complete.\n")

        # TODO:  Read the filenames of all stepwise extracted sepctra.
        if extractFullSlit:
#            nsteps = 1
            pass

        # TODO:  Read the filenames of all stepwise extracted sepctra.
        if extractStepwise:
#            nsteps = 
            pass

        # If the parameter 'shiftToRestframe' is set, check redshift of the target is a real number; else set 
        # 'shiftToRestframe' to 'False'.
        if shiftToRestframe:
            try:
                redshift = float(redshift)  # Test if it a real number
                logger.info("Specified redshift is a real number.  The science spectra will be shifted to the rest") 
                logger.info("frame.")
            except:
                logger.warning("Specified redshift is not a real number, but the parameter 'shiftToRestframe' is set")
                logger.warning("to 'True'.  Setting 'shiftToRestframe' to False.  The science spectra will not be")
                logger.warning("shifted to the rest frame.")
                shiftToRestframe = False

        # Record the right number of order expected according to the GNIRS XD configuration.
        if 'LB_SXD' in scipath:
            orders = [3, 4, 5]
        elif 'LB_LXD' in scipath:
            orders = [3, 4, 5, 6, 7, 8]
        elif 'SB_SXD' in scipath:
            orders = [3, 4, 5, 6, 7, 8]
        else:
            logger.error("#####################################################################################")
            logger.error("#####################################################################################")
            logger.error("#                                                                                   #")
            logger.error("#   ERROR in combining XD orders: unknown GNIRS XD configuration. Exiting script.   #")
            logger.error("#                                                                                   #")
            logger.error("#####################################################################################")
            logger.error("#####################################################################################\n")
            raise SystemExit

        ###########################################################################
        ##                                                                       ##
        ##                 COMPLETE - OBSERVATION SPECIFIC SETUP                 ##
        ##         BEGIN COMBINING XD ORDERS SPECTRA FOR AN OBSERVATION          ##
        ##                                                                       ##
        ###########################################################################

        # First, get the "good" regions of the orders.  These are the range of pixels in each order that should be 
        # incorporated into the final, merged spectrum.  These are best educated guesses and can be edited by user.
        regions = {}
        for order, r in orderScalingRegions:
            regions[int(order)] = r

        if not offsetCorrectionMethod:  # Not correcting for offsets between orders
            # Simply copy the flux calibrated science spectra as final spectra so that all files have same name later
            for i in range(len(orders)):
                iraf.copy(input=sci_fluxCalibrated[i], output=finalPrefix+sci_fluxCalibrated[i], verbose='yes')

        if offsetCorrectionMethod == 'manual':
            # Allow the user to attempt to adjust the relative fluxes of each order before combining them
            specplotInputList = makeIRAFinputList(sci_fluxCalibrated, orders, regions)
            logger.info("Entering iraf.specplot so that the user can adjust scaling for different orders as necessary.")
            iraf.specplot(spectra=specplotInputList, apertures='', bands='', dispaxis=1, nsum=1, autolayout='no', 
                autoscale='yes', frac=1.0, units='', transform='none', scale=1., offset=0.,step=0, ptype='1', 
                labels='user', ulabels='', xlpos=1.02, ylpos=0.0, sysid='yes', yscale='yes', title='', xlabel='', 
                ylabel='', xmin='INDEF', ymin='INDEF', xmax='INDEF', ymax='INDEF', logfile=orderOffsetLog, 
                graphics='stdgraph', cursor='')
            orderOffsets = open(orderOffsetLog,'r')
            offsets = orderOffsets.readlines()
            for i in range(len(orders)):
            # Locate the scaling of the order that is written to offsets.log, and multiply by it
                scale = float(offsets[i+5].split()[5])
                iraf.imarith(operand1=sci_fluxCalibrated[i], op='*', operand2=scale, 
                    result=orderScaledPrefix+sci_fluxCalibrated[i], title='', divzero=0.0, hparams='', pixtype='', 
                    calctype='', verbose='yes', noact='no')

        # Construct lists of files to go into odcombine and the gnirs_unc.joinorders_noresampling
        # For each extraction, combine the orders into a single spectrum and write out fits and text files
        
        # NOTE:  You may simply use iraf.odcombine to combine the orders, but odcombine resamples to a linear 
        # dispersion function, oversampling some of the spectrum.  To avoid this, you can use the 
        # joinorders_noresampling routine (by Daniel Duschel Rutra, UFRGS) to join the orders that uses 
        # output from odcombine to find wavelength limits of each order, so we still need to run odcombine.
        iraf.onedspec.interp='linear'

        # Set the name of the science as the name of the order combined output spectrum
        finalSpectrum = sciName + '_src.fits'
        finalSpectrum_fullslit = sciName + '_src_FS.fits'
        finalSpectrum_stepwise = sciName + '_step'

        sci_scaledOrders = []
        sci_scaledOrders.append(sorted(glob.glob(scipath + '/' + orderScaledPrefix + fluxCalibPrefix + \
            dividedTelContinuumPrefix + telluricPrefix + extractRegularPrefix + nofits(combinedsrc) + \
            '_order*_MEF.fits')))
        odcombineInputList = makeIRAFinputList(sci_scaledOrders, orders, regions)
        if orderResampling:
            combineOrdersXD(','.join(odcombineInputList), finalSpectrum, finalpath, overwrite)
            iraf.wspectext(input=finalpath+finalSpectrum, output=nofits(finalpath+finalSpectrum)+'.txt', header='no', 
                wformat='', mode='al')
        else:
            combineOrdersXD(odcombineInputList, orderResampledSrc, scipath, overwrite)
            gnirs_unc.joinorders_noresampling(inlist=sci_scaledOrders, merged_spec=orderResampledSrc, \
                outspec=finalpath+finalSpectrum, writefits=True, snrlist=None)

        if extractFullSlit:
            sci_scaledOrders = []
            sci_scaledOrders.append(sorted(glob.glob(scipath + '/' + orderScaledPrefix + fluxCalibPrefix + \
                dividedTelContinuumPrefix + telluricPrefix + extractFullSlitPrefix + nofits(combinedsrc) + \
                '_order*_MEF.fits')))
            odcombineInputList = makeIRAFinputList(sci_scaledOrders, orders, regions)
            if orderResampling:
                combineOrdersXD(','.join(odcombineInputList), finalSpectrum_fullslit, finalpath, overwrite)
                iraf.wspectext(input=finalpath+finalSpectrum_fullslit, output=nofits(finalpath+finalSpectrum_fullslit)+'.txt', 
                    header='no', wformat='', mode='al')
            else:
                combineOrdersXD(odcombineInputList, orderResampledSrc, scipath, overwrite)
                gnirs_unc.joinorders_noresampling(inlist=sci_scaledOrders, merged_spec=orderResampledSrc, \
                    outspec=finalpath+finalSpectrum_fullslit, writefits=True, snrlist=None)
            pass

        if extractStepwise:
            sci_scaledOrders = []
            sci_scaledOrders.append(sorted(glob.glob(scipath + '/' + orderScaledPrefix + fluxCalibPrefix + \
                dividedTelContinuumPrefix + telluricPrefix + extractStepwisePrefix + nofits(combinedsrc) + \
                '_order*_MEF.fits')))
            odcombineInputList = makeIRAFinputList(sci_scaledOrders, orders, regions)
            for k in range (1,steps):
                iraf.odcombine(input=step_dict['step'+str(k)],output='odcombine_output_step'+str(k),headers='',bpmasks='',rejmasks='',nrejmasks='',expmasks='',sigmas='',logfile='STDOUT',apertures='',group='all',first='no',w1='INDEF',w2='INDEF',dw='INDEF',nw='INDEF',log='no',combine='average',reject='none',outtype='real',outlimits='',smaskformat='bpmspectrum',smasktype='none',smaskvalue=0.0,blank=0.0,scale='none',zero='none',weight='none',statsec='',expname='',lthreshold='INDEF',hthreshold='INDEF',nlow=1,nhigh=1,nkeep=1,mclip='yes',lsigma=3.0,hsigma=3.0,rdnoise='0.0',gain='1.0',snoise='0.0',sigscale=0.1,pclip=-0.5,grow=0.0,offsets='physical',masktype='none',maskvalue=0.0,mode='al')
                gnirs_unc.joinorders_noresampling(temp,'odcombine_output_step'+str(k)+'.fits','../PRODUCTS/'+target_name+'_step'+str(k)+'.fits',writefits=True,snrlist=None)
                #iraf.wspectext(input='../PRODUCTS/'+target_name+'_step'+str(k), output='../PRODUCTS/'+target_name+'_step'+str(k)+'.txt',header='no',wformat='',mode='al')
            pass

        # Combine orders for the sky spectrum
        if calculateSpectrumSNR:
            odcombineInputList = makeIRAFinputList(sky_extracted_spectra, orders, regions)
            if orderResampling:
                combineOrdersXD(','.join(odcombineInputList), finalSpectrum_fullslit, finalpath, overwrite)
                iraf.wspectext(input=finalpath+finalSpectrum_fullslit, output=nofits(finalpath+finalSpectrum_fullslit)+'.txt', 
                    header='no', wformat='', mode='al')
            else:
                combineOrdersXD(odcombineInputList, orderResampledSrc, scipath, overwrite)
                gnirs_unc.joinorders_noresampling(inlist=sci_scaledOrders, merged_spec=orderResampledSrc, \
                    outspec=finalpath+finalSpectrum_fullslit, writefits=True, snrlist=None)
            iraf.odcombine(skylist,output='odcombine_sky',headers='',bpmasks='',rejmasks='',nrejmasks='',expmasks='',sigmas='',logfile='STDOUT',apertures='',group='all',first='no',w1='INDEF',w2='INDEF',dw='INDEF',nw='INDEF',log='no',combine='average',reject='none',outtype='real',outlimits='',smaskformat='bpmspectrum',smasktype='none',smaskvalue=0.0,blank=0.0,scale='none',zero='none',weight='none',statsec='',expname='',lthreshold='INDEF',hthreshold='INDEF',nlow=1,nhigh=1,nkeep=1,mclip='yes',lsigma=3.0,hsigma=3.0,rdnoise='0.0',gain='1.0',snoise='0.0',sigscale=0.1,pclip=-0.5,grow=0.0,offsets='physical',masktype='none',maskvalue=0.0,mode='al')
            gnirs_unc.joinorders_noresampling(temp,'odcombine_sky.fits','sky.fits',writefits=True,snrlist=None)
            #iraf.wspectext(input='sky',output='sky.txt',header='no',wformat='',mode='al')

            
        #Shift the spectra to rest (if desired by user)
        if doshift == "yes":
            #Standard extraction
            iraf.dopcor (input='../PRODUCTS/'+target_name, output='../PRODUCTS/'+target_name+'_rest', redshift=redshift, isvelocity='no', add='no', dispersion='yes', flux='no', factor='3.0', apertures='', verbose='no', mode='al')
            iraf.wspectext(input='../PRODUCTS/'+target_name+'_rest',output='../PRODUCTS/'+target_name+'_rest.txt',header='no',wformat='',mode='al')
            #Sky spectrum
            if error_spectrum == 'yes':
                iraf.dopcor (input='sky', output='sky_rest', redshift=redshift, isvelocity='no', add='no', dispersion='yes', flux='no', factor='3.0', apertures='', verbose='no', mode='al')
                iraf.wspectext(input='sky_rest',output='sky_rest.txt',header='no',wformat='',mode='al')
            if extras == "yes":
                #Full-slit extraction
                iraf.dopcor (input='../PRODUCTS/'+target_name+'_fullslit', output='../PRODUCTS/'+target_name+'_fullslit_rest', redshift=redshift, isvelocity='no', add='no', dispersion='yes', flux='no', factor='3.0', apertures='', verbose='no', mode='al')
                iraf.wspectext(input='../PRODUCTS/'+target_name+'_fullslit_rest',output='../PRODUCTS/'+target_name+'_fullslit_rest.txt',header='no',wformat='',mode='al')
                #Step-by-step extractions
                for k in range (1,steps):
                    iraf.dopcor (input='../PRODUCTS/'+target_name+'_step'+str(k), output='../PRODUCTS/'+target_name+'_step'+str(k)+'_rest', redshift=redshift, isvelocity='no', add='no', dispersion='yes', flux='no', factor='3.0', apertures='', verbose='no', mode='al')
                    iraf.wspectext(input='../PRODUCTS/'+target_name+'_step'+str(k)+'_rest', output='../PRODUCTS/'+target_name+'_step'+str(k)+'_rest.txt',header='no',wformat='',mode='al')
            

        #Thinking it would be a good idea to plot the separate orders so the user can judge if there are any unacceptable offsets and edit the regions used for combining, if they like
        #Requires writing out relevant files as text
        for j in range (1,7):
            if os.path.isfile('final6.txt'):
                iraf.delete (files="final"+str(j)+'.txt',verify="no",Stdout="/dev/null") 
            iraf.wspectext(input="final"+str(j),output="final"+str(j)+'.txt',header='no',wformat='',mode='al')
            if extras == "yes":
                if os.path.isfile('flamfull6.txt'):
                    iraf.delete (files="flamfull"+str(j)+'.txt',verify="no",Stdout="/dev/null") 
                iraf.wspectext(input="flamfull"+str(j),output="flamfull"+str(j)+'.txt',header='no',wformat='',mode='al')
                for k in range (1,steps):
                    if os.path.isfile('flamstep'+str(steps-1)+'_6.txt'):
                        iraf.delete (files='flamstep'+str(k)+'_'+str(j)+'.txt',verify="no",Stdout="/dev/null") 
                    iraf.wspectext(input='flamstep'+str(k)+'_'+str(j),output='flamstep'+str(k)+'_'+str(j)+'.txt',header='no',wformat='',mode='al')
                    
        #This is a useful plot that users should look at        
        make_orders_fig ("final", "../PRODUCTS/orders.pdf")
        if extras == 'yes':
            make_orders_fig ("flamfull", "../PRODUCTS/orders_fullslit.pdf")
            for k in range (1,steps):
                make_orders_fig ("flamstep"+str(k)+"_", "../PRODUCTS/orders_step"+str(k)+".pdf")

        file.close()

        logger.info("##############################################################################")
        logger.info("#                                                                            #")
        logger.info("#  COMPLETE - Combining 2D spectra completed for                             #")
        logger.info("#  %s", scipath)
        logger.info("#                                                                            #")
        logger.info("##############################################################################\n")

    # Return to directory script was begun from.
    os.chdir(path)
    iraf.chdir(path)

    return


##################################################################################################################
#                                                     ROUTINES                                                   #
##################################################################################################################

def makeIRAFinputList(inputlist, orders, regions):
    """
    Make a list of comma-separated ONLY files along with regions used for scaling orders to be passed as a list to 
    IRAF commands.
    inputlist: Python list of input files (MEFs with PHU in extension [0] and data in extension[1])
    orders: List of orders
    regions: Python dictionary of good regions corresponding to different orders
    """
    logger = log.getLogger('gnirsCombineOrdersXD.makeIRAFinputList')

    outputlist = []
    outputlist_regions = []
    for step in inputlist:
        for i in range(len(orders)):
            outputlist.append(inputlist[i] + '[1][' + regions[orders[i]] + ']')
    return outputlist

#----------------------------------------------------------------------------------------------------------------------#

def combineOrdersXD(inlist, combinedspec, combinedspec_path, overwrite):
    """
    Combine spectral orders.
    """
    logger = log.getLogger('gnirsCombineOrdersXD.combineOrdersXD')

    if os.path.exists(combinedspec_path+'/'+combinedspec):
        if overwrite:
            logger.warning("Removing old %s", combinedspec_path+'/'+combinedspec)
            os.remove(combinedspec_path+'/'+combinedspec)
        else:
            logger.warning("Output exists and -overwrite not set - using the existing output for further processing.")
            return
    iraf.odcombine(input=odcombineInputList, output=combinedspec_path+'/'+combinedspec, headers='', bpmasks='', 
        rejmask='', nrejmasks='', expmasks='', sigmas='', logfile=logger.root.handlers[0].baseFilename, apertures='', 
        group='all', first='no', w1='INDEF', w2='INDEF', dw='INDEF', nw='INDEF', log='no', combine='average', 
        reject='none', outtype='real', outlimits='', smaskformat='bpmspectrum', smasktype='none', smaskvalue=0.0, 
        blank=0.0, scale='none', zero='none', weight='none', statsec='', expname='', lthreshold='INDEF', 
        hthreshold='INDEF', nlow=1, nhigh=1, nkeep=1, mclip='yes', lsigma=3.0, hsigma=3.0, rdnoise='0.0', gain='1.0', 
        snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, offsets='physical', masktype='none', maskvalue=0.0, 
        mode='al')

#----------------------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
