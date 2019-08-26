#!/usr/bin/env python

from astropy.io import fits
import ConfigParser
import glob
import gnirsHeaders
import gnirs_unc
#import matplotlib
import log
import matplotlib.pyplot as plt
import numpy as np
import os
from pyraf import iraf
import utils
#matplotlib.use('Agg')


# ----------------------------------------------------------------------------------------------------------------------
def start(configfile):
    """
    Look for offsets between orders and then combines orders.
    
    Note:  
      1. No order scaling currently done on full-slit or stepwise extractions.
      2. offsets=manual will only allow user to scale orders in regular extraction (no full-slit or stepwise)
    """
    logger = log.getLogger('gnirsCombineOrdersXD.start')

    ########################################################################
    #                                                                      #
    #                  BEGIN - GENERAL COMBINE 2D SETUP                    #
    #                                                                      #
    ########################################################################

    path = os.getcwd()  # Store current working directory for later use.

    logger.info('#################################################')
    logger.info('#                                               #')
    logger.info('#       Start Combining GNIRS XD Orders         #')
    logger.info('#                                               #')
    logger.info('#################################################')

# Set up/prepare IRAF.
    iraf.gemini()
    iraf.gemtools()
    iraf.gnirs()
    iraf.onedspec()
    iraf.imutil()

    # Reset to default parameters the used IRAF tasks.
    iraf.unlearn(iraf.gemini, iraf.gemtools,iraf.gnirs, iraf.imutil)

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
    manualMode = config.getboolean('defaults', 'manualMode')
    overwrite = config.getboolean('defaults', 'overwrite')
    
    # config required for flux calibration
    calculate_snr = config.getboolean('gnirsPipeline', 'Calculate_SNR')
    extractFullSlit = config.getboolean('extractSpectra1D', 'extractFullSlit')
    extractStepwise = config.getboolean('extractSpectra1D', 'extractStepwise')
    extractStepSize = config.getfloat('extractSpectra1D', 'extractStepSize')

    combinedsrc = config.get('runtimeFilenames', 'combinedsrc')
    combinedsky = config.get('runtimeFilenames', 'combinedsky')
    extractRegularPrefix = config.get('runtimeFilenames', 'extractRegularPrefix')
    extractFullSlitPrefix = config.get('runtimeFilenames', 'extractFullSlitPrefix')
    extractStepwisePrefix = config.get('runtimeFilenames', 'extractStepwisePrefix')
    hLinePrefix = config.get('runtimeFilenames', 'hLinePrefix')
    dividedTelContinuumPrefix = config.get('runtimeFilenames', 'dividedTelContinuumPrefix')
    telluricPrefix = config.get('runtimeFilenames', 'telluricPrefix')
    bb_unscaled = config.get('runtimeFilenames', 'bb_unscaled')
    bb_scaled = config.get('runtimeFilenames', 'bb_scaled')
    fluxCalibPrefix = config.get('runtimeFilenames', 'fluxCalibPrefix')
    orderOffsetLog = config.get('runtimeFilenames', 'orderOffsetLog')
    finalPrefix = config.get('runtimeFilenames', 'finalPrefix')
    orderResampledSrc = config.get('runtimeFilenames', 'orderResampledSrc')
    orderResampledSky = config.get('runtimeFilenames', 'orderResampledSky')

    # combineOrdersXD specific config
    redshift = config.get('combineOrdersXD', 'redshift')  # import this as string for check later
    shiftToRestframe = config.getboolean('combineOrdersXD', 'shiftToRestframe')
    offsetCorrectionMethod = config.get('combineOrdersXD', 'offsetCorrectionMethod')
    orderScalingRegions = config.items('orderScalingRegions')
    orderResampling = config.getboolean('combineOrdersXD', 'orderResampling')

    #########################################################################
    #                                                                       #
    #                 COMPLETE - GENERAL COMBINE 2D SETUP                   #
    #                                                                       #
    #########################################################################

    for scipath in config.options('ScienceDirectories'):

        if not config.getboolean("ScienceDirectories", scipath):  # only process directories that are marked True
            logger.info('Skipping flux calibaration in %s', scipath)
            continue

        #########################################################################
        #                                                                       #
        #                  BEGIN - OBSERVATION SPECIFIC SETUP                   #
        #                                                                       #
        #########################################################################

        scipath += '/Intermediate'
        logger.info("Moving to science directory: %s\n", scipath)
        iraf.chdir(scipath)
        finalpath = '../Final/'

        orders = utils.get_orders(scipath)

        logger.info('Checking if input science spectra exist...')
        prefix = fluxCalibPrefix + dividedTelContinuumPrefix + telluricPrefix + extractRegularPrefix
        sci_spectra = ['%s%s_order%d_MEF.fits' % (prefix, utils.nofits(combinedsrc), o) for o in orders]
        logger.debug('sci_spectra: %s', sci_spectra)
        exists = [os.path.exists(f) for f in sci_spectra]
        logger.debug('Exists: %s', exists)
        if all(exists):
            logger.info("Found flux calibrated science source spectra.")
            sci_header_info = gnirsHeaders.info(sci_spectra[0])
            sciName = sci_header_info[os.path.basename(sci_spectra[0])]['OBJECT']
        else:
            logger.warning("Could not find the flux calibrated science spectra.")
            logger.warning("Please run gnirsFluxCalibrate.py or manually put them in %s.", scipath)
            raise SystemExit

        if calculate_snr:  # Check if extracted sky spectra available
            prefix = extractRegularPrefix
            sky_spectra = ['%s%s_order%d_MEF.fits' % (prefix, utils.nofits(combinedsky), o) for o in orders]
            exists = [os.path.exists(f) for f in sky_spectra]
            logger.debug('Exists: %s', exists)
            if all(exists):
                logger.info('Found the extracted sky spectra.')
            else:
                logger.warning("Could not find the extracted sky spectra in %s.", scipath)
                logger.warning("Cannot calculate the S/N.")
                calculate_snr = False

        # TODO:  Read the filenames of all stepwise extracted sepctra.
        if extractFullSlit:
            # nsteps = 1
            pass

        # TODO:  Read the filenames of all stepwise extracted sepctra.
        if extractStepwise:
            # nsteps =
            pass

        # If 'shiftToRestframe' is set, check if redshift is a real number; else set 'shiftToRestframe' to 'False'.
        if shiftToRestframe:
            try:
                redshift = float(redshift)
                logger.info("The science spectra will be shifted to the restframe.")
            except:
                logger.warning("Cannot convert the redshift in the config file to a float.")
                logger.warning("The science spectra cannot be shifted to the rest frame.")
                shiftToRestframe = False


        #########################################################################
        #                                                                       #
        #                 COMPLETE - OBSERVATION SPECIFIC SETUP                 #
        #         BEGIN COMBINING XD ORDERS SPECTRA FOR AN OBSERVATION          #
        #                                                                       #
        #########################################################################

        # First, get the "goodbits"(RM).  These are the range of pixels in each order that should be
        # incorporated into the final, merged spectrum.  These are educated guesses and can be edited by user.
        regions = {}
        for order, r in orderScalingRegions:
            regions[int(order)] = r
        logger.debug('orderScalingRegions: %s', regions)

        if config.getboolean('interactive', 'combine_orders'):
            # Allow the user to attempt to adjust the relative fluxes of each order before combining them
            logger.info("Running specplot so that you can adjust scaling for different orders as necessary...")
            #specplotInputList = makeIRAFinputList(sci_spectra, orders, regions)

            inlist = []
            for i in range(len(orders)):
                inlist.append('%s[1][%s]' % (sci_spectra[i], regions[orders[i]]))
            logger.debug('inlist: %s', inlist)

            iraf.specplot(
                spectra=inlist, apertures='', bands='', dispaxis=1, nsum=1, autolayout='no',
                autoscale='yes', frac=1.0, units='', transform='none', scale=1., offset=0.,step=0, ptype='1',
                labels='user', ulabels='', xlpos=1.02, ylpos=0.0, sysid='yes', yscale='yes', title='', xlabel='',
                ylabel='', xmin='INDEF', ymin='INDEF', xmax='INDEF', ymax='INDEF', logfile=orderOffsetLog,
                graphics='stdgraph', cursor='')

            with open(orderOffsetLog, 'r') as f:
                offsets = f.readlines()
            logger.debug('Spectplot logfile: %s', orderOffsetLog)

            for i in range(len(orders)):
                # Locate the scaling of the order that is written to offsets.log, and multiply by it
                scale = float(offsets[i+5].split()[5])
                logger.debug('scale: %s', scale)
                iraf.imarith(
                    operand1=sci_spectra[i], op='*', operand2=scale,
                    result=finalPrefix+sci_spectra[i], title='', divzero=0.0, hparams='', pixtype='',
                    calctype='', verbose='yes', noact='no')

        else:   # Do not attempt to amke any correction for offsets between the orders
            for s in sci_spectra:
                iraf.copy(input=s, output=finalPrefix+s, verbose='yes')

        # Construct lists of files to go into odcombine and the gnirs_unc.joinorders_noresampling
        # For each extraction, combine the orders into a single spectrum and write out fits and text files

        # NOTE:  You may simply use iraf.odcombine to combine the orders, but odcombine resamples to a linear 
        # dispersion function, oversampling some of the spectrum.  To avoid this, you can use the 
        # joinorders_noresampling routine (by Daniel Duschel Rutra, UFRGS) to join the orders that uses 
        # output from odcombine to find wavelength limits of each order, so we still need to run odcombine.
        iraf.onedspec.interp = 'linear'

        # Set the name of the science as the name of the order combined output spectrum
        finalSpectrum = sciName + '_src.fits'
        finalSpectrum_fullslit = sciName + '_src_FS.fits'
        finalSpectrum_stepwise = sciName + '_step'
        prefix = finalPrefix + fluxCalibPrefix + dividedTelContinuumPrefix + telluricPrefix + extractRegularPrefix
        logger.debug('prefix: %s', prefix)
        logger.debug('combinedsrc: %s', combinedsrc)

        # sci_scaledOrders = sorted(glob.glob(prefix + utils.nofits(combinedsrc) + '_order*_MEF.fits'))
        # odcombineInputList = makeIRAFinputList(sci_scaledOrders, orders, regions)

        inlist = []
        for o in orders:
            inlist.append('%s%s_order%d_MEF[1][%s]' % (prefix, utils.nofits(combinedsrc), o, regions[o]))
        logger.debug('inlist: %s', inlist)

        if orderResampling:
            combineOrdersXD(','.join(inlist), finalpath+finalSpectrum, overwrite)
            iraf.wspectext(input=finalpath+finalSpectrum, output=utils.nofits(finalpath+finalSpectrum)+'.txt',
                header='no', wformat='', mode='al')
        else:
            combineOrdersXD(','.join(inlist), orderResampledSrc, overwrite)

            inlist = []
            for o in orders:
                inlist.append('%s%s_order%d_MEF[1]' % (prefix, utils.nofits(combinedsrc), o))
            logger.debug('inlist: %s', inlist)

            gnirs_unc.joinorders_noresampling(
                inlist=inlist, merged_spec=orderResampledSrc,
                outspec=finalpath+finalSpectrum, writefits=True, snrlist=None)


        raiseSystemExit('STOP!')  # Andy, you stopped here.


        if extractFullSlit:
            sci_scaledOrders = sorted(glob.glob(finalPrefix + fluxCalibPrefix +
                                                dividedTelContinuumPrefix + telluricPrefix + extractFullSlitPrefix + utils.nofits(combinedsrc) +
                                                '_order*_MEF.fits'))
            odcombineInputList = makeIRAFinputList(sci_scaledOrders, orders, regions)
            if orderResampling:
                combineOrdersXD(','.join(odcombineInputList), finalpath+finalSpectrum_fullslit, overwrite)
                iraf.wspectext(input=finalpath+finalSpectrum_fullslit, output=utils.nofits(finalpath+finalSpectrum_fullslit)+'.txt',
                               header='no', wformat='', mode='al')
            else:
                combineOrdersXD(odcombineInputList, orderResampledSrc, overwrite)
                gnirs_unc.joinorders_noresampling(inlist=sci_scaledOrders, merged_spec=orderResampledSrc,
                                                  outspec=finalpath+finalSpectrum_fullslit, writefits=True, snrlist=None)
            pass

        if extractStepwise:
            sci_scaledOrders = sorted(glob.glob(scipath + '/' + finalPrefix + fluxCalibPrefix +
                                                dividedTelContinuumPrefix + telluricPrefix + extractStepwisePrefix + utils.nofits(combinedsrc) +
                                                '_order*_MEF.fits'))
            odcombineInputList = makeIRAFinputList(sci_scaledOrders, orders, regions)
            for k in range (1,steps):
                iraf.odcombine(input=step_dict['step'+str(k)],output='odcombine_output_step'+str(k),headers='',bpmasks='',rejmasks='',nrejmasks='',expmasks='',sigmas='',logfile='STDOUT',apertures='',group='all',first='no',w1='INDEF',w2='INDEF',dw='INDEF',nw='INDEF',log='no',combine='average',reject='none',outtype='real',outlimits='',smaskformat='bpmspectrum',smasktype='none',smaskvalue=0.0,blank=0.0,scale='none',zero='none',weight='none',statsec='',expname='',lthreshold='INDEF',hthreshold='INDEF',nlow=1,nhigh=1,nkeep=1,mclip='yes',lsigma=3.0,hsigma=3.0,rdnoise='0.0',gain='1.0',snoise='0.0',sigscale=0.1,pclip=-0.5,grow=0.0,offsets='physical',masktype='none',maskvalue=0.0,mode='al')
                gnirs_unc.joinorders_noresampling(temp,'odcombine_output_step'+str(k)+'.fits','../PRODUCTS/'+target_name+'_step'+str(k)+'.fits',writefits=True,snrlist=None)
                #iraf.wspectext(input='../PRODUCTS/'+target_name+'_step'+str(k), output='../PRODUCTS/'+target_name+'_step'+str(k)+'.txt',header='no',wformat='',mode='al')
            pass

        # Combine orders for the sky spectrum
        if calculate_snr:
            odcombineInputList = makeIRAFinputList(sky_spectra, orders, regions)
            if orderResampling:
                combineOrdersXD(','.join(odcombineInputList), finalpath+finalSpectrum_fullslit, overwrite)
                iraf.wspectext(input=finalpath+finalSpectrum_fullslit, output=utils.nofits(finalpath+finalSpectrum_fullslit)+'.txt',
                               header='no', wformat='', mode='al')
            else:
                combineOrdersXD(odcombineInputList, orderResampledSrc, overwrite)
                gnirs_unc.joinorders_noresampling(inlist=sci_scaledOrders, merged_spec=orderResampledSrc,
                                                  outspec=finalpath+finalSpectrum_fullslit, writefits=True, snrlist=None)
            iraf.odcombine(skylist,output='odcombine_sky', headers='', bpmasks='', rejmasks='', nrejmasks='', expmasks='', sigmas='',logfile='STDOUT',apertures='',group='all',first='no',w1='INDEF',w2='INDEF',dw='INDEF',nw='INDEF',log='no',combine='average',reject='none',outtype='real',outlimits='',smaskformat='bpmspectrum',smasktype='none',smaskvalue=0.0,blank=0.0,scale='none',zero='none',weight='none',statsec='',expname='',lthreshold='INDEF',hthreshold='INDEF',nlow=1,nhigh=1,nkeep=1,mclip='yes',lsigma=3.0,hsigma=3.0,rdnoise='0.0',gain='1.0',snoise='0.0',sigscale=0.1,pclip=-0.5,grow=0.0,offsets='physical',masktype='none',maskvalue=0.0,mode='al')
            gnirs_unc.joinorders_noresampling(temp,'odcombine_sky.fits', 'sky.fits', writefits=True, snrlist=None)
            #iraf.wspectext(input='sky', output='sky.txt', header='no', wformat='', mode='al')


        # Shift the spectra to rest wavelength if desired
        if shiftToRestframe:
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


        # Thinking it would be a good idea to plot the separate orders so the user can judge if there are any unacceptable offsets and edit the regions used for combining, if they like
        # Requires writing out relevant files as text
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

        # This is a useful plot that users should look at
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
        logger.info("##############################################################################")

    iraf.chdir(path)  # Return to directory script was begun from.

    return


# ----------------------------------------------------------------------------------------------------------------------
def makeIRAFinputList_OLD(inputlist, orders, regions):
    """
    Make a list of comma-separated ONLY files along with regions used for scaling orders to be passed as a list to 
    IRAF commands.
    inputlist: Python list of input files (MEFs with PHU in extension [0] and data in extension[1])
    orders: List of orders
    regions: Python dictionary of good regions corresponding to different orders
    """
    logger = log.getLogger('gnirsCombineOrdersXD.makeIRAFinputList')
    logger.debug('inputlist: %s', inputlist)
    logger.debug('order: %s', orders)
    logger.debug('regions: %s', regions)

    if len(inputlist) == 0:
        logger.error('Input list is empty.')
        raise SystemExit

    outputlist = []
    for step in inputlist:
        for i in range(len(orders)):
            outputlist.append(inputlist[i] + '[1][' + regions[orders[i]] + ']')
    logger.debug('outputlist: %s', outputlist)
    return outputlist


# ----------------------------------------------------------------------------------------------------------------------
def combineOrdersXD(inlist, output, overwrite):
    """
    Combine spectral orders.
    """
    logger = log.getLogger('combineOrdersXD')
    logger.debug('inlist: %s', inlist)
    logger.debug('output: %s', output)

    if os.path.exists(output):
        if overwrite:
            logger.warning("Removing old %s", output)
            os.remove(output)
        else:
            logger.warning("Output exists and overwrite flat not set.  Using existing output for further processing.")
            return

    iraf.odcombine(
        input=inlist, output=output, headers='', bpmasks='', rejmask='', nrejmasks='', expmasks='',
        sigmas='', logfile=logger.root.handlers[0].baseFilename, apertures='', group='all', first='no',
        w1='INDEF', w2='INDEF', dw='INDEF', nw='INDEF', log='no', combine='average', reject='none', outtype='real',
        outlimits='', smaskformat='bpmspectrum', smasktype='none', smaskvalue=0.0, blank=0.0, scale='none',
        zero='none', weight='none', statsec='', expname='', lthreshold='INDEF', hthreshold='INDEF', nlow=1, nhigh=1,
        nkeep=1, mclip='yes', lsigma=3.0, hsigma=3.0, rdnoise='0.0', gain='1.0', snoise='0.0', sigscale=0.1,
        pclip=-0.5, grow=0.0, offsets='physical', masktype='none', maskvalue=0.0, mode='al')

    return


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
