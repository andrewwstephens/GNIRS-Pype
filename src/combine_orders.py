#!/usr/bin/env python

import ConfigParser
import gnirsHeaders
import gnirs_unc
import log
import os
from pyraf import iraf
import utils


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
    iraf.unlearn(iraf.gemini, iraf.gemtools,iraf.gnirs, iraf.imutil)

    # Set clobber to 'yes' for the script. This still does not make the gemini tasks overwritewrite files, so: YOU WILL 
    # LIKELY HAVE TO REMOVE FILES IF YOU RE_RUN THE SCRIPT.
    us_clobber=iraf.envget("clobber")
    iraf.reset(clobber='yes')
    
    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)

    manualMode = config.getboolean('defaults', 'manualMode')
    overwrite = config.getboolean('defaults', 'overwrite')
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
        sci_spectra = utils.make_list(prefix + utils.nofits(combinedsrc), orders=orders, suffix='_MEF.fits')

        if utils.exists(sci_spectra, overwrite=False):
            logger.info("Found flux calibrated science spectra.")
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
                logger.info('Found all the extracted sky spectra.')
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
            prefix = fluxCalibPrefix + dividedTelContinuumPrefix + telluricPrefix + extractRegularPrefix
            iraf.specplot(
                spectra=utils.make_list(prefix + utils.nofits(combinedsrc), regions=regions),
                apertures='', bands='', dispaxis=1, nsum=1, autolayout='no',
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
                    result=finalPrefix + sci_spectra[i],
                    title='', divzero=0.0, hparams='', pixtype='', calctype='', verbose='yes', noact='no')

        else:   # Do not attempt to make any correction for offsets between the orders
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
        finalSpectrum_fullslit = sciName + '_src_fullslit.fits'
        finalSpectrum_stepwise = sciName + '_step'
        finalSky = sciName + '_sky.fits'

        # Standard Extraction ------------------------------------------------------------------------------------------
        logger.info('Combining orders...')
        prefix = finalPrefix + fluxCalibPrefix + dividedTelContinuumPrefix + telluricPrefix + extractRegularPrefix

        if orderResampling:

            odcombine(
                inlist=utils.make_list(prefix + utils.nofits(combinedsrc), regions=regions),
                output=finalpath + finalSpectrum,
                overwrite=overwrite)

            iraf.wspectext(
                input=finalpath+finalSpectrum,
                output=utils.nofits(finalpath+finalSpectrum)+'.txt',
                header='no', wformat='', mode='al')

        else:  # TODO: TEST THIS

            odcombine(
                inlist=utils.make_list(prefix + utils.nofits(combinedsrc), regions=regions),
                output=orderResampledSrc,
                overwrite=overwrite)

            gnirs_unc.joinorders_noresampling(
                inlist=utils.make_list(prefix + utils.nofits(combinedsrc), orders=orders),
                merged_spec=orderResampledSrc,
                outspec=finalpath+finalSpectrum,
                writefits=True, snrlist=None)

        if extractFullSlit:  # Full-Slit Extraction --------------------------------------------------------------------
            logger.info('Combining orders of full-slit extraction...')
            prefix = finalPrefix + fluxCalibPrefix + dividedTelContinuumPrefix + telluricPrefix + extractFullSlitPrefix

            if orderResampling:
                odcombine(
                    inlist=utils.make_list(prefix + utils.nofits(combinedsrc), regions=regions),
                    output=finalpath + finalSpectrum_fullslit,
                    overwrite=overwrite)
                iraf.wspectext(
                    input=finalpath+finalSpectrum_fullslit,
                    output=utils.nofits(finalpath+finalSpectrum_fullslit)+'.txt',
                    header='no', wformat='', mode='al')
            else:
                odcombine(
                    inlist=utils.make_list(prefix + utils.nofits(combinedsrc), regions=regions),
                    output=orderResampledSrc,
                    overwrite=overwrite)

                gnirs_unc.joinorders_noresampling(
                    inlist=utils.make_list(prefix + utils.nofits(combinedsrc), orders=orders),
                    merged_spec=orderResampledSrc,
                    outspec=finalpath+finalSpectrum_fullslit,
                    writefits=True, snrlist=None)

        if extractStepwise:  # Step-wise Extraction --------------------------------------------------------------------
            logger.info('Combining orders of full-slit extraction...')
            prefix = finalPrefix + fluxCalibPrefix + dividedTelContinuumPrefix + telluricPrefix + extractStepwisePrefix
            logger.error('Stepwise extraction not yet supported')  # TODO: Add support or stepwise extraction
            raise SystemExit

        if calculate_snr:      # Combine orders of the SKY spectrum
            prefix = extractRegularPrefix

            if orderResampling:
                odcombine(
                    inlist=utils.make_list(prefix + utils.nofits(combinedsky),regions=regions),
                    output=finalpath + finalSky,
                    overwrite=overwrite)
                iraf.wspectext(
                    input=finalpath+finalSky,
                    output=utils.nofits(finalpath+finalSky)+'.txt',
                    header='no', wformat='', mode='al')

            else:
                odcombine(
                    inlist=utils.make_list(prefix + utils.nofits(combinedsky), regions=regions),
                    output=orderResampledSky,
                    overwrite=overwrite)

                gnirs_unc.joinorders_noresampling(
                    inlist=utils.make_list(prefix + utils.nofits(combinedsrc), orders=orders),
                    merged_spec=orderResampledSky,
                    outspec=finalpath+finalSky,
                    writefits=True, snrlist=None)

        if shiftToRestframe:  # Shift spectra to the rest wavelength if desired

            # Standard extraction --------------------------------------------------------------------------------------
            iraf.dopcor(
                input=finalpath + finalSpectrum,
                output=finalpath + finalSpectrum + '_rest',
                redshift=redshift, isvelocity='no', add='no', dispersion='yes',
                flux='no', factor='3.0', apertures='', verbose='no')

            iraf.wspectext(
                input=finalpath + finalSpectrum + '_rest',
                output=finalpath + finalSpectrum + '_rest.txt',
                header='no', wformat='')

            if calculate_snr:  # Sky spectrum

                iraf.dopcor(
                    input='sky',
                    output='sky_rest',
                    redshift=redshift, isvelocity='no', add='no', dispersion='yes',
                    flux='no', factor='3.0', apertures='', verbose='no')

                iraf.wspectext(
                    input='sky_rest',
                    output='sky_rest.txt',
                    header='no', wformat='')

            if extractFullSlit:  # Full-slit Extraction ----------------------------------------------------------------

                iraf.dopcor(
                    input=finalpath + finalSpectrum + '_fullslit',
                    output=finalpath + finalSpectrum + '_fullslit_rest',
                    redshift=redshift, isvelocity='no', add='no', dispersion='yes', flux='no',
                    factor='3.0', apertures='', verbose='no')

                iraf.wspectext(
                    input=finalpath + finalSpectrum + '_fullslit_rest',
                    output=finalpath + finalSpectrum + '_fullslit_rest.txt',
                    header='no', wformat='')

            if extractStepwise:  # Step-wise Extraction ----------------------------------------------------------------
                for k in range(1, steps):

                    iraf.dopcor(
                        input=finalpath + finalSpectrum + '_step' + str(k),
                        output=finalpath + finalSpectrum + '_step' + str(k) + '_rest',
                        redshift=redshift, isvelocity='no', add='no', dispersion='yes', flux='no',
                        factor='3.0', apertures='', verbose='no')

                    iraf.wspectext(
                        input=finalpath + finalSpectrum + '_step' + str(k) + '_rest',
                        output=finalpath + finalSpectrum + '_step' + str(k) + '_rest.txt',
                        header='no', wformat='')

        # Write ASCII output of the individual orders for later plotting:
        prefix = finalPrefix + fluxCalibPrefix + dividedTelContinuumPrefix + telluricPrefix + extractRegularPrefix
        output = utils.make_list(prefix + utils.nofits(combinedsrc), orders=orders, suffix='.txt')
        if not utils.exists(output, overwrite):
            for f in utils.make_list(prefix + utils.nofits(combinedsrc), orders=orders, suffix=''):
                logger.debug('%s -> %s', f + '_MEF[1]', f + '.txt')
                iraf.wspectext(
                    input=f + '_MEF[1]',
                    output=f + '.txt',
                    header='no', wformat='')

        if extractFullSlit:
            prefix = finalPrefix + fluxCalibPrefix + dividedTelContinuumPrefix + telluricPrefix + extractFullSlitPrefix
            output = utils.make_list(prefix + utils.nofits(combinedsrc), orders=orders, suffix='.txt')
            if not utils.exists(output, overwrite):
                for f in utils.make_list(prefix + utils.nofits(combinedsrc), orders=orders, suffix=''):
                    logger.debug('%s -> %s', f + '_MEF[1]', f + '.txt')
                    iraf.wspectext(
                        input=f + '_MEF[1]',
                        output=f + '.txt',
                        header='no', wformat='')

        if extractStepwise:
            for k in range(1, steps):
                prefix = finalPrefix + fluxCalibPrefix + dividedTelContinuumPrefix + telluricPrefix + extractStepwisePrefix
                output = utils.make_list(prefix + utils.nofits(combinedsrc), orders=orders, suffix='.txt')
                if not utils.exists(output, overwrite):
                    for f in utils.make_list(prefix + utils.nofits(combinedsrc), orders=orders, suffix=''):
                        logger.debug('%s -> %s', f + '_MEF[1]', f + '.txt')
                        iraf.wspectext(
                            input=f + '_MEF[1]',
                            output=f + '.txt',
                            header='no', wformat='')

        logger.info("##############################################################################")
        logger.info("#                                                                            #")
        logger.info("#  COMPLETE - Combining 2D spectra completed for                             #")
        logger.info("#  %s", scipath)
        logger.info("#                                                                            #")
        logger.info("##############################################################################")

    iraf.chdir(path)  # Return to directory script was begun from.

    return


# ----------------------------------------------------------------------------------------------------------------------
def odcombine(inlist, output, overwrite=False):
    """
    Combine spectral orders with IRAF odcombine
    :param inlist: list of input files
    :param output: output fits file (string)
    :param overwrite:
    :return:
    """
    logger = log.getLogger('odcombine')
    logger.debug('inlist: %s', inlist)
    logger.debug('output: %s', output)
    logger.debug('overwrite: %s', overwrite)

    if utils.exists([output], overwrite):
        logger.warning("Output exists and overwrite = False so skipping this step.")
        return

    iraf.odcombine(
        input=','.join(inlist), output=output, headers='', bpmasks='', rejmask='', nrejmasks='', expmasks='',
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
