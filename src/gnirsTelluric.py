#!/usr/bin/env python

from astropy.io import fits
import ConfigParser
import glob
import log
import matplotlib.pyplot as plt
import os
from pyraf import iraf
import utils


# ----------------------------------------------------------------------------------------------------------------------
def start(configfile):
    """
    Do a Telluric correction using the IRAF TELLURIC task.
    """
    logger = log.getLogger('gnirsTelluric')

    path = os.getcwd()  # Store current working directory for later use


    # Set up/prepare IRAF.
    iraf.gemini()
    iraf.gemtools()
    iraf.gnirs()

    # Reset to default parameters the used IRAF tasks.
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs)

    # Prepare the IRAF package for GNIRS.
    # NSHEADERS lists the header parameters used by the various tasks in the GNIRS package (excluding headers values 
    # which have values fixed by IRAF or FITS conventions).
    iraf.nsheaders("gnirs", logfile=logger.root.handlers[0].baseFilename)

    # Set clobber to 'yes' for the script. This still does not make the gemini tasks overwrite files, so: YOU WILL 
    # LIKELY HAVE TO REMOVE FILES IF YOU RE_RUN THE SCRIPT.
    us_clobber=iraf.envget("clobber")
    iraf.reset(clobber='yes')
    
    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)

    manualMode = config.getboolean('defaults', 'manualMode')
    overwrite = config.getboolean('defaults', 'overwrite')
    runtimedata = config.get('defaults', 'runtimeData')
    hLineInter = config.getboolean('interactive', 'hLineInter')
    continuumInter = config.getboolean('interactive', 'continuumInter')
    telluricInter = config.getboolean('interactive', 'telluricInter')
    tempInter = config.getboolean('interactive', 'tempInter')
    extractFullSlit = config.getboolean('extractSpectra1D','extractFullSlit')
    extractStepwise = config.getboolean('extractSpectra1D', 'extractStepwise')
    extractStepSize = config.getfloat('extractSpectra1D','extractStepSize')

    combinedsrc = config.get('runtimeFilenames', 'combinedsrc')
    extractRegularPrefix = config.get('runtimeFilenames', 'extractRegularPrefix')
    extractFullSlitPrefix = config.get('runtimeFilenames', 'extractFullSlitPrefix')
    extractStepwiseTracePrefix = config.get('runtimeFilenames', 'extractStepwiseTracePrefix')
    extractStepwisePrefix = config.get('runtimeFilenames', 'extractStepwisePrefix')
    hLinePrefix = config.get('runtimeFilenames', 'hLinePrefix')
    fitTelContinuumPrefix = config.get('runtimeFilenames', 'fitTelContinuumPrefix')
    dividedTelContinuumPrefix = config.get('runtimeFilenames', 'dividedTelContinuumPrefix')
    telluricPrefix = config.get('runtimeFilenames', 'telluricPrefix')

    startstep = config.getint('telluricCorrection', 'Start')
    stopstep = config.getint('telluricCorrection', 'Stop')
    hLineMethod = config.get('telluricCorrection', 'hLineMethod')

    hLineRegions = utils.dictify(config.items('hLineRegions'), fmt='int')
    continuumRegions = utils.dictify(config.items('continuumRegions'), fmt='int')
    telluricRegions = utils.dictify(config.items('telluricRegions'), fmt='int')
    telluricFitOrders = utils.dictify(config.items('TelluricFitOrders'), fmt='int')

    logger.debug('hLineRegions: %s', hLineRegions)
    logger.debug('continuumRegions: %s', continuumRegions)
    logger.debug('telluricRegions: %s', telluricRegions)
    logger.debug('telluricFitOrders: %s', telluricFitOrders)

    for scipath in config.options("ScienceDirectories"):

        if not config.getboolean("ScienceDirectories", scipath):  # only process directories that are marked True
            logger.debug('Skipping %s', scipath)
            continue

        logger.info(' --------------------- ')
        logger.info('| Telluric Correction |')
        logger.info(' --------------------- ')

        scipath += '/Intermediate'
        logger.info("%s", scipath)
        iraf.chdir(scipath)

        telpath = '../Telluric/Intermediate'
        logger.info("Telluric directory: %s", telpath)
        logger.info("Runtime data path: %s", runtimedata)

        sci_spec = scipath + '/' + extractRegularPrefix + combinedsrc
        tel_spec = telpath + '/' + extractRegularPrefix + combinedsrc
        vega_spec = runtimedata + 'vega_ext.fits'
        utils.requires([sci_spec, tel_spec, vega_spec])

        sci_airmass = fits.getheader(sci_spec)['AIRMASS']

        # TODO: support step-wise extraction

        orders = utils.get_orders(scipath)

        # Define output filenames with prefixes added at different stages of script but without the '.fits':
        telluric_hLineCorrected = utils.nofits(telpath + '/' + hLinePrefix + os.path.basename(tel_spec))
        telluric_fitContinuum = utils.nofits(telpath + '/' + fitTelContinuumPrefix + os.path.basename(telluric_hLineCorrected))
        telluric_dividedContinuum = utils.nofits(telpath + '/' + dividedTelContinuumPrefix + os.path.basename(telluric_hLineCorrected))
        science_dividedTelluricLines = utils.nofits(telluricPrefix + os.path.basename(sci_spec))
        science_correctedTelluric = utils.nofits(dividedTelContinuumPrefix + os.path.basename(science_dividedTelluricLines))

        for valindex in range(startstep, stopstep + 1):
            logger.debug('valindex = %d', valindex)

            if valindex == 1:
                logger.info(" ------------------------------- ")
                logger.info("| STEP 1: Remove hydrogen lines |")
                logger.info(" ------------------------------- ")
                # Output: H line corrected Telluric 1D spectrum
                utils.pause(manualMode)

                hLineRemoval(tel_spec, telluric_hLineCorrected, hLineInter, orders, hLineMethod, hLineRegions,
                             vega_spec, tempInter, telpath + '/telluric_hLineInfo', overwrite)

            elif valindex == 2:
                logger.info(" -------------------------------- ")
                logger.info("| STEP 2: Fit Telluric continuum |")
                logger.info(" -------------------------------- ")
                # Output: 1D Fit to the telluric continuum
                utils.pause(manualMode)

                fitTelluricContinuum(telluric_hLineCorrected, telluric_fitContinuum, continuumInter, orders,
                                     continuumRegions, telluricFitOrders, tempInter, overwrite)

            elif valindex == 3:
                logger.info(" ------------------------------------------ ")
                logger.info("| STEP 3: Divide Telluric by its continuum |")
                logger.info(" ------------------------------------------ ")
                # Output: Continuum-divided, H line removed telluric 1D source spectra
                utils.pause(manualMode)

                divideTelluricContinuum(telluric_hLineCorrected, telluric_fitContinuum,
                                        telluric_dividedContinuum, orders, overwrite)

            elif valindex == 4:
                logger.info(" -------------------------------------------- ")
                logger.info("| STEP 4: Divide the science by the Telluric |")
                logger.info(" -------------------------------------------- ")
                # Output: Telluric line removed science 1D source spectra
                utils.pause(manualMode)

                telluricCorrection(sci_spec, telluric_dividedContinuum, science_dividedTelluricLines, telluricInter,
                                   orders, telluricRegions, 'science_telluricInfo.txt', overwrite)

            elif valindex == 5:
                logger.info(" -------------------------------------------- ")
                logger.info("| STEP 5: reintroduceTelluricContinuum       |")
                logger.info(" -------------------------------------------- ")
                # Output: Continuum-divided, telluric corrected science 1D source spectra.
                utils.pause(manualMode)

                # TODO: Give this step a better description

                reintroduceTelluricContinuum(science_dividedTelluricLines, telluric_fitContinuum,
                                             science_correctedTelluric, orders, overwrite)

        logger.info(" --------------------------------------------------- ")
        logger.info("| Telluric correction complete for this observation |")
        logger.info(" --------------------------------------------------- ")

    iraf.chdir(path)  # Return to the original directory

    return


# ----------------------------------------------------------------------------------------------------------------------
def hLineRemoval(infile, outroot, hLineInter, orders, hLineMethod, hLineRegions, vegafile, plot, inforoot, overwrite):
    """
    Remove hydrogen absorption lines from the extracted 1D Telluric spectrum.  Output filename prefix is 'h'.
    Reads:   Extracted 1D telluric spectrum.
    Writes:  H line corrected 1D telluric source spectrum.
    """
    logger = log.getLogger('hLineRemoval')

    utils.requires([infile])

    if utils.exists([outroot + '_order%d.fits' % i for i in range(1, len(orders) + 1)], overwrite):
        logger.info('Telluric standard already has the hydrogen lines removed.  Skipping this step.')
        return


    # Would it make more sense to move the loop over orders into each of the functions?
    # That would pretty much eliminate the need for this function (hlineremoval).

    airmass = fits.getheader(infile)['AIRMASS']

    for i in range(len(orders)):
        logger.debug('Processing order %d...', orders[i])
        extension = i+1

        inspec = infile + '[SCI,' + str(extension) + ']'
        outspec = outroot + '_order%d.fits' % orders[i]
        vegaspec = vegafile + '[' + str(extension) + ']'
        infofile = inforoot + '_order%d.dat' % orders[i]

        logger.debug('Adding the AIRMASS keyword...')
        iraf.hedit(images=inspec, fields='AIRMASS', value=airmass, add='yes', addonly='no', delete='no',
                   verify='no', show='no', update='yes')

        if hLineMethod == 'None':
            logger.warning("Skipping H line removal for order %d.", orders[i])
            logger.info("Copying files to have the right names for later use in the pipeline.")
            iraf.imcopy(input=inspec, output=outspec, verbose='yes')

        elif 'vega' in hLineMethod.lower():  # 'vega' or 'Vega_Tweak' or 'Tweak-Vega'

            logger.info("Removing H lines from order %d...", orders[i])
            vega(inspec, vegaspec, outspec, hLineInter, hLineRegions[orders[i]], infofile, overwrite)






        # if 'tweak' in hLineMethod.lower():  # Vega-Tweak or LineFit-Tweak





        # TODO(Viraja):  Update the functions for the rest of the H line removal options (commented out below).
        # NOTE: Untested because interactive scripted iraf tasks are broken... Once ready, add them to the
        # part of the script where the output did not already exists and the task had to be run to generate
        # files for the first time.

        elif hLineMethod == 'lineFitAuto':
            logger.error('Auto line fitting not yet implemented')
            # lineFitAuto(combined_extracted_1d_spectra, grating)

        elif hLineMethod == 'lineFitManual':
            logger.error('Manual line fitting not yet implemented')
            # lineFitManual(combined_extracted_1d_spectra+'[sci,1]', grating)

        elif hLineMethod == 'vega_tweak':
            logger.error('Vega tweak not yet implemented')
            # First, do H line removal using the 'vega' method automatically, and then give user a chance
            # to interact with the spectrum
            # vega(combined_extracted_1d_spectra, path, hLineInter, telluric_shift_scale_record, overwrite)
            # lineFitManual("final_tel_no_hLines_no_norm", grating)

        elif hLineMethod == 'lineFit_tweak':
            logger.error('Line fit tweaking not yet implemented')
            # First, do H line removal using the 'lineFitAuto' method automatically, then give user a
            # chance to interact with the spectrum
            # lineFitAuto(combined_extracted_1d_spectra,grating)
            # lineFitManual("final_tel_no_hLines_no_norm", grating)

        else:
            logger.error("Unrecognized H line removal method: %s", hLineMethod)
            raise SystemExit




        # Create a MEF setting the original primary header as extension [0]
        # iraf.wmef(input=telluric_hLineRemovalOutput_SEF, output=telluric_hLineRemovalOutput_MEF, extnames='',
        #    phu=tel_src_extracted_spectrum, verbose='yes', mode='al')
        
        # TODO(Viraja):  Temporarily, pause to display and take a look at the telluric 1D source spectra before and 
        # after H line correction. I am not sure if this can be more useful than just looking at the spectra. If 
        # this has to be done, it would require adding another science extension to the MEF created in the previous 
        # step using iraf.wmef task BEFORE starting the H line correction in this loop.

        if plot:
            # Plot the telluric 1D source spectrum with and without H line correction
            uncorrected = fits.open(inspec).data
            corrected = fits.open(outspec).data
            plt.title('Telluric 1D Spectrum Before and After H-Line Removal')
            plt.plot(uncorrected, 'k', label='Telluric w/o H-line correction')
            plt.plot(corrected, 'r', label='Telluric w/ H-line correction')
            plt.legend(loc='best')
            plt.show()

    return


# ----------------------------------------------------------------------------------------------------------------------
def vega(inputtelspectrum, inputcalspectrum, outputtelspectrum, interactive, sample, infofile, overwrite):
    # Use IRAF TELLURIC to remove H absorption lines by scaling and shifting a Vega spectrum.
    logger = log.getLogger('vega')
    logger.debug('inputtelspectrum: %s', inputtelspectrum)
    logger.debug('inputcalspectrum: %s', inputcalspectrum)
    logger.debug('sample: %s', sample)
    logger.debug('infofile: %s', infofile)

    if utils.exists([outputtelspectrum, infofile], overwrite):
        logger.warning('Corrected spectrum already exists.  Skipping this step.')
        return

    results = iraf.telluric(
        input=inputtelspectrum,
        output=outputtelspectrum,
        cal=inputcalspectrum,
        ignoreaps='yes', xcorr='yes', tweakrms='yes', interactive=interactive, sample=sample, threshold=0.1, lag=3,
        shift=0., dshift=0.05, scale=1., dscale=0.05, offset=0, smooth=1, cursor='', answer='yes', mode='al', Stdout=1)

    logger.debug('Results: %s', results)

    # What is this used for?
    # Record shift and scale results for future reference in the script
    # NOTE: If function "hLineremoval" is run more than once, the function call to "vega" will overwrite the old file.
    with open(infofile, 'w') as f:
        f.write(str(results) + '\n')

    # Parse the Telluric output.  The normalization is usually last [-1], unless there is a warning about pixels
    # outside limits at the end, then the normlization is second to last [-2]:
    if 'limits' in results[-1].split()[-1]:
        normalization = results[-2].split()[-1]
    else:
        normalization = results[-1].split()[-1]
    logger.info('Normalization: %s', normalization)
    
    logger.debug('Undoing the normalization introduced by iraf.telluric...')
    iraf.imarith(
        operand1=outputtelspectrum, operand2=normalization, op='/', result=outputtelspectrum, title='',
        divzero=0.0, hparams='', pixtype='', calctype='', verbose='yes', noact='no', mode='al')

    '''
    # Comment from nifsTelluric.py -- There are subtle bugs in iraf mean imarith does not work. So we use an 
    # astropy/numpy solution. Open the image and the scalar we will be dividing it by.
    # Viraja:  I am not sure what "... subtle bugs in iraf mean imarith does not work" means. The comment is not clear.
    # TODO(Viraja):  This commented section should be incorporated later in all relevant places of this script if the 
    # imarith way of removing normalization adopted here from XDGNIRS is found to not work properly during further 
    # testing.  As of August 2019, iraf.imarith seems to work fine.
    operand1 = fits.open(outputtelspectrum)[0].data
    operand2 = float(normalization)
    # Create a new data array
    if operand2 != 0:
        operand1 = operand1 / operand2
    else:
        operand1 = 1
        # Viraja:  Why is operand1 set to 1? I would imagine it is set to itself as the divisor is 0.
    '''

    return


# ----------------------------------------------------------------------------------------------------------------------
def fitTelluricContinuum(telluric_hLineCorrectedSpectrum, telluric_fitContinuum, continuumInter, orders,
                         continuumRegions, fitTelContinuumOrders, tempInter, overwrite):
    """
    Fit the continua of the H line corrected telluric 1D source spectra using IRAF CONTINUUM task to normalize them.
    """
    logger = log.getLogger('fitTelluricContinuum')

    sample = {}
    for order, s in continuumRegions:
        sample[int(order)] = s

    # Continuum fitting order depends on the spectral order.
    for i in range(len(orders)):
        extension = i+1
        telluric_fitContinuumInput = telluric_hLineCorrectedSpectrum+'_order'+str(orders[i])+'_MEF.fits[1]'
        telluric_fitContinuumOutput = telluric_fitContinuum+'_order'+str(orders[i])+'.fits'
    
        if os.path.exists(telluric_fitContinuumOutput):
            if overwrite:
                logger.warning("Removing old %s", telluric_fitContinuumOutput)
                os.remove(telluric_fitContinuumOutput)
            else:
                logger.warning("Output exists and -overwrite not set - skipping telluric continuum fitting for order")
                logger.warning("%d of the H line corrected spectra.", orders[i])
                continue

        # Here, iraf.continuum is used with type='fit' to write out the continuum fit rather than type='ratio'.

        # NOTE:  As of August 2019, the help documentation for iraf.continuum does not show 'ask' as a
        # non-optional input parameter; however, <lpar continuum> lists ask as the third non-optional paranter
        # to be given with the task. If the input is a list of spectra as against an individual image and 'ask'
        # is set to 'yes', the user will be asked, while proceeding to the next spectrum in the input list,
        # whether they would like to run the task interactively. A "YES" would run the task interactively for
        # all spectra in the input list.

        # NOTE: In TelluricXD.py in XDGNIRS, parameter 'logfiles', where to write the power series coefficients,
        # was NULL (""). In nifsTelluric.py in NIFTY, parameter 'logfiles' was set to the main logging file for
        # the data reduction.

        logger.info("Fitting telluric continuum for order %d of the H line corrected spectra.", orders[i])
        iraf.continuum(input=telluric_fitContinuumInput, output=telluric_fitContinuumOutput, ask='yes', lines='*',
            bands='1', type="fit", replace='no', wavescale='yes', logscale='no', override='no', listonly='no',
            logfiles=logger.root.handlers[0].baseFilename, inter=continuumInter, sample=sample[orders[i]], 
            naverage=1, func='spline3', order=fitTelContinuumOrders[i], low_reject=1.0, high_reject=3.0, niterate=2, 
            grow=1.0, markrej='yes', graphics='stdgraph', cursor='', mode='ql')

        if tempInter:
            # Plot H line corrected telluric 1D source spectrum without continuum fitting and the continuum fit
            telluric_withoutNormalization = fits.open(telluric_hLineCorrectedSpectrum+'_order'+str(orders[i])+'_MEF.fits')[1].data
            telluric_continuumFit = fits.open(telluric_fitContinuumOutput)[0].data
            plt.title('H-Line Corrected Telluric 1D Source Spectrum and Continuum Fit for Normalization')
            plt.plot(telluric_withoutNormalization, 'k', label='Telluric w/o normalization')
            plt.plot(telluric_continuumFit, 'r', label='Telluric continuum fit')
            plt.legend(loc='best')
            plt.show()


# ----------------------------------------------------------------------------------------------------------------------
def divideTelluricContinuum(telluric_hLineCorrectedSpectrum, telluric_fitContinuum, telluric_dividedContinuum, orders,
    overwrite):
    """
    Divide the H line corrected telluric 1D source spectra by its continuum fit to normalize it.
    """
    logger = log.getLogger('divideTelluricContinuum')

    for i in range(len(orders)):
        extension = i+1
        telluric_divideContinuumInput = telluric_hLineCorrectedSpectrum+'_order'+str(orders[i])+'_MEF.fits[1]'
        telluricContinuumFit = telluric_fitContinuum+'_order'+str(orders[i])+'.fits'
        telluric_divideContinuumOutput_SEF = telluric_dividedContinuum+'_order'+str(orders[i])+'_SEF.fits'
        telluric_divideContinuumOutput_MEF = telluric_dividedContinuum+'_order'+str(orders[i])+'_MEF.fits'

        oldfiles = glob.glob(telluric_dividedContinuum+'_order'+str(orders[i])+'*.fits')
        if len(oldfiles) > 0:
            if overwrite:
                logger.warning('Removing old %s_order'+str(orders[i])+'*.fits files.', telluric_dividedContinuum)
                [os.remove(filename) for filename in oldfiles]
            else:
                logger.warning("Output exists and -overwrite not set - skipping division for order %d of", orders[i])
                logger.warning("the H line corrected telluric 1D source spectra by the telluric continuum fit.")
                continue

        # Divide the H line corrected telluric 1D source spectra by the continuum fit.
        logger.info("Dividing order %d of the H line corrected telluric 1D source spectra by the telluric", orders[i])
        logger.info("continuum fit.")

        iraf.imarith(operand1=telluric_divideContinuumInput, operand2=telluricContinuumFit, op='/',
            result=telluric_divideContinuumOutput_SEF, title='', divzero=0.0, hparams='', pixtype='', calctype='',
            verbose='yes', noact='no', mode='al')

        # Creating a MEF setting the original primary header as extension [0]
        iraf.wmef(input=telluric_divideContinuumOutput_SEF, output=telluric_divideContinuumOutput_MEF,
            extnames='', phu=telluric_hLineCorrectedSpectrum+'_order'+str(orders[i])+'_MEF', verbose='yes',
            mode='al')


# ----------------------------------------------------------------------------------------------------------------------
def telluricCorrection(sci_src_extracted_spectrum, telluric_dividedContinuum, science_dividedTelluricLines,
    telluricInter, orders, telluricRegions, sci_tel_infofile, overwrite):
    """
    Run IRAF telluric task on each order
    :param sci_src_extracted_spectrum: input multi-extension fits file name
    :param telluric_dividedContinuum: input calibration spectrum file name
    :param science_dividedTelluricLines: output spectrum file name
    :param telluricInter: interactive parameter setting for iraf.telluric
    :param orders: list of orders
    :param sci_airmass: float airmass
    :param telluricRegions: list of regions to fit in each order [(order3,region1), (order4,region2), ...]
    :param sci_tel_infofile: filename for output results that are normally printed to the terminal
    :param overwrite: overwrite parameter setting
    """
    logger = log.getLogger('telluricCorrection')

    sci_airmass = fits.getheader(sci_src_extracted_spectrum)['AIRMASS']

    sample = {}
    for order, s in telluricRegions:
        sample[int(order)] = s

    for i in range(len(orders)):
        extension = i+1
        science_telluricCorrectionInput = os.path.basename(sci_src_extracted_spectrum)+'[SCI,'+str(extension)+']'
        calibrationInput = telluric_dividedContinuum+'_order'+str(orders[i])+'_MEF.fits[1]'
        telluric_telluricCorrectionOutput_SEF = science_dividedTelluricLines+'_order'+str(orders[i])+'_SEF.fits'
        telluric_telluricCorrectionOutput_MEF = science_dividedTelluricLines+'_order'+str(orders[i])+'_MEF.fits'
        
        oldfiles = glob.glob(science_dividedTelluricLines+'_order'+str(orders[i])+'*.fits')
        if len(oldfiles) > 0:
            if overwrite:
                logger.warning('Removing %s_order'+str(orders[i])+'*.fits files.', science_dividedTelluricLines)
                [os.remove(filename) for filename in oldfiles]
            else:
                logger.warning("Output exists and overwrite flag not set.")
                logger.warning("Skipping division of order %d of the telluric corrected science 1D source", orders[i])
                logger.warning("spectra by order %d of the H line corrected, continuum-divided telluric 1D", orders[i])
                logger.warning("source spectra.")
                continue

        science_telluricInfo = iraf.telluric(input=science_telluricCorrectionInput, 
            output=telluric_telluricCorrectionOutput_SEF, cal=calibrationInput, airmass=sci_airmass, answer='yes',
            ignoreaps='yes', xcorr='yes', tweakrms='yes', interactive=telluricInter, sample=sample[orders[i]], 
            threshold=0.1, lag=10, shift=0.0, dshift=1.0, scale=1.0, dscale=0.2, offset=1, smooth=1, cursor='', 
            Stdout=1)

        with open(sci_tel_infofile,'w') as f:
            f.write(str(science_telluricInfo)+'\n')

        if "limits" in science_telluricInfo[-1].split()[-1]:
            index = -2
        else:
            index = -1
        normalization = science_telluricInfo[index].split()[-1]
        scale = float(science_telluricInfo[index].split()[-4].replace(',', ''))
        shift = float(science_telluricInfo[index].split()[-7].replace(',', ''))
        logger.info('Shift: %.2f pixels', shift)
        if abs(float(shift) > 1.0):
            logger.warning("Shift shift > 1 pixel!")

        logger.info("Scale: %.3f", scale)
        if abs(scale - 1.0) > 0.1:
            logger.warning("Scaling is > 10%")

        # Remove normalization introduced by iraf.telluric
        iraf.imarith(operand1=telluric_telluricCorrectionOutput_SEF, operand2=normalization, op='/',
            result=telluric_telluricCorrectionOutput_SEF, title='', divzero=0.0, hparams='', pixtype='',
            calctype='', verbose='yes', noact='no', mode='al')

        # Creating a MEF setting the original primary header as extension [0]
        iraf.wmef(input=telluric_telluricCorrectionOutput_SEF, output=telluric_telluricCorrectionOutput_MEF,
            extnames='', phu=sci_src_extracted_spectrum, verbose='yes', mode='al')

    return


# ----------------------------------------------------------------------------------------------------------------------
def reintroduceTelluricContinuum(science_dividedTelluricLines, telluric_fitContinuum, science_correctedTelluric, \
    orders, overwrite):
    """
    Re-introduce telluric continuum shape into telluric lines removed science 1D source spectra.
    """
    logger = log.getLogger('reintroduceTelluricContinuum')

    # Start with regularly extracted science spectra (not using full slit or stepwise extraction methods)
    for i in range(len(orders)):
        extension = i+1
        science_reintroduceTelluricContinuumInput = science_dividedTelluricLines+'_order'+str(orders[i])+'_MEF.fits[1]'
        telluricContinuumFit = telluric_fitContinuum+'_order'+str(orders[i])+'.fits'
        science_reintroduceTelluricContinuumOutput_SEF = science_correctedTelluric+'_order'+str(orders[i])+'_SEF.fits'
        science_reintroduceTelluricContinuumOutput_MEF = science_correctedTelluric+'_order'+str(orders[i])+'_MEF.fits'
        
        oldfiles = glob.glob(science_correctedTelluric+'_order'+str(orders[i])+'*.fits')
        if len(oldfiles) > 0:
            if overwrite:
                logger.warning('Removing %s_order'+str(orders[i])+'*.fits files.', science_correctedTelluric)
                [os.remove(f) for f in oldfiles]
            else:
                logger.warning("Output exists and -overwrite not set - skipping re-introduction of telluric continuum")
                logger.warning("into order %d of the telluric corrected removed science 1D source spectra.", orders[i])
                continue

        logger.info("Re-introducing telluric continuum into order %d of the telluric lines removed science", orders[i])
        logger.info("1D source spectra.\n")
        iraf.imarith(operand1=science_reintroduceTelluricContinuumInput, operand2=telluricContinuumFit, op="/",
            result=science_reintroduceTelluricContinuumOutput_SEF, title='', divzero=0.0, hparams='', pixtype='',
            calctype='', verbose='yes', noact='no', mode='al')

        # Creating a MEF setting the original primary header as extension [0]
        iraf.wmef(input=science_reintroduceTelluricContinuumOutput_SEF, output=science_reintroduceTelluricContinuumOutput_MEF, extnames='',
            phu=science_dividedTelluricLines+'_order'+str(orders[i])+'_MEF', verbose='yes', mode='al')

# ----------------------------------------------------------------------------------------------------------------------
'''
# TODO(nat): linefitAuto and linefitManual could be useful at some point.
def lineFitAuto(spectrum, grating):
    """
    Automatically fit Lorentz profiles to lines defined in existing cur* files. Go to x position in cursor file and use 
    space bar to find spectrum at each of those points.
    """
    logger = log.getLogger('lineFitAuto')

    specpos = iraf.bplot(images=spectrum+'[SCI,1]', cursor='cur'+grating, Stdout=1, StdoutG='/dev/null')
    specpose = str(specpos).split("'x,y,z(x):")
    nextcur = 'nextcur'+grating+'.txt'
    # Write line x,y info to file containing Lorentz fitting commands for bplot
    write_line_positions(nextcur, specpos)
    iraf.delete('final_tel_no_hLines_no_norm.fits',ver="no",go_ahead='yes',Stderr='/dev/null')
    # Fit and subtract Lorentz profiles. Might as well write output to file.
    iraf.bplot(images=spectrum+'[sci,1]',cursor='nextcur'+grating+'.txt', new_image='final_tel_no_hLines_no_norm', overwrite="yes",StdoutG='/dev/null',Stdout='Lorentz'+grating)

# ----------------------------------------------------------------------------------------------------------------------

def lineFitManual(spectrum, grating):
    """ 
    Enter splot so the user can fit and subtract lorentz (or, rather any) profiles.
    """
    logger = log.getLogger('lineFitManual')

    iraf.splot(images=spectrum, new_image='final_tel_no_hLines_no_norm', save_file='../PRODUCTS/lorentz_hLines.txt', overwrite='yes')
    # it's easy to forget to use the 'i' key to actually write out the line-free spectrum, so check that it exists:
    # with the 'tweak' options, the line-free spectrum will already exists, so this lets the user simply 'q' and move on w/o editing (too bad if they edit and forget to hit 'i'...)
    while True:
        try:
            with open("final_tel_no_hLines_no_norm.fits") as f: pass
            break
        except IOError as e:
            logger.info("It looks as if you didn't use the i key to write out the lineless spectrum. We'll have to try again. --> Re-entering splot")
            iraf.splot(images=spectrum, new_image='final_tel_no_hLines_no_norm', save_file='../PRODUCTS/lorentz_hLines.txt', overwrite='yes')
'''
# ----------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs-pype.cfg')
