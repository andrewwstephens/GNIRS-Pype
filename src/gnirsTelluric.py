#!/usr/bin/env python

from astropy.io import fits
import ConfigParser
import glob
import log
import matplotlib.pyplot as plt
import os
from pyraf import iraf


# ----------------------------------------------------------------------------------------------------------------------
def start(configfile):
    """
    Do a Telluric correction using the IRAF TELLURIC task.
    """
    logger = log.getLogger('gnirsTelluric.start')

    path = os.getcwd()  # Store current working directory for later use

    logger.info('#################################################')
    logger.info('#                                               #')
    logger.info('#       Start GNIRS Telluric Correction         #')
    logger.info('#                                               #')
    logger.info('#################################################\n')

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
    
    # telluricCorrection specific config
    start = config.getint('telluricCorrection', 'Start')
    stop = config.getint('telluricCorrection', 'Stop')
    hLineMethod = config.get('telluricCorrection', 'hLineMethod')
    hLineRegions = config.items('hLineRegions')
    continuumRegions = config.items('continuumRegions')
    telluricRegions = config.items('telluricRegions')

    for scipath in config.options("ScienceDirectories"):

        if not config.getboolean("ScienceDirectories", scipath):  # only process directories that are marked True
            logger.debug('Skipping %s', scipath)
            continue

        ###########################################################################
        ##                                                                       ##
        ##                  BEGIN - OBSERVATION SPECIFIC SETUP                   ##
        ##                                                                       ##
        ###########################################################################

        scipath += '/Intermediate'
        logger.info("Moving to science directory: %s", scipath)
        iraf.chdir(scipath)

        # Get symbolic paths to the std and tel directories in the sci directory and the runtime data directory
        # Relative path/link expected to be at the top level of every sci directory
        telpath = '../Telluric/Intermediate'
        logger.info("Telluric directory: %s\n", telpath)
        logger.info("Runtime data path: %s\n", runtimedata)
        
        # Check if required runtime data files and extracted science source spectra available in their respective paths
        logger.info("Checking if required runtime data files and extracted spectra available in %s", runtimedata)
        logger.info("and %s, respectively.\n", scipath)

        vega_spectrum = runtimedata + 'vega_ext.fits'
        if os.path.exists(vega_spectrum):
            logger.info("Required vega spectrum available.")
        else:
            logger.warning("Required vega spectrum not available. Please provide the vega spectrum in")
            logger.warning("%s", runtimedata)
            logger.warning("Exiting script.\n")
            raise SystemExit

        sci_src_extracted_spectrum = scipath + '/' + extractRegularPrefix + combinedsrc
        if os.path.exists(sci_src_extracted_spectrum):
            logger.info("Required extracted science source spectrum available.")
            sci_airmass = fits.open(sci_src_extracted_spectrum)[0].header['AIRMASS']
        else:
            logger.warning("Required extracted science source spectrum not available.")
            logger.warning("Please run gnirsExtarctSpectra1D.py to create the extracted spectrum or provide it")
            logger.warning("manually in %s", scipath)
            logger.warning("Exiting script.\n")
            raise SystemExit

        tel_src_extracted_spectrum = telpath + '/' + extractRegularPrefix + combinedsrc
        if os.path.exists(tel_src_extracted_spectrum):
            logger.info("Required extracted telluric 1D source spectrum available.")
            tel_airmass = fits.open(tel_src_extracted_spectrum)[0].header['AIRMASS']
        else:
            logger.warning("Required extracted telluric 1D source spectrum not available.")
            logger.warning("Please run gnirsExtarctSpectra1D.py to create the extracted spectrum or provide it")
            logger.warning("manually in %s", telpath)
            logger.warning("Exiting script.\n")
            raise SystemExit

        # TODO:  if extractionStepwise:

        logger.info("Required extracted spectra check complete.\n")

        # Record the number of orders expected according to the GNIRS XD configuration.
        if 'LB_SXD' in scipath:
            orders = [3, 4, 5]
            fitTelContinuumOrders = [config.getint('telluricCorrection','fitTelContinuum_order3'), \
                config.getint('telluricCorrection','fitTelContinuum_order4'), \
                config.getint('telluricCorrection','fitTelContinuum_order5')]
        elif 'LB_LXD' in scipath:
            orders = [3, 4, 5, 6, 7, 8]
            fitTelContinuumOrders = [config.getint('telluricCorrection','fitTelContinuum_order3'), \
                config.getint('telluricCorrection','fitTelContinuum_order4'), \
                config.getint('telluricCorrection','fitTelContinuum_order5'), \
                config.getint('telluricCorrection','fitTelContinuum_order6'), \
                config.getint('telluricCorrection','fitTelContinuum_order7'), \
                config.getint('telluricCorrection','fitTelContinuum_order8')]
        elif 'SB_SXD' in scipath:
            orders = [3, 4, 5, 6, 7, 8]
            fitTelContinuumOrders = [config.getint('telluricCorrection','fitTelContinuum_order3'), \
                config.getint('telluricCorrection','fitTelContinuum_order4'), \
                config.getint('telluricCorrection','fitTelContinuum_order5'), \
                config.getint('telluricCorrection','fitTelContinuum_order6'), \
                config.getint('telluricCorrection','fitTelContinuum_order7'), \
                config.getint('telluricCorrection','fitTelContinuum_order8')]
        else:
            logger.error("#############################################################################")
            logger.error("#                                                                           #")
            logger.error("#     ERROR in telluric: unknown GNIRS XD configuration. Exiting script.    #")
            logger.error("#                                                                           #")
            logger.error("#############################################################################\n")
            raise SystemExit

        # Define output filenames with prefixes added at different stages of script but without the '.fits':
        telluric_hLineCorrectedSpectrum = nofits(telpath + '/' + hLinePrefix + \
                os.path.basename(tel_src_extracted_spectrum))
        telluric_fitContinuum = nofits(telpath + '/' + fitTelContinuumPrefix + \
            os.path.basename(telluric_hLineCorrectedSpectrum))
        telluric_dividedContinuum = nofits(telpath + '/' + dividedTelContinuumPrefix + \
            os.path.basename(telluric_hLineCorrectedSpectrum))
        science_dividedTelluricLines = nofits(telluricPrefix + os.path.basename(sci_src_extracted_spectrum))
        science_correctedTelluric = nofits(dividedTelContinuumPrefix + os.path.basename(science_dividedTelluricLines))

        #########################################################################
        #                                                                       #
        #                 COMPLETE - OBSERVATION SPECIFIC SETUP                 #
        #             BEGIN TELLURIC CORRECTION FOR AN OBSERVATION              #
        #                                                                       #
        #########################################################################

        valindex = start
        while valindex > stop or valindex < 1 or stop > 5:
            logger.warning("#####################################################################")
            logger.warning("#                                                                   #")
            logger.warning("#   WARNING in telluric: invalid start/stop values of Telluric      #")
            logger.warning("#                        correction steps.                          #")
            logger.warning("#                                                                   #")
            logger.warning("#####################################################################\n")
            valindex = int(raw_input("Please enter a valid start value (1 to 5, default 1): "))
            stop = int(raw_input("Please enter a valid stop value (1 to 5, default 5): "))

        while valindex <= stop:

            ###########################################################################
            #  STEP 1: H line removal.                                                #
            #  Output: H line corrected telluric 1D source spectra.                   #
            ###########################################################################

            if valindex == 1:
                if manualMode:
                    a = raw_input("About to enter step 1: H line removal.")

                hLineRemoval(tel_src_extracted_spectrum, telluric_hLineCorrectedSpectrum, hLineInter, orders,
                    hLineMethod, hLineRegions, tel_airmass, vega_spectrum, tempInter,
                    telpath+'/telluric_hLineInfo.txt', overwrite)

                logger.info("##################################################################")
                logger.info("#                                                                #")
                logger.info("#       STEP 1: H line removal - COMPLETED                       #")
                logger.info("#                                                                #")
                logger.info("##################################################################\n")

            ###########################################################################
            #  STEP 2: Fit telluric continuum.                                        #
            #  Output: 1D Fit to the telluric continuum.                              #
            ###########################################################################

            elif valindex == 2:
                if manualMode:
                    a = raw_input("About to enter step 2: Fit telluric continuum.")

                fitTelluricContinuum(telluric_hLineCorrectedSpectrum, telluric_fitContinuum, continuumInter, orders,
                    continuumRegions, fitTelContinuumOrders, tempInter, overwrite)

                logger.info("##################################################################")
                logger.info("#                                                                #")
                logger.info("#       STEP 2: Fit telluric continuum - COMPLETED               #")
                logger.info("#                                                                #")
                logger.info("##################################################################\n")

            ###########################################################################################
            #  STEP 3: Division of telluric by the telluric continuum.                                #
            #  Output: Continuum-divided, H line removed telluric 1D source spectra.                  #
            ###########################################################################################

            elif valindex == 3:
                if manualMode:
                    a = raw_input("About to enter step 3: Division of telluric by the telluric continuum.")

                divideTelluricContinuum(telluric_hLineCorrectedSpectrum, telluric_fitContinuum, \
                    telluric_dividedContinuum, orders, overwrite)

                logger.info("##############################################################################")
                logger.info("#                                                                            #")
                logger.info("#    STEP 3: Division of telluric by the telluric continuum  - COMPLETED     #")
                logger.info("#                                                                            #")
                logger.info("##############################################################################\n")

            ###########################################################################
            #  STEP 4: Telluric line removal.                                         #
            #  Output: Telluric line removed science 1D source spectra.               #
            ###########################################################################

            elif valindex == 4:
                if manualMode:
                    a = raw_input("About to enter step 4: Telluric line removal.")

                telluricCorrection(sci_src_extracted_spectrum, telluric_dividedContinuum, science_dividedTelluricLines, \
                    telluricInter, orders, sci_airmass, telluricRegions, 'science_telluricInfo.txt', overwrite)

                logger.info("##################################################################")
                logger.info("#                                                                #")
                logger.info("#       STEP 4: Telluric line removal - COMPLETED                #")
                logger.info("#                                                                #")
                logger.info("##################################################################\n")

            ###############################################################################
            #  STEP 5: Division of science by the telluric continuum.                     #
            #  Output: Continuum-divided, telluric corrected science 1D source spectra.   #
            ###############################################################################

            elif valindex == 5:
                if manualMode:
                    a = raw_input("About to enter step 5: Division of science by the telluric continuum.")

                reintroduceTelluricContinuum(science_dividedTelluricLines, telluric_fitContinuum, \
                    science_correctedTelluric, orders, overwrite)

                logger.info("##################################################################################")
                logger.info("#                                                                                #")
                logger.info("#       STEP 5: Division of telluric by the telluric continuum - COMPLETED       #")
                logger.info("#                                                                                #")
                logger.info("##################################################################################\n")

            else:
                logger.error("###########################################################################")
                logger.error("#                                                                         #")
                logger.error("#      ERROR in telluric: %d is not valid. Exiting script.", valindex)
                logger.error("#                                                                         #")
                logger.error("###########################################################################\n")
                raise SystemExit

            valindex += 1

        logger.info("##############################################################################")
        logger.info("#                                                                            #")
        logger.info("#  COMPLETE - Telluric correction completed for                              #")
        logger.info("#  %s", scipath)
        logger.info("#                                                                            #")
        logger.info("##############################################################################")


    iraf.chdir(path)  # Return to the original directory

    return


# ----------------------------------------------------------------------------------------------------------------------
def hLineRemoval(tel_src_extracted_spectrum, telluric_hLineCorrectedSpectrum, hLineInter, orders, hLineMethod,
    hLineRegions, tel_airmass, vega_spectrum, tempInter, tel_hline_infofile, overwrite):
    """
    Remove hydrogen (H) absorption lines from the extracted 1D telluric source spectrum. Output filename prefix is 'h'.

    Reads:
        Extracted 1D telluric source spectrum.

    Writes:
        H line corrected 1D telluric source spectrum.
    """
    logger = log.getLogger('gnirsTelluric.hLineRemoval')

    sample = {}
    for order, s in hLineRegions:
        sample[int(order)] = s

    for i in range(len(orders)):

        extension = i+1
        telluric_hLineRemovalInput = tel_src_extracted_spectrum+'[SCI,'+str(extension)+']'
        calibrationInput = vega_spectrum + '['+str(extension)+']'
        telluric_hLineRemovalOutput_SEF = telluric_hLineCorrectedSpectrum+'_order'+str(orders[i])+'_SEF.fits'
        telluric_hLineRemovalOutput_MEF = telluric_hLineCorrectedSpectrum+'_order'+str(orders[i])+'_MEF.fits'

        oldfiles = glob.glob(telluric_hLineCorrectedSpectrum+'_order'+str(orders[i])+'*.fits')
        if len(oldfiles) > 1:
            if overwrite:
                logger.warning('Removing old %s_order'+str(orders[i])+'*.fits files.', telluric_hLineCorrectedSpectrum)
                [os.remove(filename) for filename in oldfiles]
            else:
                logger.warning("Output exists and -overwrite not set - skipping H line removal for order %d", orders[i])
                logger.warning("of the telluric 1D source spectra.")
                continue

        iraf.hedit(images=telluric_hLineRemovalInput, fields='AIRMASS', value=tel_airmass, add='yes',
            addonly='no', delete='no', verify='no', show='no', update='yes', mode='al')

        if not hLineMethod:
            logger.warning("Parameter 'hLineMethod' in the configuration file is %s.", hLineMethod)
            logger.warning("Skipping H line removal for order %d of the telluric 1D source spectra.", orders[i])
            logger.info("Copying files to have the right names for later use in the pipeline.")
            iraf.imcopy(input=telluric_hLineRemovalInput, output=telluric_hLineRemovalOutput_SEF, verbose='yes')
        else:
            logger.info("Removing H lines from order %d of the telluric 1D source spectra.", orders[i])
            if hLineMethod == 'vega':
                vega(telluric_hLineRemovalInput, calibrationInput, telluric_hLineRemovalOutput_SEF, hLineInter,
                    tel_airmass, sample[orders[i]], tel_hline_infofile, overwrite)

            # TODO(Viraja):  Update the functions for the rest of the H line removal options (commented out below).
            # NOTE: Untested because interactive scripted iraf tasks are broken... Once ready, add them to the 
            # part of the script where the output did not already exists and the task had to be run to generate
            # files for the first time. 
            elif hLineMethod == 'lineFitAuto':
#                    lineFitAuto(combined_extracted_1d_spectra, grating)
                pass

            elif hLineMethod == 'lineFitManual':
#                    lineFitManual(combined_extracted_1d_spectra+'[sci,1]', grating)
                pass

            elif hLineMethod == 'vega_tweak':
                # First, do H line removal using the 'vega' method automatically, and then give user a chance  
                # to interact with the spectrum
#                    vega(combined_extracted_1d_spectra, path, hLineInter, telluric_shift_scale_record, overwrite)
#                    lineFitManual("final_tel_no_hLines_no_norm", grating)
                pass

            elif hLineMethod == 'lineFit_tweak':
                # First, do H line removal using the 'lineFitAuto' method automatically, then give user a 
                # chance to interact with the spectrum
#                    lineFitAuto(combined_extracted_1d_spectra,grating)
#                    lineFitManual("final_tel_no_hLines_no_norm", grating)
                pass
            
            else:
                logger.error("Unrecognized H line removal method encountered in the configuration file.")
                logger.error("Exiting script.\n")
                raise SystemExit

        # Creating a MEF setting the original primary header as extension [0]
        iraf.wmef(input=telluric_hLineRemovalOutput_SEF, output=telluric_hLineRemovalOutput_MEF, extnames='',
            phu=tel_src_extracted_spectrum, verbose='yes', mode='al')
        
        # TODO(Viraja):  Temporarily, pause to display and take a look at the telluric 1D source spectra before and 
        # after H line correction. I am not sure if this can be more useful than just looking at the spectra. If 
        # this has to be done, it would require adding another science extension to the MEF created in the previous 
        # step using iraf.wmef task BEFORE starting the H line correction in this loop.
        if tempInter:
            # Plot the telluric 1D source spectrum with and without H line correction
            telluric_hLineUncorrected = fits.open(tel_src_extracted_spectrum)[i*3+1].data
            telluric_hLineCorrected = fits.open(telluric_hLineRemovalOutput_MEF)[1].data
            plt.title('Telluric 1D Source Spectrum Before and After H-Line Removal')
            plt.plot(telluric_hLineUncorrected, 'k', label='Telluric w/o H-line correction')
            plt.plot(telluric_hLineCorrected, 'r', label='Telluric w/ H-line correction')
            plt.legend(loc='best')
            plt.show()

#---------------------------------------------------------------------------------------------------------------------#

def vega(inputtelspectrum, inputcalspectrum, outputtelspectrum, hLineInter, tel_airmass, sample, tel_hline_infofile, \
    overwrite):
    """
    Use IRAF TELLURIC task to remove H absorption lines from the telluric 1D source spectrum using the 'vega' method, 
    then divide normalization introduced by TELLURIC with IMARITH.
    """
    logger = log.getLogger('gnirsTelluric.vega')

    # For each order, this will be done interactively if parameter 'hLineInter' in the configuration file is 'yes'
    telluric_hLineInfo = iraf.telluric(input=inputtelspectrum, output=outputtelspectrum, cal=inputcalspectrum,
        ignoreaps='yes', xcorr='yes', tweakrms='yes', interactive=hLineInter, sample=sample, threshold=0.1,
        lag=3, shift=0., dshift=0.05, scale=1., dscale=0.05, offset=0, smooth=1, cursor='', airmass=tel_airmass,
        answer='yes', mode='al', Stdout=1)

    # Record shift and scale info for future reference in the script
    # NOTE: If function "hLineremoval" is ran more than once, the function call to "vega" will automatically overwrite 
    # the existing tel_hline_infofile.
    tel_hline_info = open(tel_hline_infofile,'w')
    tel_hline_info.write(str(telluric_hLineInfo)+'\n')
    tel_hline_info.close()
    
    # This loop identifies telluric output containing warning about pixels outside calibration limits (different 
    # formatting)  Viraja:  This comment is from HlinesXD.py in XDGNIRS. I am not sure what the "different formatting"
    # means.
    if 'limits' in telluric_hLineInfo[-1].split()[-1]:
        normalization = telluric_hLineInfo[-2].split()[-1]
    else:
        normalization = telluric_hLineInfo[-1].split()[-1]
    
    # Divide normalization introduced by iraf.telluric
    iraf.imarith(operand1=outputtelspectrum, operand2=normalization, op='/', result=outputtelspectrum, title='',
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


# ----------------------------------------------------------------------------------------------------------------------
def fitTelluricContinuum(telluric_hLineCorrectedSpectrum, telluric_fitContinuum, continuumInter, orders,
                         continuumRegions, fitTelContinuumOrders, tempInter, overwrite):
    """
    Fit the continua of the H line corrected telluric 1D source spectra using IRAF CONTINUUM task to normalize them.
    """
    logger = log.getLogger('gnirsTelluric.fitTelluricContinuum')

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

#---------------------------------------------------------------------------------------------------------------------#

def divideTelluricContinuum(telluric_hLineCorrectedSpectrum, telluric_fitContinuum, telluric_dividedContinuum, orders,\
    overwrite):
    """
    Divide the H line corrected telluric 1D source spectra by its continuum fit to normalize it.
    """
    logger = log.getLogger('gnirsTelluric.divideTelluricContinuum')

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
        
#---------------------------------------------------------------------------------------------------------------------#

def telluricCorrection(sci_src_extracted_spectrum, telluric_dividedContinuum, science_dividedTelluricLines, \
    telluricInter, orders, sci_airmass, telluricRegions, sci_tel_infofile, overwrite):
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
    logger = log.getLogger('gnirsTelluric.telluricCorrection')

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

#---------------------------------------------------------------------------------------------------------------------#

def reintroduceTelluricContinuum(science_dividedTelluricLines, telluric_fitContinuum, science_correctedTelluric, \
    orders, overwrite):
    """
    Re-introduce telluric continuum shape into telluric lines removed science 1D source spectra.
    """
    logger = log.getLogger('gnirsTelluric.reintroduceTelluricContinuum')

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

#---------------------------------------------------------------------------------------------------------------------#
'''
# TODO(nat): linefitAuto and linefitManual could be useful at some point.
def lineFitAuto(spectrum, grating):
    """
    Automatically fit Lorentz profiles to lines defined in existing cur* files. Go to x position in cursor file and use 
    space bar to find spectrum at each of those points.
    """
    logger = log.getLogger('gnirsTelluric.lineFitAuto')

    specpos = iraf.bplot(images=spectrum+'[SCI,1]', cursor='cur'+grating, Stdout=1, StdoutG='/dev/null')
    specpose = str(specpos).split("'x,y,z(x):")
    nextcur = 'nextcur'+grating+'.txt'
    # Write line x,y info to file containing Lorentz fitting commands for bplot
    write_line_positions(nextcur, specpos)
    iraf.delete('final_tel_no_hLines_no_norm.fits',ver="no",go_ahead='yes',Stderr='/dev/null')
    # Fit and subtract Lorentz profiles. Might as well write output to file.
    iraf.bplot(images=spectrum+'[sci,1]',cursor='nextcur'+grating+'.txt', new_image='final_tel_no_hLines_no_norm', overwrite="yes",StdoutG='/dev/null',Stdout='Lorentz'+grating)

#---------------------------------------------------------------------------------------------------------------------#

def lineFitManual(spectrum, grating):
    """ 
    Enter splot so the user can fit and subtract lorentz (or, rather any) profiles.
    """
    logger = log.getLogger('gnirsTelluric.lineFitManual')

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
#---------------------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
