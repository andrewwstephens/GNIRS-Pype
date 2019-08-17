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

import os, sys, glob, log, ConfigParser, gnirsHeaders
from pyraf import iraf
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

#---------------------------------------------------------------------------------------------------------------------#

def start(configfile):
    """
    Do a flux calibration.
    """
    logger = log.getLogger('gnirsFluxCalibrate.start')

    # Store current working directory for later use.
    path = os.getcwd()

    logger.info('#################################################')
    logger.info('#                                               #')
    logger.info('#         Start GNIRS Flux Calibration          #')
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
    combinedsrc = config.get('runtimeFilenames', 'combinedsrc')
    combinedsky = config.get('runtimeFilenames', 'combinedsky')
    bb_unscaled = config.get('runtimeFilenames', 'bb_unscaled')
    bb_scaled = config.get('runtimeFilenames','bb_scaled')
    extractFullSlit = config.getboolean('extractSpectra1D','extractFullSlit')
    extractStepwise = config.getboolean('extractSpectra1D','extractStepwise')
    extractStepSize = config.getfloat('extractSpectra1D','extractStepSize')
    dividedTelluricPrefix = config.get('runtimeFilenames', 'dividedTelluricPrefix')
    dividedSciencePrefix_extractREG = config.get('runtimeFilenames', 'dividedSciencePrefix_extractREG')
    dividedSciencePrefix_extractFS = config.get('runtimeFilenames', 'dividedSciencePrefix_extractFS')
    dividedSciencePrefix_extractSW = config.get('runtimeFilenames', 'dividedSciencePrefix_extractSW')
    # fluxCalibration specific config
    fluxCalibrationMethod = config.get('fluxCalibration','fluxCalibrationMethod')


    for scipath in config.options("ScienceDirectories"):
        
        if not config.getboolean("ScienceDirectories", scipath):  # only process directories that are marked True
            logger.info('Skipping flux calibaration in %s', scipath)
            continue

        ###########################################################################
        ##                                                                       ##
        ##                  BEGIN - OBSERVATION SPECIFIC SETUP                   ##
        ##                                                                       ##
        ###########################################################################

        scipath += '/Intermediate'
        logger.info("Moving to science directory: %s", scipath)
        os.chdir(scipath)
        iraf.chdir(scipath)

        # Get symbolic paths to the std and tel directories in the sci directory and the runtime data directory
        # Relative path/link expected to be at the top level of every sci directory
        stdpath = '../Standard/Intermediate' 
        telpath = '../Telluric/Intermediate'
        runtimedatapath = '../../runtimeData'

        if os.path.exists(stdpath):
            logger.info("Standard directory: %s", stdpath)
            if fluxCalibrationMethod == 'fluxcalibrator':
                logger.info("Using 'fluxcalibrator' method for flux calibration in %s", scipath)
            elif fluxCalibrationMethod == 'telluricapproximate':
                logger.warning("Standard directory for flux calibrator available, but flux calibration method set to")
                logger.warning("'telluricapproximate'. This may affect your science.")
                fluxCalibrationMethod = raw_input("Confirm <fluxcalibrator> or <telluricapproximate> for flux calibration:")
            else:
                logger.error("#######################################################################################")
                logger.error("#######################################################################################")
                logger.error("#                                                                                     #")
                logger.error("#   ERROR in get flux calibration: unknown flux calibration method. Exiting script.   #")
                logger.error("#                                                                                     #")
                logger.error("#######################################################################################")
                logger.error("#######################################################################################\n")
                raise SystemExit
        elif os.path.exists(telpath):
            logger.info("Standard directory does not exist.")
            if fluxCalibrationMethod == 'fluxcalibrator':
                logger.error("Standard directory for flux calibrator not available in %s and flux", scipath)
                logger.error("calibration method set to 'fluxcalibrator'. Expects 'telluricapproximate'.")
                fluxCalibrationMethod = raw_input("Confirm <fluxcalibrator> or <telluricapproximate> for flux calibration:")
            elif fluxCalibrationMethod == 'telluricapproximate':
                logger.info("Moving on to use the telluric to derive the approximate flux calibration.\n")
                stdpath = telpath
                logger.info("Standard directory: %s\n", stdpath)
            else:
                logger.error("#######################################################################################")
                logger.error("#######################################################################################")
                logger.error("#                                                                                     #")
                logger.error("#   ERROR in get flux calibration: unknown flux calibration method. Exiting script.   #")
                logger.error("#                                                                                     #")
                logger.error("#######################################################################################")
                logger.error("#######################################################################################\n")
                raise SystemExit
        else:
            logger.warning("Parameter 'fluxCalibration' is set to 'yes', but no standard and telluric data available.")
            logger.warning("Skipping flux calibration in %s", scipath)
            continue

        # Check if required sci and std source spectra available in their respective paths
        logger.info("Checking if required science and standard source spectra available in %s",scipath)
        logger.info("and %s, respectively.\n", stdpath)

        sci_telluricCorrected = sorted(glob.glob(scipath + '/' + dividedSciencePrefix_extractREG + \
            nofits(combinedsrc) + '_order*_MEF.fits'))
        if len(sci_telluricCorrected) > 0:
            logger.info("Required telluric corrected science source spectra available.")
            sci_header_info = gnirsHeaders.info(sci_telluricCorrected[0])
        else:
            logger.warning("Required telluric corrected science source spectra not available.")
            logger.warning("Please run gnirsTelluric.py to create the telluric corected spectra or provide them")
            logger.warning("manually in %s.", scipath)
            logger.warning("Exiting script.\n")
            raise SystemExit

        if fluxCalibrationMethod == 'fluxcalibrator':
            # Get all telluric corrected standard source spectra files
            # The following codes lines are assuming that the standard spectra are reduced and telluric corrected in 
            # the same way as the science and are in a directory with a symbolic link to the top level processed data 
            # directory (similar to the telluric symbolic link in the science directory)
            '''
            std_telluricCorrected = sorted(glob.glob(stdpath + '/' + dividedStandardPrefix + nofits(combinedsrc) + \
                '_order*_MEF.fits'))
            if len(std_telluricCorrected) > 0:
                logger.info("Required telluric corrected standard source spectra available.")
                std_header_info = gnirsHeaders.info(std_telluricCorrected[0])  
            else:
                logger.warning("Required telluric corrected standard source spectra not available.")
                logger.warning("Please run gnirsTelluric.py to create the telluric corected spectra or provide them")
                logger.warning("manually in %s.", stdpath)
                logger.warning("Exiting script.\n")
                raise SystemExit
            '''
            pass
        
        if fluxCalibrationMethod == 'telluricapproximate':
            tel_dividedContinuum = sorted(glob.glob(telpath + '/' + dividedTelluricPrefix + nofits(combinedsrc) + \
                '_order*_MEF.fits'))
            if len(tel_dividedContinuum) > 0:
                logger.info("Required continuum divided telluric source spectra available.")
                tel_header_info = gnirsHeaders.info(tel_dividedContinuum[0])
            else:
                logger.warning("Required continuum divided telluric spectra not available.")
                logger.warning("Please run gnirsTelluric.py to create the continuum divided spectra or provide them")
                logger.warning("manually in %s.", scipath)
                logger.warning("Exiting script.\n")
                raise SystemExit

        # TODO:  if extractionStepwise:

        logger.info("Required telluric corrected science and standard source spectra check complete.")
        
        stdName = tel_header_info[os.path.basename(tel_dividedContinuum[0])]['OBJECT']

        # Record the right number of order expected according to the GNIRS XD configuration.
        if 'LB_SXD' in scipath:
            orders = [3, 4, 5]
            stdMagnitudes = [config.get(stdName,'stdMagnitude_order3'), config.get(stdName,'stdMagnitude_order4'), \
                config.get(stdName,'stdMagnitude_order5')]
            zeroMagnitudeFluxes = [config.getfloat(stdName,'zeroMagnitudeFlux_order3'), \
                config.getfloat(stdName,'zeroMagnitudeFlux_order4'), config.getfloat(stdName,'zeroMagnitudeFlux_order5')]
        elif 'LB_LXD' in scipath:
            orders = [3, 4, 5, 6, 7, 8]
            stdMagnitudes = [config.getfloat(stdName,'stdMagnitude_order3'), config.getfloat(stdName,'stdMagnitude_order4'), \
                config.getfloat(stdName,'stdMagnitude_order5'), config.getfloat(stdName,'stdMagnitude_order6'), \
                config.getfloat(stdName,'stdMagnitude_order7'), config.getfloat(stdName,'stdMagnitude_order8')]
            zeroMagnitudeFluxes = [config.getfloat(stdName,'zeroMagnitudeFlux_order3'), \
                config.getfloat(stdName,'zeroMagnitudeFlux_order4'), config.getfloat(stdName,'zeroMagnitudeFlux_order5'), \
                config.getfloat(stdName,'zeroMagnitudeFlux_order6'), config.getfloat(stdName,'zeroMagnitudeFlux_order7'), \
                config.getfloat(stdName,'zeroMagnitudeFlux_order8')]
        elif 'SB_SXD' in scipath:
            orders = [3, 4, 5, 6, 7, 8]
            stdMagnitudes = [config.getfloat(stdName,'stdMagnitude_order3'), config.getfloat(stdName,'stdMagnitude_order4'), \
                config.getfloat(stdName,'stdMagnitude_order5'), config.getfloat(stdName,'stdMagnitude_order6'), \
                config.getfloat(stdName,'stdMagnitude_order7'), config.getfloat(stdName,'stdMagnitude_order8')]
            zeroMagnitudeFluxes = [config.getfloat(stdName,'zeroMagnitudeFlux_order3'), \
                config.getfloat(stdName,'zeroMagnitudeFlux_order4'), config.getfloat(stdName,'zeroMagnitudeFlux_order5'), \
                config.getfloat(stdName,'zeroMagnitudeFlux_order6'), config.getfloat(stdName,'zeroMagnitudeFlux_order7'), \
                config.getfloat(stdName,'zeroMagnitudeFlux_order8')]
        else:
            logger.error("###################################################################################")
            logger.error("###################################################################################")
            logger.error("#                                                                                 #")
            logger.error("#   ERROR in get telluric info: unknown GNIRS XD configuration. Exiting script.   #")
            logger.error("#                                                                                 #")
            logger.error("###################################################################################")
            logger.error("###################################################################################\n")
            raise SystemExit

        ###########################################################################
        ##                                                                       ##
        ##                 COMPLETE - OBSERVATION SPECIFIC SETUP                 ##
        ##               BEGIN FLUX CALIBRATION FOR AN OBSERVATION               ##
        ##                                                                       ##
        ###########################################################################

        if fluxCalibrationMethod == 'fluxcalibrator':
            # This block is not currently set up.
            pass

        elif fluxCalibrationMethod == 'telluricapproximate':
            """
            - Convert magnitude to flux density for the telluric.  
            - Derived spectrum (FLambda) for the telluric for each order.

            - Make a blackbody using the telluric temperature.
            - A blackbody at the telluric temperature for each order.

            - Get the scale for the blackbody to the telluric derived spectrum.
            - Blackbody scale (float) at each order.

            - Scale the blackbody to the telluric derived spectrum.
            - Scaled blackbody at each order.
            """

            if manualMode:
                a = raw_input("About to enter ")

            stdTemperature = config.getfloat(stdName,'stdTemperature')

            # EXPTIME keyword is the "Exposure time (s) for sum of all coadds"
            sciExptime = sci_header_info[os.path.basename(sci_telluricCorrected[0])]['EXPTIME']
            logger.debug("sciExptime: %f", sciExptime)
            stdExptime = tel_header_info[os.path.basename(tel_dividedContinuum[0])]['EXPTIME']
            
            for i in range(len(orders)):
                
                logger.info("Flux-calibrating order %d of %s.", extension, os.path.basename(sci_telluricCorrected[i]))

                logger.info("Converting magnitude to flux density for the telluric.")
                # if the magnitude value in the configuration file is not empty
                if stdMagnitudes[i]:  
                    # Convert magnitude to erg/cm2/s/A, for a rough flux scaling
                    # Account for standard/science exposure times
                    logger.info("Performing absolute flux calibration accounting for the ratio of standard to science")
                    logger.info("exposure times.")
                    flambda = 10**(- float(stdMagnitudes[i])/2.5) * zeroMagnitudeFluxes[i] * (stdExptime / sciExptime)
                    absolute_fluxcalib = True
                # if the magnitude value in the configuration file is empty, no absolute flux calibration is performed
                else:  
                    logger.info("Performing relative flux calibration accounting for the ratio of standard to science")
                    logger.info("exposure times.")
                    flambda = 1 * (stdExptime / sciExptime)
                    absolute_fluxcalib = False
                logger.info("Completed converting magnitude to flux density for the telluric.\n")

                logger.info("Making a blackbody.")
                # First find the start and end wavelengths of the spectral order
                waveReferencePixel = iraf.hselect(images=tel_dividedContinuum[i]+'['+str(orders[i])+']', fields='CRPIX1', \
                    expr='yes',  missing='INDEF', mode='al', Stdout=1)
                logger.debug("waveReferencePixel: %d", waveReferencePixel)
                waveReferencePixel = float(waveReferencePixel[0].replace("'",""))
                logger.debug("waveReferencePixel: %d", waveReferencePixel)

                waveReferenceValue = iraf.hselect(images=tel_dividedContinuum[i]+'['+str(orders[i])+']', fields='CRVAL1', \
                    expr='yes',  missing='INDEF', mode='al', Stdout=1)
                waveReferenceValue = float(waveReferenceVlue[0].replace("'",""))
                
                waveDelt = iraf.hselect(images=tel_dividedContinuum[i]+'['+str(orders[i])+']', fields='CD1_1', \
                    expr='yes', missing='INDEF', mode='al', Stdout=1)
                waveDelt = float(waveDelt[0].replace("'",""))
                
                ndimensions = iraf.hselect(images=tel_dividedContinuum[i]+'['+str(orders[i])+']', fields='NAXIS1', \
                    expr='yes', missing='INDEF', mode='al', Stdout=1)
                ndimensions = float(ndimensions[0].replace("'",""))
                
                waveStart = waveReferenceValue - (waveReferencePixel-1) * waveDelt
                waveEnd = waveStart + (ndimensions * float(waveDelt[0].replace("'","")))
                # Then make a blackbody
                iraf.mk1dspec(input=bb_unscaled+str(orders[i]), output=bb_unscaled+str(orders[i]), ap=1, rv=0.0, \
                    z='no', title='', ncols=waveDelt, naps=1, header='', wstart = waveStart, wend=waveEnd, continuum=1000,\
                    slope=0.0, temperature=stdTemperature, fnu='no', lines='', nlines=0, profile='gaussian', \
                    peak=-0.5, gfwhm=20.0, lfwhm=20.0, seed=1, comments='yes', mode='ql')
                logger.info("Completed making a blackbody.\n")

                logger.info("Scaling the blackbody to the telluric magnitude.")
                if 3 <= extension <= 5:
                    # Roughly scale blackbody for orders 3-5 to the respective magnitudes of the telluric
                    logger.info("Calculating the blackbody scale for order %d.", orders[i])
                    
                    meanCounts = iraf.imstat(images=bb_unscaled+str(orders[i]), fields="mean", lower='INDEF', \
                        upper='INDEF', nclip=0, lsigma=3.0, usigma=3.0, binwidth=0.1, format='yes', cache='no', \
                        mode='al', Stdout=1)

                    blackbodyScaleFactor = flambda / float(meanCounts[1].replace("'",""))
                    logger.info("The blackbody scale factor for order %d is %f", orders[i], blackbodyScaleFactor)
                else:
                    # Scale blackbody for orders 6-8 to the previous order's scaled blackbody

                    # Ideally, waveEnd should be > waveStart.  If it is not, then there is something with the data that 
                    # is not necessarily a big problem, scientifically, but beeds some handling.  So, check for this
                    # and exit if there is no overlap between any two orders
                    
                    # First, find region of overlap with the previous order 
                    logger.info("Calculating the region of overlap of order %d with the previous order.", orders[i])
                    
                    # waveStart_current is the first wavelength of the current order (i.e, the short-wavelength end of 
                    # the region of overlap)
                    waveStart = iraf.hselect(images=bb_unscaled+str(orders[i]), fields='CRVAL1', expr='yes', missing='INDEF', \
                        mode='al', Stdout=1)
                    waveStart = float(waveStart[0].replace("'",""))
                    
                    # waveEnd_previous is the last wavelength of the previous order (i.e, the long-wavelength end of 
                    # the region of overlap)
                    waveTemp = iraf.hselect(bb_unscaled+str(orders[i-1]), field='CRVAL1', expr='yes',  missing='INDEF', \
                        mode='al', Stdout=1)
                    waveDelt = iraf.hselect(bb_unscaled+str(orders[i-1]), field='CDELT1', expr='yes',  missing='INDEF', \
                        mode='al', Stdout=1)
                    waveEnd = float(waveTemp[0].replace("'","")) + (1022 * float(waveDelt[0].replace("'","")))

                    if waveEnd < waveStart:
                        logging.error("Orders %d and %d do not overlap in wavelength.", orders[i-1], [orderi])
                        loggeing.error("This is unusual and suggests that the grating was not at the expected position.")
                        logging.error("This may not be a problem for the scientific use of the data, but the script")
                        logging.error("cannot handle this and is not able to flux calibrate the spectral orders.") 
                        logging.error("Please plot the calibrated arc spectrum (plotted with different orders) to see")
                        logging.error("if the data cover the wavelength range you need.")
                        logging.error("Exiting script.\n")
                        raise SystemExit
                    else:
                        # Find the mean in the overlapping wavelength region (using the scaled blackbody of the 
                        # previous order
                        iraf.scopy(input=bb_unscaled+str(orders[i]), output='temp'+bb_unscaled, w1=waveStart, w2=waveEnd, 
                            apertures='', bands='', beams='', apmodulus=0, format='multispec', renumber='no', offset=0, 
                            clobber='no', merge='no', rebin='yes', verbose='no',mode='ql')
                        iraf.scopy(input=bb_scaled+str(i), output='temp'+bb_scaled, w1=waveStart, w2=waveEnd,
                            apertures='', bands='', beams='', apmodulus=0, format='multispec', renumber='no', offset=0,
                            clobber='no', merge='no', rebin='yes', verbose='no', mode='ql')
                        meanCounts_temp_bb_unscaled = iraf.imstat(images='temp'+bb_unscaled+str(orders[i]), fields="mean",
                            lower='INDEF', upper='INDEF', nclip=0, lsigma=3.0, usigma=3.0, binwidth=0.1, format='yes', 
                            cache='no', mode='al',Stdout=1)
                        meanCounts_temp_bb_unscaled = float(meanCounts_temp_bb_unscaled[1].replace("'",""))
                        meanCounts_temp_bb_scaled = iraf.imstat(images='temp'+bb_scaled, fields='mean', lower='INDEF', 
                            upper='INDEF', nclip=0, lsigma=3.0, usigma=3.0, binwidth=0.1, format='yes', cache='no', 
                            mode='al', Stdout=1)
                        meanCounts_temp_bb_scaled = float(meanCounts_temp_bb_scaled[1].replace("'",""))
                        
                        # Scale current blackbody to that for the previous order using ratio of means
                        blackbodyScaleFactor = meanCounts_temp_bb_scaled / meanCounts_temp_bb_unscaled
                        logger.info("The blackbody scale factor for order %d is %f", orders[i], blackbodyScaleFactor)

                iraf.imarith(operand1=bb_unscaled + str(orders[i]), op="*", operand2=blackbodyScaleFactor,
                    result=bb_scaled + str(extensorders[i]ion), title='',divzero=0.0, hparams='', pixtype='', calctype='',
                    verbose='yes', noact='no', mode='al')
                logger.info("Completed scaling the blackbody to the telluric magnitude.")

                logger.info("Applying the blackbody to telluric magnitude scaling fator to telluric corrected science")
                logger.info("spectra.")
                fluxCalibrationInput = os.path.basename(sci_telluricCorrected[i])+'['+str(orders[i])+']'
                fluxCalibrationOutput_SEF = fluxCalibPrefix_extractREG + os.path.basename(sci_telluricCorrected[i])
                fluxCalibrationOutput_MEF = fluxCalibPrefix_extractREG + os.path.basename(sci_telluricCorrected[i])
                iraf.imarith(operand1=fluxCalibrationInput, op="*", operand2=bb_scaled + str(orders[i]),
                    result=fluxCalibrationOutput_SEF, title='', divzero=0.0, hparams='', pixtype='', calctype='', 
                    verbose='yes', noact='no', mode='al')
                # Record flux density units in headers
                fluxcalib_units_in_headers(fluxCalibrationOutput_SEF, absolute_fluxcalib)

                iraf.wmef(input=fluxCalibrationOutput_SEF, output=fluxCalibrationOutput_MEF, extnames='',
                    phu=os.path.basename(sci_telluricCorrected[i]), verbose='yes', mode='al')

                if extractFullSlit:
                    '''
                    fluxCalibrationInput = os.path.basename()+'[SCI,'+str(orders[i])+']'
                    fluxCalibrationOutput_SEF = fluxCalibPrefix_extractFS + os.path.basename()
                    fluxCalibrationOutput_MEF = fluxCalibPrefix_extractFS + os.path.basename()
                    iraf.imarith(operand1=fluxCalibrationInput, op="*", operand2=bb_scaled + str(orders[i]), 
                        result=fluxCalibrationOutput_SEF, title='', divzero=0.0, hparams='', pixtype='', 
                        calctype='', verbose='yes', noact='no', mode='al')
                    fluxcalib_units_in_headers("flamfull"+str(j), absolute_fluxcalib)
                    '''
                    pass
                if extractStepwise:
                    '''
                    fluxCalibrationInput = os.path.basename()+'[SCI,'+str(orders[i])+']'
                    fluxCalibrationOutput_SEF = fluxCalibPrefix_extractFS + os.path.basename()
                    fluxCalibrationOutput_MEF = fluxCalibPrefix_extractFS + os.path.basename()
                    for k in range(1, steps):
                        iraf.imarith(operand1=fluxCalibrationInput, op="*", operand2=bb_scaled + str(orders[i]), 
                            result=fluxCalibrationOutput_SEF, title='', divzero=0.0, hparams='', pixtype='',
                            calctype='', verbose='yes', noact='no', mode='al')
                        fluxcalib_units_in_headers(fluxCalibrationInput, absolute_fluxcalib)
                    '''
                    pass
                logger.info("Completed applying the blackbody to telluric magnitude scaling fator to telluric")
                logger.info("corrected science spectra.")

                logger.info("Completed flux calibrating order %d of %s.", orders[i], os.path.basename(sci_telluricCorrected[i]))

        else:
            logger.error("#######################################################################################")
            logger.error("#######################################################################################")
            logger.error("#                                                                                     #")
            logger.error("#     ERROR in flux calibration: unknown flux calibration method. Exiting script.     #")
            logger.error("#                                                                                     #")
            logger.error("#######################################################################################")
            logger.error("#######################################################################################\n")
            raise SystemExit

        logger.info("##############################################################################")
        logger.info("#                                                                            #")
        logger.info("#  COMPLETE - Flux calibration completed for                                 #")
        logger.info("#  %s", scipath                                                                )
        logger.info("#                                                                            #")
        logger.info("##############################################################################\n")

    # Return to directory script was begun from.
    os.chdir(path)

    return

##################################################################################################################
#                                                     ROUTINES                                                   #
##################################################################################################################

def nofits(filename):
    return filename.replace('.fits', '')

#---------------------------------------------------------------------------------------------------------------------#

def fluxcalib_units_in_headers(image, abs_fluxcal):
    """
    Add appropriate flux calibration units (whether absolute or relative flux calibration) to the image headers after
    flux calibration is done. 
    """
    logger = log.getLogger('gnirsFluxCalibration.fluxcalib_units_in_headers')

    #This is so we know whether we did absolute or relative flux cal
    if abs_fluxcal:
		iraf.hedit(images=image, fields='FUNITS', value='erg/cm^2/s/A', add='yes', addonly='no', delete='no',
            verify='no', show='no', update='yes')
    else:
        iraf.hedit(images=image, fields='FUNITS', value='Flambda, relative', add='yes', addonly='no', delete='no',
            verify='no', show='no', update='yes')

#---------------------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
