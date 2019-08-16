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

import os, sys, log, ConfigParser, gnirsHeaders
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
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs,iraf.onedspec,iraf.imutil)

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
    combinedsrc = config.get('Prefixes_and_Filenames', 'combinedsrc')
    combinedsky = config.get('Prefixes_and_Filenames', 'combinedsky')
    extractionStepwise = config.getboolean('extractSpectra1D','extractionStepwise')
    extractionStepSize = config.getboolean('extractSpectra1D','extractionStepSize')
    calculateSpectrumSNR = config.getboolean('gnirsPipeline','calculateSpectrumSNR')
    dividedTelluricPrefix = config.get('Prefixes_and_Filenames', 'dividedTelluricPrefix')
    dividedSciencePrefix_extractREG = config.get('Prefixes_and_Filenames', 'dividedSciencePrefix_extractREG')
    dividedSciencePrefix_extractFS = config.get('Prefixes_and_Filenames', 'dividedSciencePrefix_extractFS')
    dividedSciencePrefix_extractSW = config.get('Prefixes_and_Filenames', 'dividedSciencePrefix_extractSW')
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
        iraf.chdir(scipath)

        # Get symbolic paths to the std and tel directories in the sci directory and the runtime data directory
        stdpath = scipath + '/Standard/Intermediate'  # relative path/link expected to be at the top level of every sci directory
        telpath = scipath + '/Telluric/Intermediate'
        runtimedatapath = scipath + '/../runtimeData'

        # Record the right number of order expected according to the GNIRS XD configuration.
        if 'LB_SXD' in scipath:
            orders = [3, 4, 5]
        elif 'LB_LXD' in scipath:
            orders = [3, 4, 5, 6, 7, 8]
        elif 'SB_SXD' in scipath:
            orders = [3, 4, 5, 6, 7, 8]
        else:
            logger.error("###################################################################################")
            logger.error("###################################################################################")
            logger.error("#                                                                                 #")
            logger.error("#   ERROR in get telluric info: unknown GNIRS XD configuration. Exiting script.   #")
            logger.error("#                                                                                 #")
            logger.error("###################################################################################")
            logger.error("###################################################################################\n")
            raise SystemExit

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
            sci_telluricCorrected = sorted(glob.glob(stdpath + '/' + dividedStandardPrefix + nofits(combinedsrc) + \
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
            telluric_dividedContinuum = sorted(glob.glob(telpath + '/' + dividedTelluricPrefix + nofits(combinedsrc) + \
                '_order*_MEF.fits'))
            if len(sci_telluricCorrected) > 0:
                logger.info("Required continuum divided telluric source spectra available.")
                tel_header_info = fits.open(telluric_dividedContinuum).header
            else:
                logger.warning("Required continuum divided telluric spectra not available.")
                logger.warning("Please run gnirsTelluric.py to create the continuum divided spectra or provide them")
                logger.warning("manually in %s.", scipath)
                logger.warning("Exiting script.\n")
                raise SystemExit

        # TODO:  if extractionStepwise:

        logger.info("Required telluric corrected science and standard source spectra check complete.")
        
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
            
            stdName = tel_header_info[os.path.basename(tel_dividedContinuum[0])]['OBJECT']
            
            stdMagnitudeK = config.get(stdName,'stdMagnitudeK')
            stdMagnitudeH = config.get(stdName,'stdMagnitudeH')
            stdMagnitudeJ = config.get(stdName,'stdMagnitudeJ')
            stdMagniudes = [stdMagnitudeK, stdMagnitudeH, stdMagnitudeJ, stdMagnitudeJ, stdMagnitudeJ, stdMagnitudeJ]
            
            # Zero magnitude fluxes for bands corresponding to different orders [K, H, J, J, J, J] in erg/s/cm^2/A
            fluxZeroPoints = [4.283E-10, 1.13E-9, 3.129E-10, 3.129E-10, 3.129E-10, 3.129E-10]

            # EXPTIME keyword is the "Exposure time (s) for sum of all coadds"
            sciExptime = sci_header_info[os.path.basename(sci_telluricCorrected[0])]['EXPTIME']
            stdExptime = tel_header_info[os.path.basename(tel_dividedContinuum[0])]['EXPTIME']
            
            for i in range(len(orders)):
                
                extension = i+1
                stdTemperature = config.getfloat(stdName,'stdTemperature')
           
                logger.info("Converting magnitude to flux density for the telluric.")
                makeFLambda(stdMagnitudes[i], sciExptime, stdExptime, fluxZeroPoints[i], overwrite)
                logger.info("Completed converting magnitude to flux density for the telluric.")

                logger.info("Making a blackbody.")
                makeBlackBody(overwrite)
                logger.info("Completed making a blackbody.")

                logger.info("Calculating the blackbody scale to the telluric magnitude.")
                makeBlackBodyScale(overwrite)
                logger.info("Completed calculating the blackbody scale to the telluric magnitude.")

                logger.info("Scaling the blackbody to the telluric magnitude.")
                scaleBlackBody(overwrite)
                logger.info("Scaling the blackbody to the telluric magnitude.")

                logger.info("Applying the blackbody to telluric magnitude scaling fator to science.")
                multiplyByBlackBody(overwrite)
                logger.info("Completed applying the blackbody to telluric magnitude scaling fator to science.")
        
        else:
            logger.error("#######################################################################################")
            logger.error("#######################################################################################")
            logger.error("#                                                                                     #")
            logger.error("#   ERROR in get flux calibration: unknown flux calibration method. Exiting script.   #")
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
    iraf.chdir(path)

    return

##################################################################################################################
#                                                     ROUTINES                                                   #
##################################################################################################################

def nofits(filename):
    return filename.replace('.fits', '')

#---------------------------------------------------------------------------------------------------------------------#

def makeFLambda(magnitude, sci_exptime, std_exptime, flux_zero_point, overwrite):
    """
    Tasks:
    - Multiply magnitude expression by appropriate constant for the bandpass of the spectral order.
    - Multiply by the ratio of exposure times of the telluric to the science.
    - If no magnitude, set fLambda to 1. In this case, no absolute flux calibration will be performed.
    
    Returns:
    - fLambda: floating point constant.
    """
    logger = log.getLogger('gnirsTelluric.makeFLambda')
    
    if magnitude:  # if the magnitude value in the configuration file is not empty
        # Convert magnitude to erg/cm2/s/A, for a rough flux scaling
        # Account for standard/science exposure times
        logger.info("Performing absolute flux calibration accounting for the ratio of standard to science exposure")
        logger.info("times.")
        flambda = 10**(- float(magnitude)/2.5) * flux_zero_point * (std_exptime / sci_exptime)
        perform_absolute_fluxcalib = True
    except:  # if the magnitude value in the configuration file is empty, no absolute flux calibration is performed
        logger.info("Performing relative flux calibration accounting for the ratio of standard to science exposure")
        logger.info("times.")
        flambda = 1 * (std_exptime / sci_exptime)
        perform_absolute_fluxcalib = False

#---------------------------------------------------------------------------------------------------------------------#

def makeBlackBody(rawFrame, grating, log, overwrite):
    """
    - From the science header read the keywords for the wavelength calibration solution to convert from pixel to 
      wavelength scale.
    - Make scale factor: mean of black body overwrite fLambda.
    - Multiply blackbody spectrum by scale factor.
    Creates:
        - Unscaled blackbody, bbody.fits
        - A scaled 1D blackbody spectrum, scaledBlackBody.fits[0]
    """
    # Find the start and end wavelengths of the blackbody from our cube header.
    target_header = astropy.io.fits.open('../products_uncorrected/ctfbrsn'+rawFrame+'.fits')
    wstart = target_header[1].header['CRVAL3']
    wdelt = target_header[1].header['CD3_3']
    wend = wstart + (2040 * wdelt)
    crpix3 = target_header[1].header['CRPIX3']
    # Find the standard star temperature from 0_std_starRAWNAME.txt
    try:
        with open("0_std_star"+rawFrame+".txt", "r") as f:
            lines = f.read()
        # ['k', 'K', '7.615', '9700', 'h', 'H', '7.636', '9700', 'j', 'J', '7.686', '9700', 'j', 'J', '7.686', '9700']
        lines = lines.split()
        # Mag is entry after the grating, but may also be N/A. Check for that.
        for i in range(len(lines)):
            if grating in lines[i]:
                standardStarSpecTemperature = lines[i+2]
                logger.info("Read a standard star teff of " + str(standardStarSpecTemperature))
    except IOError:
        logger.info("No std_starRAWNAME.txt file found; setting to spec temperature to 9700K for a rough flux scaling")
        standardStarSpecTemperature = 9700
    if crpix3 != 1.:
        logger.info("WARNING in Reduce: CRPIX of wavelength axis not equal to one. Exiting flux calibration.")
        raise SystemExit
    # Make a blackbody for each of the 2040 NIFS spectral pixels.
    if os.path.exists("3_BBody"+rawFrame+".fits"):
        if overwrite:
            os.remove("3_BBody"+rawFrame+".fits")
            iraf.chdir(os.getcwd())
            iraf.mk1dspec(input="3_BBody"+rawFrame,output="",title='',ncols=2040,naps=1,header='',wstart=wstart,wend=wend,temperature=standardStarSpecTemperature)
            logger.info("\nMade a blackbody in 3_BBody{}.fits".format(rawFrame))
        else:
            logger.info("\nOutput exists and -overwrite not set - skipping production of unscaled black body")
    else:
        iraf.chdir(os.getcwd())
        iraf.mk1dspec(input="3_BBody"+rawFrame,output="",title='',ncols=2040,naps=1,header='',wstart=wstart,wend=wend,temperature=standardStarSpecTemperature)
        logger.info("\nMade a blackbody in 3_BBody{}.fits".format(rawFrame))

def makeBlackBodyScale(rawFrame, log, overwrite):
    """
    For now, scale the black body by the ratio of the black body mean flux to fLambda.
    """
    # Get the mean of the unscaled blackbody.
    # FOR SOME REASON, iraf.imstat is having problems opening the image. So I am using Numpy for now.
    #mean = iraf.imstat(images="3_BBody"+rawFrame+".fits", fields="mean", lower='INDEF', upper='INDEF', nclip=0, lsigma=3.0, usigma=3.0, binwidth=0.1, format='yes', cache='no', mode='al',Stdout=1)
    unscaled = astropy.io.fits.open("3_BBody"+rawFrame+".fits")
    data = unscaled[0].data
    mean = float(np.mean(data))
    #mean = float(mean[1].replace("'",""))

    # Get fLambda.
    with open("2_fLambda"+rawFrame+".txt", "r") as f:
        fLambda = f.read()
    fLambda = float(fLambda.split()[1])

    # Create the scale factor.
    bbodyScaleFactor = fLambda / mean

    # Write it to a file.
    if os.path.exists("4_bbodyScaleFactor"+rawFrame+".txt"):
        if overwrite:
            os.remove("4_bbodyScaleFactor"+rawFrame+".txt")
            logger.info("\nFound a scale factor of {}; writing to 4_bbodyScaleFactor{}.txt".format(bbodyScaleFactor, rawFrame))
            with open("4_bbodyScaleFactor"+rawFrame+".txt", "w") as f:
                f.write("bbodyScaleFactor: {} \n".format(bbodyScaleFactor))
        else:
            logger.info("\nOutput exists and -overwrite not set - skipping writing of 4_bbodyScaleFactor{}.txt".format(rawFrame))
    else:
        logger.info("\nFound a scale factor of {}; writing to 4_bbodyScaleFactor{}.txt".format(bbodyScaleFactor, rawFrame))
        with open("4_bbodyScaleFactor"+rawFrame+".txt", "w") as f:
            f.write("bbodyScaleFactor: {} \n".format(bbodyScaleFactor))

def scaleBlackBody(rawFrame, log, overwrite):
    """
    Scale a black body.
    """
    # Get scale factor from file.
    with open("4_bbodyScaleFactor"+rawFrame+".txt", "r") as f:
        line = f.read()
    bbodyScaleFactor = float(line.split()[1])
    if os.path.exists("5_scaledBBody"+rawFrame+".fits"):
        if overwrite:
            os.remove("5_scaledBBody"+rawFrame+".fits")
            # A bug involving iraf.gemini() causes imarith to fail here. Use astropy unless you fixed it.
            #iraf.imarith(operand1="3_BBody"+rawFrame, op="*", operand2=bbodyScaleFactor, result="5_scaledBBody"+rawFrame,title='',divzero=0.0,hparams='',pixtype='',calctype='',verbose='no',noact='no',mode='al')
            operand1 = astropy.io.fits.open("3_BBody"+rawFrame+".fits")[0].data
            operand2 = bbodyScaleFactor
            multiplied = operand1 * operand2
            hdu = astropy.io.fits.PrimaryHDU(multiplied)
            hdu.writeto("5_scaledBBody"+rawFrame+".fits")

            logger.info("\nCreated a scaled blackbody, 5_scaledBBody{}.fits".format(rawFrame))
        else:
            logger.info("\nOutput exists and -overwrite not set - skipping production of scaled black body")
    else:
        # A bug involving iraf.gemini() causes imarith to fail here. Use astropy unless you fixed it.
        #iraf.imarith(operand1="3_BBody"+rawFrame, op="*", operand2=bbodyScaleFactor, result="5_scaledBBody"+rawFrame,title='',divzero=0.0,hparams='',pixtype='',calctype='',verbose='no',noact='no',mode='al')
        operand1 = astropy.io.fits.open("3_BBody"+rawFrame+".fits")[0].data
        operand2 = bbodyScaleFactor
        multiplied = operand1 * operand2
        hdu = astropy.io.fits.PrimaryHDU(multiplied)
        hdu.writeto("5_scaledBBody"+rawFrame+".fits")
    # We now have a scaled blackbody, scaledBlackBody.fits

def multiplyByBlackBody(rawFrame, log, overwrite):
    """
    - Multiply each slice of continuum multiplied telluric corrected cube
      by the scaled black body.

    Creates:
        - Flux calibrated cube, "factfbrsn"+scienceObjectName+".fits"
    """
    # Open the telluric corrected, continuum multiplied, un-fluxcalibrated data cube.
    cube = astropy.io.fits.open('1_continuum'+rawFrame+'.fits')
    # Open the scaled blackbody. We will multiply the cube by this.
    scaledBlackBody = astropy.io.fits.open("5_scaledBBody"+rawFrame+".fits")

    if os.path.exists("factfbrsn"+rawFrame+'.fits'):
        if overwrite:
            os.remove('factfbrsn'+rawFrame+'.fits')
            # Divide each spectrum in the cubedata array by the telluric correction spectrum.
            for i in range(cube[1].header['NAXIS2']):         # NAXIS2 is the y axis of the final cube.
                for j in range(cube[1].header['NAXIS1']):     # NAXIS1 is the x axis of the final cube.
                    cube[1].data[:,i,j] *= (scaledBlackBody[0].data)
            # Write the corrected cube to a new file.
            cube.writeto('factfbrsn'+rawFrame+'.fits', output_verify='ignore')
        else:
            logger.info("\nOutput exists and -overwrite not set - skipping division of telluric corrected cube by scaled black body")
    else:
        # Divide each spectrum in the cubedata array by the telluric correction spectrum.
        for i in range(cube[1].header['NAXIS2']):         # NAXIS2 is the y axis of the final cube.
            for j in range(cube[1].header['NAXIS1']):     # NAXIS1 is the x axis of the final cube.
                cube[1].data[:,i,j] *= (scaledBlackBody[0].data)
        # Write the corrected cube to a new file.
        cube.writeto('factfbrsn'+rawFrame+'.fits', output_verify='ignore')
'''