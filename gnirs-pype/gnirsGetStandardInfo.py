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

import os, glob, log, ConfigParser, pkg_resources, gnirsHeaders
from pyraf import iraf
from astroquery.simbad import Simbad
from astropy.table import Table
from astropy.io import fits
from astropy import units as u
import astropy.coordinates as coord


def start(configfile):
    """
    Get standard star parameters from the configuration file or SIMBAD.
    """
    logger = log.getLogger('gnirsGetStandardInfo.start')

    # Store current working directory for later use.
    path = os.getcwd()

    logger.info('#################################################')
    logger.info('#                                               #')
    logger.info('#         Start GNIRS Get Standard Info         #')
    logger.info('#                                               #')
    logger.info('#################################################\n')

    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)
    # Read general config
    manualMode = config.getboolean('defaults','manualMode')
    overwrite = config.getboolean('defaults','overwrite')
    # config required for getting standard info
    extractionStepwise = config.getboolean('extractSpectra1D','extractStepwise')
    combinedsrc = config.get('runtimeFilenames', 'combinedsrc')
    combinedsky = config.get('runtimeFilenames', 'combinedsky')
    dividedTelluricPrefix = config.get('runtimeFilenames', 'dividedTelluricPrefix')
    dividedSciencePrefix_extractREG = config.get('runtimeFilenames', 'dividedSciencePrefix_extractREG')
    dividedSciencePrefix_extractFS = config.get('runtimeFilenames', 'dividedSciencePrefix_extractFS')
    dividedSciencePrefix_extractSW = config.get('runtimeFilenames', 'dividedSciencePrefix_extractSW')
    referenceStarTableFilename = config.get('fluxCalibration','referenceStarTableFilename')

    # Get symbolic paths to the std and tel directories in the sci directory and the runtime data directory
    stdpath = '../Standard/Intermediate'  # relative path/link expected to be at the top level of every sci directory
    telpath = '../Telluric/Intermediate'
    runtimedatapath = '../../runtimeData'

    for scipath in config.options("ScienceDirectories"):
        
        if not config.getboolean("ScienceDirectories", scipath):  # only process directories that are marked True
            logger.info('Skipping getting standard info in %s', scipath)
            continue

        ###########################################################################
        ##                                                                       ##
        ##                  BEGIN - OBSERVATION SPECIFIC SETUP                   ##
        ##                                                                       ##
        ###########################################################################

        scipath += '/Intermediate'
        logger.info("Moving to science directory: %s", scipath)
        iraf.chdir(scipath)

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
            logger.error("#   ERROR in get standard info: unknown GNIRS XD configuration. Exiting script.   #")
            logger.error("#                                                                                 #")
            logger.error("###################################################################################")
            logger.error("###################################################################################\n")
            raise SystemExit

        if os.path.exists(stdpath):
            logger.info("Standard directory: %s", stdpath)
            logger.info("Skipping getting standard info as standard directory found in %s", scipath)
            continue
        elif os.path.exists(telpath):
            logger.info("Standard directory does not exist.")
            logger.info("Moving on to use the telluric to derive the approximate flux calibration.\n")
            stdpath = telpath
            logger.info("Standard directory: %s\n", stdpath)
        else:
            logger.error("#######################################################################################")
            logger.error("#######################################################################################")
            logger.error("#                                                                                     #")
            logger.error("#   ERROR in get standard info: Found no standard or telluric data. Exiting script.   #")
            logger.error("#                                                                                     #")
            logger.error("#######################################################################################")
            logger.error("#######################################################################################\n")
            raise SystemExit

        logger.info("Runtime data path: %s\n", runtimedatapath)

        # Check if required runtime data files and telluric corrected science source spectra available in their
        # respective paths
        logger.info("Checking if required runtime data files and continuum divided telluric source spectra available")
        logger.info("in %s and %s, respectively.\n", runtimedatapath, stdpath)

        if os.path.exists(runtimedatapath + '/' + referenceStarTableFilename):
            logger.info("Required reference star table of spectral types and temperatures available.")
            referenceStarTable = open(runtimedatapath + '/' + referenceStarTableFilename, "r").readlines()
        else:
            logger.warning("Required reference star table of spectral types and temperatures not available.")
            logger.warning("Please provide the star table in %s.", runtimedatapath)
            logger.warning("Exiting script.\n")
            raise SystemExit
        '''
        sci_telluricCorrected = sorted(glob.glob(scipath + '/' + dividedSciencePrefix_extractREG + \
            nofits(combinedsrc) + '_order*_MEF.fits'))
        if len(sci_telluricCorrected) > 0:
            logger.info("Required telluric corrected science source spectra available.")
            sci_header = fits.open(sci_telluricCorrected[0])[1].header  
        else:
            logger.warning("Required telluric corrected science source spectra not available.")
            logger.warning("Please run gnirsTelluric.py to create the telluric corected spectra or provide them")
            logger.warning("manually in %s.", scipath)
            logger.warning("Exiting script.\n")
            raise SystemExit
        '''
        tel_dividedContinuum = sorted(glob.glob(stdpath + '/' + dividedTelluricPrefix + nofits(combinedsrc) + \
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

        logger.info("Required runtime data files and continuum divided telluric spectra check complete.")
        
        ###########################################################################
        ##                                                                       ##
        ##                 COMPLETE - OBSERVATION SPECIFIC SETUP                 ##
        ##               BEGIN GET STANDARD INFO FOR AN OBSERVATION              ##
        ##                                                                       ##
        ###########################################################################
        
        stdName = tel_header_info[os.path.basename(tel_dividedContinuum[0])]['OBJECT']

        stdRA = config.get('telluricCorrection','stdRA')
        stdDEC = config.get('telluricCorrection','stdDEC')
        stdSpectralType = config.get('telluricCorrection','stdSpectralType')
        stdMagnitudeK = config.get('telluricCorrection','stdMagnitudeK')
        stdMagnitudeH = config.get('telluricCorrection','stdMagnitudeH')
        stdMagnitudeJ = config.get('telluricCorrection','stdMagnitudeJ')
        stdTemperature = config.get('telluricCorrection','stdTemperature')



            getTelluricInfo(tel_header, stdRA, stdDEC, stdSpectralType, stdMagnitudeK, \
                stdTemperature, reference_startable, telpath+'/telluric_info.txt', overwrite)

            logger.info("##################################################################")
            logger.info("#                                                                #")
            logger.info("#       STEP 1: Get telluric information - COMPLETED             #")
            logger.info("#                                                                #")
            logger.info("##################################################################\n")
            
            #############################################################################
            ##  STEP 2: Convert magnitude to flux density for the telluric.            ##
            ##  Output: Derived spectrum (FLambda) for the telluric.                   ##
            #############################################################################

            if valindex == 2:
                makeFLambda(rawFrame, grating, log, over)

                logger.info("##################################################################")
                logger.info("#                                                                #")
                logger.info("#       STEP 1: Get telluric information - COMPLETED             #")
                logger.info("#                                                                #")
                logger.info("##################################################################\n")

            #############################################################################
            ##  STEP 1: Get telluric information by querrying SIMBAD if not obtained   ## 
            ##          from the configuration file.                                   ##
            ##  Output: An ascii file containing telluric information.                 ##
            #############################################################################

            if valindex == 3:
                makeBlackBody(rawFrame, grating, log, over)

                logger.info("##################################################################")
                logger.info("#                                                                #")
                logger.info("#       STEP 1: Get telluric information - COMPLETED             #")
                logger.info("#                                                                #")
                logger.info("##################################################################\n")

            #############################################################################
            ##  STEP 1: Get telluric information by querrying SIMBAD if not obtained   ## 
            ##          from the configuration file.                                   ##
            ##  Output: An ascii file containing telluric information.                 ##
            #############################################################################

            if valindex == 4:
                makeBlackBodyScale(rawFrame, log, over)

                logger.info("##################################################################")
                logger.info("#                                                                #")
                logger.info("#       STEP 1: Get telluric information - COMPLETED             #")
                logger.info("#                                                                #")
                logger.info("##################################################################\n")

            #############################################################################
            ##  STEP 1: Get telluric information by querrying SIMBAD if not obtained   ## 
            ##          from the configuration file.                                   ##
            ##  Output: An ascii file containing telluric information.                 ##
            #############################################################################

            if valindex == 5:
                scaleBlackBody(rawFrame, log, over)

                logger.info("##################################################################")
                logger.info("#                                                                #")
                logger.info("#       STEP 1: Get telluric information - COMPLETED             #")
                logger.info("#                                                                #")
                logger.info("##################################################################\n")

            #############################################################################
            ##  STEP 1: Get telluric information by querrying SIMBAD if not obtained   ## 
            ##          from the configuration file.                                   ##
            ##  Output: An ascii file containing telluric information.                 ##
            #############################################################################

            if valindex == 6:
                multiplyByBlackBody(rawFrame, log, over)

                logger.info("##################################################################")
                logger.info("#                                                                #")
                logger.info("#       STEP 1: Get telluric information - COMPLETED             #")
                logger.info("#                                                                #")
                logger.info("##################################################################\n")
            
            else:
                logger.error("##############################################################################")
                logger.error("##############################################################################")
                logger.error("#                                                                            #")
                logger.error("#    ERROR in flux calibrate: %d is not valid. Exiting script.", valindex      )
                logger.error("#                                                                            #")
                logger.error("##############################################################################")
                logger.error("##############################################################################")
                raise SystemExit
            
            valindex += 1

        logger.info("##############################################################################")
        logger.info("#                                                                            #")
        logger.info("#  COMPLETE - Flux calibration completed for                                 #")
        logger.info("#  %s", scipath                                                                )
        logger.info("#                                                                            #")
        logger.info("##############################################################################\n")

    # Return to directory script was begun from.
    os.chdir(path)

    return
'''
##################################################################################################################
#                                                     ROUTINES                                                   #
##################################################################################################################

def nofits(filename):
    return filename.replace('.fits', '')

# ----------------------------------------------------------------------------------------------------------------------
'''
def getTelluricInfo(tel_header, stdRA, stdDEC, stdSpectralType, stdMagnitudeK, stdTemperature, \
    reference_startable, tel_infofile, overwrite):
    """
    Find telluric star spectral type, temperature, and K, H, and J magnitudes. Based on XDGNIRS code. Modified from
    nifsTelluric.py code.

    Executes a SIMBAD query and parses the resulting html to find spectal type, temperature, and K, H, and J magnitudes.

    Reads:
        Science extension [1] header of the order 0 continuum divided telluric source spectrum.
    
    Writes:
        Telluric information file in the telluric directory path.
    """
    logger = log.getLogger('gnirsTelluric.getTelluricInfo')

    # This code structure checks if output files already exist. If output files exist and overwrite is specified, all 
    # output is overwritten.
    if os.path.exists(tel_infofile):
        if overwrite:
            logger.warning("Removing old %s", tel_infofile)
            os.remove(tel_infofile)
        else:
            logger.info("Output exists and -overwrite not set - skipping getting telluric information.\n")
            return

        # If user did not specify a stdMagnitudeK, stdTemperature, stdRA, or stdDEC, get stdRA
        # and stdDEC from the extracted telluric spectrum header. Use SIMBAD to look up stdSpectralType, 
        # stdMagnitudeK, and stdTemperature.
        if (not telluricMagnitudeK or not telluricMagnitudeH or not telluricMagnitudeJ or not stdTemperature) \
            and (not stdRA or not stdDEC):
            if not stdRA:
                stdRA = tel_header['RA']
            if not stdDEC:
                stdDEC = tel_header['DEC']
            telExptime = str(tel_header['EXPTIME'])
        else:
            # Get the telluric exposure time anyways
            telExptime = str(tel_header['EXPTIME'])

        # Check to see if a spectral type or temperature has been given
        if stdTemperature:
            findSpectralType = False
            findTemperature = False
        else:
            findSpectralType = True
            findTemperature = True

        if telluricMagnitudeK:
            findMagnitudeK = False
        else:
            findMagnitudeK = True

        if telluricMagnitudeH:
            findMagnitudeH = False
        else:
            findMagnitudeH = True
        
        if telluricMagnitudeJ:
            findMagnitudeJ = False
        else:
            findMagnitudeJ = True

        if findSpectralType or findTemperature or findMagnitudeK or findMagnitudeH or findMagnitudeJ:
            # Construct URL based on telluric coordinates and execute SIMBAD query to find the spectral type
            Simbad.add_votable_fields('flux(K)', 'flux(J)', 'flux(H)', 'sp')
            simbad_startable = Simbad.query_region(coord.SkyCoord(ra=stdRA, dec=stdDEC, \
                unit=(u.deg, u.deg), frame='fk5'), radius=0.1 * u.deg)
            # Viraja:  How are the RA and DEC formatted in XDpiped.csh??

            if findSpectralType:
                # Get spectral type -- only the first 3 characters (strip off end of types like AIVn as they are not
                # in the 'reference_startable.txt')
                stdSpectralType = simbad_startable['SP_TYPE'][0][0:3]
            else:
                logger.error("Cannot locate the spectral type of the telluric in the table generated by the")
                logger.error("SIMBAD query. Please update the parameter 'stdSpectralType' in the")
                logger.error("configuration file.")
                raise SystemExit
            
            if findMagnitudeK:
                telluricMagnitudeK = str(simbad_startable['FLUX_K'][0])
            else:
                logger.error("Cannot find a the K magnitude for the telluric in the table generated by the")
                logger.error("SIMBAD query.")
                logger.error("Please manually update the parameter 'telluricMagnitudeK' in the configuration file.")
                logger.error("Exiting script.\n")
                raise SystemExit

            if findMagnitudeH:
                telluricMagnitudeH = str(simbad_startable['FLUX_H'][0])
            else:
                logger.error("Cannot find a the H magnitude for the telluric in the table generated by the")
                logger.error("SIMBAD query.")
                logger.error("Please manually update the parameter 'telluricMagnitudeH' in the configuration file.")
                logger.error("Exiting script.\n")
                raise SystemExit

            if findMagnitudeJ:
                telluricMagnitudeJ = str(simbad_startable['FLUX_J'][0])
            else:
                logger.error("Cannot find a the J magnitude for the telluric in the table generated by the")
                logger.error("SIMBAD query.")
                logger.error("Please manually update the parameter 'telluricMagnitudeJ' in the configuration file.")
                logger.error("Exiting script.\n")
                raise SystemExit

            if findTemperature:
                # Find temperature for the spectral type in 'reference_startable.txt'
                count = 0
                for line in reference_startable:
                    if '#' in line:
                        continue
                    else:
                        if stdSpectralType in line.split()[0]:
                            stdTemperature = line.split()[1]
                            count = 0
                            break
                        else:
                            count += 1
                if count > 0:  ## Viraja:  I am wondering why this condition is given and why is this an error??
                    logger.error("Cannot find a temperature for spectral type %s of the telluric.", stdSpectralType)
                    logger.error("Please update the parameter 'stdTemperature' in the configuration file.")
                    raise SystemExit

        tel_info = open(tel_infofile,'w')
        tel_info.write('k K ' + telluricMagnitudeK + ' ' + stdTemperature + '\n')
        tel_info.write('h H ' + telluricMagnitudeH + ' ' + stdTemperature + '\n')
        tel_info.write('j J ' + telluricMagnitudeJ + ' ' + stdTemperature + '\n')
        tel_info.write('j J ' + telluricMagnitudeJ + ' ' + stdTemperature + '\n')
        tel_info.write('j J ' + telluricMagnitudeJ + ' ' + stdTemperature + '\n')
        tel_info.write('j J ' + telluricMagnitudeJ + ' ' + stdTemperature + '\n')
        tel_info.close()
        
        tel_info = open(tel_infofile,'r').readlines()
        logger.info("Contents of %s:", tel_infofile)
        for line in tel_info:
            logger.info("%s", line.split('\n'))

#---------------------------------------------------------------------------------------------------------------------#

def makeFLambda(rawFrame, grating, log, over):
    """
    Tasks:
    - Multiply magnitude expression by appropriate constant for the bandpass of the spectral order.
    - Multiply by the ratio of exposure times of the telluric to the science.
    - If no magnitude, set fLambda to 1. In this case, no absolute flux calibration will be performed.
    
    Returns:
    - fLambda: floating point constant.
    """
    logger = log.getLogger('gnirsTelluric.makeFLambda')
    
    #account for standard star/science target exposure times
    #EXPTIME keyword is the "Exposure time (s) for sum of all coadds"



    # 2MASS doesn't have Z band magnitudes. Use J for rough absolute flux scaling.
    if grating == 'Z':
        grating = 'J'

    # Get standard star magnitude and exposure time
    try:
        with open("0_std_star"+rawFrame+".txt", "r") as f:
            lines = f.read()
        # Split into a list; it should then look something like this:
        # ['k', 'K', '7.615', '9700', 'h', 'H', '7.636', '9700', 'j', 'J', '7.686', '9700', 'j', 'J', '7.686', '9700']
        lines = lines.split()
        # Mag is entry after the grating, but may also be N/A. Check for that.
        for i in range(len(lines)):
            if grating in lines[i] and lines[i+1] != "N/A":
                standardStarMagnitude = lines[i+1]
                logger.info("Found standard star magnitude to be " + str(standardStarMagnitude))
        # Standard exposure time is the last thing in the list
        std_exp_time = lines[-1]
    except IOError:
        logger.info("No std_starRAWNAME.txt file found; no absolute flux cal will be attempted.")

    if grating == "K":
        constant = 4.283E-10
    elif grating == "H":
        constant = 1.13E-9
    elif grating == "J":
        constant = 3.129E-10
    elif grating == "Z":
        constant = 7.63E-9        

    # Look up the target exposure time.
    target_header = astropy.io.fits.open('../ctfbrsn'+rawFrame+'.fits')
    tgt_exp_time = target_header[0].header['EXPTIME']

    # Try actually making fLambda!
    try:
        standardStarMagnitude = float(standardStarMagnitude)
        fLambda = (10**((-1*standardStarMagnitude)/2.5)) * constant * (float(std_exp_time) / float(tgt_exp_time))
        logger.info("\nMade fLambda; doing rough absolute flux cal")
    except:
        # If no magnitude set to 1; no absolute flux cal. attempted
        fLambda = 1
        logger.info("\nCouldn't make fLambda; doing relative flux cal")
    if os.path.exists("2_fLambda"+rawFrame+".txt"):
        if over:
            os.remove("2_fLambda"+rawFrame+".txt")
            with open("2_fLambda"+rawFrame+".txt", "w") as f:
                f.write("fLambda: {}".format(fLambda))
            logger.info("\nWrote a value of {} to 2_fLambda{}.txt".format(fLambda, rawFrame))
        else:
            logger.info("\nOutput exists and -over not set - skipping write of fLambda to file")
    else:
        with open("2_fLambda"+rawFrame+".txt", "w") as f:
            f.write("fLambda: {}".format(fLambda))
        logger.info("\nWrote a value of {} to 2_fLambda{}.txt".format(fLambda, rawFrame))

def makeBlackBody(rawFrame, grating, log, over):
    """
    - From Z header information from the cube, make a black body.
    - Make scale factor: mean of black body over fLambda.
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
        if over:
            os.remove("3_BBody"+rawFrame+".fits")
            iraf.chdir(os.getcwd())
            iraf.mk1dspec(input="3_BBody"+rawFrame,output="",title='',ncols=2040,naps=1,header='',wstart=wstart,wend=wend,temperature=standardStarSpecTemperature)
            logger.info("\nMade a blackbody in 3_BBody{}.fits".format(rawFrame))
        else:
            logger.info("\nOutput exists and -over not set - skipping production of unscaled black body")
    else:
        iraf.chdir(os.getcwd())
        iraf.mk1dspec(input="3_BBody"+rawFrame,output="",title='',ncols=2040,naps=1,header='',wstart=wstart,wend=wend,temperature=standardStarSpecTemperature)
        logger.info("\nMade a blackbody in 3_BBody{}.fits".format(rawFrame))

def makeBlackBodyScale(rawFrame, log, over):
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
        if over:
            os.remove("4_bbodyScaleFactor"+rawFrame+".txt")
            logger.info("\nFound a scale factor of {}; writing to 4_bbodyScaleFactor{}.txt".format(bbodyScaleFactor, rawFrame))
            with open("4_bbodyScaleFactor"+rawFrame+".txt", "w") as f:
                f.write("bbodyScaleFactor: {} \n".format(bbodyScaleFactor))
        else:
            logger.info("\nOutput exists and -over not set - skipping writing of 4_bbodyScaleFactor{}.txt".format(rawFrame))
    else:
        logger.info("\nFound a scale factor of {}; writing to 4_bbodyScaleFactor{}.txt".format(bbodyScaleFactor, rawFrame))
        with open("4_bbodyScaleFactor"+rawFrame+".txt", "w") as f:
            f.write("bbodyScaleFactor: {} \n".format(bbodyScaleFactor))

def scaleBlackBody(rawFrame, log, over):
    """
    Scale a black body.
    """
    # Get scale factor from file.
    with open("4_bbodyScaleFactor"+rawFrame+".txt", "r") as f:
        line = f.read()
    bbodyScaleFactor = float(line.split()[1])
    if os.path.exists("5_scaledBBody"+rawFrame+".fits"):
        if over:
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
            logger.info("\nOutput exists and -over not set - skipping production of scaled black body")
    else:
        # A bug involving iraf.gemini() causes imarith to fail here. Use astropy unless you fixed it.
        #iraf.imarith(operand1="3_BBody"+rawFrame, op="*", operand2=bbodyScaleFactor, result="5_scaledBBody"+rawFrame,title='',divzero=0.0,hparams='',pixtype='',calctype='',verbose='no',noact='no',mode='al')
        operand1 = astropy.io.fits.open("3_BBody"+rawFrame+".fits")[0].data
        operand2 = bbodyScaleFactor
        multiplied = operand1 * operand2
        hdu = astropy.io.fits.PrimaryHDU(multiplied)
        hdu.writeto("5_scaledBBody"+rawFrame+".fits")
    # We now have a scaled blackbody, scaledBlackBody.fits

def multiplyByBlackBody(rawFrame, log, over):
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
        if over:
            os.remove('factfbrsn'+rawFrame+'.fits')
            # Divide each spectrum in the cubedata array by the telluric correction spectrum.
            for i in range(cube[1].header['NAXIS2']):         # NAXIS2 is the y axis of the final cube.
                for j in range(cube[1].header['NAXIS1']):     # NAXIS1 is the x axis of the final cube.
                    cube[1].data[:,i,j] *= (scaledBlackBody[0].data)
            # Write the corrected cube to a new file.
            cube.writeto('factfbrsn'+rawFrame+'.fits', output_verify='ignore')
        else:
            logger.info("\nOutput exists and -over not set - skipping division of telluric corrected cube by scaled black body")
    else:
        # Divide each spectrum in the cubedata array by the telluric correction spectrum.
        for i in range(cube[1].header['NAXIS2']):         # NAXIS2 is the y axis of the final cube.
            for j in range(cube[1].header['NAXIS1']):     # NAXIS1 is the x axis of the final cube.
                cube[1].data[:,i,j] *= (scaledBlackBody[0].data)
        # Write the corrected cube to a new file.
        cube.writeto('factfbrsn'+rawFrame+'.fits', output_verify='ignore')
'''

#---------------------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
