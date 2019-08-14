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

import os, sys, log, ConfigParser
from pyraf import iraf
from astroquery.simbad import Simbad
from astropy.table import Table
from astropy import units as u
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord


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

    # Set clobber to 'yes' for the script. This still does not make the gemini tasks overwrite files, so: YOU WILL 
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
    continuumInter = config.getboolean('interactive','continuumInter')
    extractionStepwise = config.getboolean('extractSpectra1D','extractionStepwise')
    extractionStepSize = config.getboolean('extractSpectra1D','extractionStepSize')
    calculateSpectrumSNR = config.getboolean('gnirsPipeline','calculateSpectrumSNR')
    # fluxCalibration specific config
    start = config.getint('telluricCorrection','Start')
    stop = config.getint('telluricCorrection','Stop')
    telluricRA = config.get('telluricCorrection','telluricRA')
    telluricDEC = config.get('telluricCorrection','telluricDEC')
    telluricSpectralType = config.get('telluricCorrection','telluricSpectralType')
    telluricMagnitude = config.get('telluricCorrection','telluricMagnitude')
    telluricTemperature = config.get('telluricCorrection','telluricTemperature')
    telluricBand = config.get('telluricCorrection','telluricBand')

    # Get symbolic path to the telluric directory within the science directory and the runtime data directory
    telpath = '../Telluric/Intermediate'  # relative path/link expected to be at the top level of every sci directory
    runtimedatapath = '../../runtimeData'

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

        logger.info("Telluric directory: %s\n", telpath)
        logger.info("Runtime data path: %s\n", runtimedatapath)
        
        os.chdir(scipath)

        # Check if required runtime data files and telluric corrected science source spectra available in their
        # respective paths
        logger.info("Checking if required runtime data files and telluric corrected science source spectra")
        logger.info("available in %s and %s, respectively.\n", runtimedatapath, scipath)
        
        reference_startable_filename = 'stars_spectraltypes_temperatures.txt'
        if os.path.exists(runtimedatapath+'/'+reference_startable_filename):
            logger.info("Required reference star table of spectral types and temperatures available.")
            reference_startable = open(runtimedatapath+'/'+reference_startable_filename, "r").readlines()
        else:
            logger.warning("Required reference star table of spectral types and temperatures not available.")
            logger.warning("Please provide the star table in %s. Exiting script.\n", runtimedatapath)
            raise SystemExit

        # Assigning different variable names to "src_comb.fits" and "sky_comb.fits" is not necessary with the
        # current file/directory arrangement for this pipeline (where this works as the telluric and science images
        # are not in the same directory), but it can be useful if the user decides to change the filenames of the
        # combined images.
        scisrccombimage = 'src_comb.fits'
        iraf.wmef(input=scipath + scisrccombimage[:scisrccombimage.rfind('.')] + '_MEF.fits', \
            output=scipath + 'h'scisrccombimage[:scisrccombimage.rfind('.')] + '_MEF.fits', extnames='', \
            phu=scipath + scisrccombimage, verbose='yes', mode='al')
        scitelcorrected = scipath+'/duv'+scisrccombimage
        sciHeader = fits.open(scitelcorrected).header  ## Use (file)[0].header for MEFs with PHU in extension 0.
        sci_airmass = sciHeader['AIRMASS']
        if os.path.exists(scitelcorrected):
            logger.info("Required telluric corrected science source spectrum available.")
        else:
            logger.warning("Required telluric corrected science source spectrum not available. Please run")
            logger.warning("gnirsTelluric.py to create the telluric corected spectrum or provide it manually in")
            logger.warning("%s. Exiting script.\n", scipath)
            raise SystemExit
        
        telsrccombimage = 'src_comb.fits'
        telsrcextractedspectrum = telpath+'/v'+telsrccombimage
        telHeader = fits.open(telsrcextractedspectrum)[0].header
        tel_airmass = telHeader['AIRMASS']
        if os.path.exists(telsrcextractedspectrum):
            logger.info("Required extracted telluric 1D source spectrum available.")
        else:
            logger.warning("Required extracted telluric 1D source spectrum not available. Please run")
            logger.warning("gnirsExtarctSpectra1D.py to create the extracted spectrum or provide it manually in")
            logger.warning("%s. Exiting script.\n", telpath)
            raise SystemExit

        # TODO:  if extractionStepwise:

        logger.info("Required telluric corrected spectra check complete.")

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
            logger.error("#     ERROR in flux calibrate: unknown GNIRS XD configuration. Exiting script.    #")
            logger.error("#                                                                                 #")
            logger.error("###################################################################################")
            logger.error("###################################################################################\n")
            raise SystemExit
        
        ###########################################################################
        ##                                                                       ##
        ##                 COMPLETE - OBSERVATION SPECIFIC SETUP                 ##
        ##             BEGIN TELLURIC CORRECTION FOR AN OBSERVATION              ##
        ##                                                                       ##
        ###########################################################################

        valindex = start
        while valindex > stop or valindex < 1 or stop > 6:
            logger.warning("#####################################################################")
            logger.warning("#####################################################################")
            logger.warning("#                                                                   #")
            logger.warning("#   WARNING in telluric: invalid start/stop values of telluric      #")
            logger.warning("#                        correction steps.                          #")
            logger.warning("#                                                                   #")
            logger.warning("#####################################################################")
            logger.warning("#####################################################################\n")

            valindex = int(raw_input("Please enter a valid start value (1 to 6, default 1): "))
            stop = int(raw_input("Please enter a valid stop value (1 to 6, default 6): "))

        while valindex <= stop:

            #############################################################################
            ##  STEP 1: Get telluric information by querrying SIMBAD if not obtained   ## 
            ##          from the configuration file.                                   ##
            ##  Output: An ascii file containing telluric information.                 ##
            #############################################################################

            if valindex == 1:
                getTelluricInfo(telHeader, telluricRA, telluricDEC, telluricSpectralType, telluricMagnitude, \
                    telluricTemperature, reference_startable, telpath+'/telluric_info.txt', overwrite)

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

##################################################################################################################
#                                                     ROUTINES                                                   #
##################################################################################################################

def getTelluricInfo(telHeader, telluricRA, telluricDEC, telluricSpectralType, telluricMagnitude, telluricTemperature, \
    reference_startable, telinfofile, overwrite):
    """
    Find telluric star spectral type, temperature, magnitude, and exposure time. Based on XDGNIRS code. Modified from
    nifsTelluric.py code.

    Executes a SIMBAD query and parses the resulting html to find spectal type, temperature, and magnitude.

    Reads:
        Extracted 1D telluric source spectrum.
    
    Writes:
        Telluric information file in the science directory path.
    """
    logger = log.getLogger('gnirsTelluric.getTelluricInfo')

    # This code structure checks if output files already exist. If output files exist and overwrite is specified, all 
    # output is overwritten.
    if os.path.exists(telinfofile):
        if overwrite:
            logger.warning("Removing old %s", telinfofile)
            os.remove(telinfofile)

            # If user did not specify a telluricMagnitude, telluricTemperature, telluricRA, or telluricDEC, get 
            # telluricRA and telluricDEC from the extracted telluric spectrum header. Use SIMBAD to look up 
            # telluricSpectralType, telluricMagnitude, and telluricTemperature.
            if (not telluricMagnitude or not telluricTemperature) and (not telluricRA or not telluricDEC):
                if not telluricRA:
                    telluricRA = telHeader['RA']
                if not telluricDEC:
                    telluricDEC = telHeader['DEC']
                telExptime = str(telHeader['EXPTIME'])
            else:
                # Get the telluric exposure time anyways
                telExptime = str(telHeader['EXPTIME'])

            # Check to see if a spectral type or temperature has been given
            if telluricTemperature:
                spectraltypeFind = False
                temperatureFind = False
            else:
                spectraltypeFind = True
                temperatureFind = True

            if telluricMagnitude:
                magnitudeFind = False
            else:
                magnitudeFind = True

            if spectraltypeFind or temperatureFind or magnitudeFind:
                # Construct URL based on telluric coordinates and execute SIMBAD query to find the spectral type
                Simbad.add_votable_fields('flux(K)', 'sp')  ## Viraja:  Why does it query only the K magnitude?
                simbad_startable = Simbad.query_region(coord.SkyCoord(ra=telluricRA, dec=telluricDEC, \
                    unit=(u.deg, u.deg), frame='fk5'), radius=0.1 * u.deg)
                # Viraja:  How are the RA and DEC formatted in XDpiped.csh??

                if spectraltypeFind:
                    # Get spectral type -- only the first 3 characters (strip off end of types like AIVn as they 
                    # are not in the 'reference_startable.txt')
                    telluricSpectralType = simbad_startable['SP_TYPE'][0][0:3]
                else:
                    logger.error("Cannot locate the spectral type of the telluric in the table generated by the")
                    logger.error("SIMBAD query. Please update the parameter 'telluricSpectralType' in the")
                    logger.error("configuration file.")
                    raise SystemExit
                
                if magnitudeFind:
                    telluricMagnitude = str(simbad_startable['FLUX_K'][0])
                else:
                    logger.error("Cannot find a the K magnitude for the telluric in the table generated by the")
                    logger.error("SIMBAD query. Please update the parameter 'telluricMagnitude' in the configuration")
                    logger.error("file. Exiting script.")
                    raise SystemExit

                if temperatureFind:
                    # Find temperature for the spectral type in 'reference_startable.txt'
                    count = 0
                    for line in reference_startable:
                        if '#' in line:
                            continue
                        else:
                            if telluricSpectralType in line.split()[0]:
                                telluricTemperature = line.split()[1]
                                count = 0
                                break
                            else:
                                count += 1
                    if count > 0:  ## Viraja:  I am wondering why this condition is given and why is this an error??
                        logger.error("Cannot find a temperature for spectral type %s of the telluric", telluricSpectralType)
                        logger.error("Please update the parameter 'telluricTemperature' in the configuration file.")
                        raise SystemExit

            telinfo = open(telinfofile,'w')
            telinfo.write('k K '+telluricMagnitude+' '+telluricTemperature+'\n')
            # Viraja: Why the same as the K for the rest of them?? --> comment in mag2mass.py in XDGNIRS -- "Not sure 
            # if the rest of these lines are necessary now we're not using "sensfunc" etc. for flux calibration
            telinfo.write('h H '+telluricMagnitude+' '+telluricTemperature+'\n')
            telinfo.write('j J '+telluricMagnitude+' '+telluricTemperature+'\n')
            telinfo.write('j J '+telluricMagnitude+' '+telluricTemperature+'\n')
            telinfo.write('i J '+telluricMagnitude+' '+telluricTemperature+'\n')
            telinfo.write('i J '+telluricMagnitude+' '+telluricTemperature+'\n')
            telinfo.close()
            
            telinfo = open(telinfofile,'r').readlines()
            logger.info("Contents of %s:", telinfofile)
            for line in telinfo:
                logger.info("%s", line.split('\n'))
        else:
            logger.info("Output exists and -overwrite not set - skipping getting telluric information.\n")
    else:
        if (not telluricMagnitude or not telluricTemperature) and (not telluricRA or not telluricDEC):
            if not telluricRA:
                telluricRA = telHeader['RA']
            if not telluricDEC:
                telluricDEC = telHeader['DEC']
            telExptime = str(telHeader['EXPTIME'])
        else:
            # Get the telluric exposure time regardless
            telHeader = fits.open(telsrcextractedspectrum)[0].header
            telExptime = str(telHeader['EXPTIME'])

        # Check to see if a spectral type or temperature has been given
        if telluricTemperature:
            spectraltypeFind = False
            temperatureFind = False
        else:
            spectraltypeFind = True
            temperatureFind = True

        if telluricMagnitude:
            magnitudeFind = False
        else:
            magnitudeFind = True

        if spectraltypeFind or temperatureFind or magnitudeFind:
            # Construct URL based on telluric coordinates and execute SIMBAD query to find the spectral type
            Simbad.add_votable_fields('flux(K)', 'sp')  ## Viraja:  Why does it query only the K magnitude?
            simbad_startable = Simbad.query_region(coord.SkyCoord(ra=telluricRA, dec=telluricDEC, unit=(u.deg, u.deg),\
                frame='fk5'), radius=0.1 * u.deg)  ## Viraja:  How are the RA and DEC formatted in XDpiped.csh??

            if spectraltypeFind:
                # Get spectral type -- only the first 3 characters (strip off end of types like AIVn as they are not in  
                # the 'reference_startable.txt')
                telluricSpectralType = simbad_startable['SP_TYPE'][0][0:3]
            else:
                logger.error("Cannot locate the spectral type of the telluric in the table generated by the SIMBAD")
                logger.error("query. Please update the parameter 'telluricSpectralType' in the configuration file.")
                raise SystemExit
            
            if magnitudeFind:
                telluricMagnitude = str(simbad_startable['FLUX_K'][0])
            else:
                logger.error("Cannot find a the K magnitude for the telluric in the table generated by the SIMBAD")
                logger.error("query. Please update the parameter 'telluricMagnitude' in the configuration file.")
                logger.error("Exiting script.")
                raise SystemExit

            if temperatureFind:
                # Find temperature for the spectral type in 'reference_startable.txt'
                count = 0
                for line in reference_startable:
                    if '#' in line:
                        continue
                    else:
                        if telluricSpectralType in line.split()[0]:
                            telluricTemperature = line.split()[1]
                            count = 0
                            break
                        else:
                            count += 1
                if count > 0:  ## Viraja:  I am wondering why this condition is given and why is this an error??
                    logger.error("Cannot find a temperature for spectral type %s of the telluric", telluricSpectralType)
                    logger.error("Please update the parameter 'telluricTemperature' in the configuration file.")
                    raise SystemExit

        telinfo = open(telinfofile,'w')
        telinfo.write('k K '+telluricMagnitude+' '+telluricTemperature+'\n')
        # Viraja: Why the same as the K for the rest of them?? --> comment in mag2mass.py in XDGNIRS -- "Not sure if 
        # the rest of these lines are necessary now we're not using "sensfunc" etc. for flux calibration
        telinfo.write('h H '+telluricMagnitude+' '+telluricTemperature+'\n')
        telinfo.write('j J '+telluricMagnitude+' '+telluricTemperature+'\n')
        telinfo.write('j J '+telluricMagnitude+' '+telluricTemperature+'\n')
        telinfo.write('i J '+telluricMagnitude+' '+telluricTemperature+'\n')
        telinfo.write('i J '+telluricMagnitude+' '+telluricTemperature+'\n')
        telinfo.close()
        
        telinfo = open(telinfofile,'r').readlines()
        logger.info("Contents of %s:", telinfofile)
        for line in telinfo:
            logger.info("%s", line.split('\n'))

#---------------------------------------------------------------------------------------------------------------------#

def makeFLambda(rawFrame, grating, log, over):
    """
    - Multiply magnitude expression by appropriate constant for grating.
    - Multiple by ratio of experiment times.
    - If no magnitude, set fLambda to 1. No absolute flux calibration performed.
    Returns:
        -fLambda: floating point constant.
    """
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
