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
from astroquery.simbad import Simbad
from astropy.table import Table
from astropy.io import fits
from astropy import units as u
import astropy.coordinates as coord

#---------------------------------------------------------------------------------------------------------------------#

def start(configfile):
    """
    - Checks if standard star directory exists in the science path; else, uses the telluric to do an approximate flux 
    calibration.
    - If approximate flux calibration is to be done, checks if the configuration files has telluric (now used as the 
    stanndard star) parameters (spectral type, temperature, and K, H, and J magnitudes). 
    - If parameters available, skips querrying SIMBAD; else, querries SIMBAD to get the parameters and edits them to 
    the configuration file. Gets standard star parameters (spectral type, temperature, and K, H, and J magnitudes) from
    the configuration file or executes a SIMBAD query parses the resulting table to find spectal type, temperature, and
    K, H, and J magnitudes. 

    Based on XDGNIRS code. Modified from nifsTelluric.py code.
    
    Reads:
        Configuration file and or science extension [0] header of the order 0 continuum divided telluric source 
        spectrum.
    
    Writes:
        Telluric parameters (if not already available) to the configuration file.
    """
    logger = log.getLogger('gnirsGetTelluricInfo.start')

    # Store current working directory for later use.
    path = os.getcwd()

    logger.info('#################################################')
    logger.info('#                                               #')
    logger.info('#         Start GNIRS Get Telluric Info         #')
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
            logger.info('Skipping getting telluric info in %s', scipath)
            continue

        ###########################################################################
        ##                                                                       ##
        ##                  BEGIN - OBSERVATION SPECIFIC SETUP                   ##
        ##                                                                       ##
        ###########################################################################

        scipath += '/Intermediate'
        logger.info("Moving to science directory: %s", scipath)
        os.chdir(scipath)

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
            logger.error("#   ERROR in get telluric info: Found no standard or telluric data. Exiting script.   #")
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
        ##               BEGIN GET TELLURIC INFO FOR AN OBSERVATION              ##
        ##                                                                       ##
        ###########################################################################
        
        if manualMode:
            a = raw_input("About to enter getting telluric (used as standard star) info.")

        stdName = tel_header_info[os.path.basename(tel_dividedContinuum[0])]['OBJECT']

        logger.info("Checking if there exists a section (possibly) with parameters for the telluric %s. in", stdName)
        logger.info("the configuration file.")

        parameters = ['stdRA', 'stdDEC', 'stdSpectralType', 'stdTemperature', 'stdTemperature', 'stdMagnitudeK', \
            'stdMagnitudeH', 'stdMagnitudeJ']
        if config.has_section(stdName):
            logger.info("Section for the telluric %s available in the configuration file.", stdName)
            for option in parameters:
                if option in parameters and option in config.options():
                    if overwrite:
                        logger.warning("Removing current parameter %s for %s in the configuration file.", option, stdName)
                        logger.info("Adding empty parameter %s for %s in the configuration file.", option, stdName)
                        config.remove_option(stdName,option)
                        config.set(stdName,option,'')
                    else:
                        logger.warning("Parameter %s for %s available in the configuration file and -overwrite not")
                        logger.info("set - no action taken.")
                elif option in parameters and option not in config.options():
                    logger.info("Parameter %s for %s not available in the configuration file.", option, stdName)
                    logger.info("Adding Parameter %s for %s in the configuration file.", option, stdName)
                    [config.set(stdName,option,'') for option in parameters]
        else:
            logger.info("Section for the telluric %s not found in the configuration file.", stdName)
            logger.info("Adding the telluric section with empty options in the configuration file.")
            config.add_section(stdName)
            [config.set(stdName,option,'') for option in parameters]
      
        with open(configfile, 'w') as f:
            config.write(f)

        stdRA = config.get(stdName,'stdRA')
        stdDEC = config.get(stdName,'stdDEC')
        stdSpectralType = config.get(stdName,'stdSpectralType')
        stdTemperature = config.get(stdName,'stdTemperature')
        stdMagnitudeK = config.get(stdName,'stdMagnitudeK')
        stdMagnitudeH = config.get(stdName,'stdMagnitudeH')
        stdMagnitudeJ = config.get(stdName,'stdMagnitudeJ')

        # If user did not specify stdMagnitudeK, stdMagnitudeH, stdMagnitudeJ, or stdTemperature, and stdRA or stdDEC, 
        # get stdRA and stdDEC from the science extension [1] of the continuum divided telluric spectrum header. Use 
        # SIMBAD to look up stdSpectralType, stdMagnitudeK, stdMagnitudeH, stdMagnitudeJ, and stdTemperature.
        if (not stdMagnitudeK or not stdMagnitudeH or not stdMagnitudeJ or not stdTemperature) \
            and (not stdRA or not stdDEC):
            if not stdRA:
                stdRA = tel_header_info[os.path.basename(tel_dividedContinuum[0])]['RA']
            if not stdDEC:
                stdDEC = tel_header_info[os.path.basename(tel_dividedContinuum[0])]['DEC']

        # Check to see if a spectral type or temperature has been given
        if stdTemperature:
            findSpectralType = False
            findTemperature = False
        else:
            findSpectralType = True
            findTemperature = True

        if stdMagnitudeK:
            findMagnitudeK = False
        else:
            findMagnitudeK = True

        if stdMagnitudeH:
            findMagnitudeH = False
        else:
            findMagnitudeH = True
        
        if stdMagnitudeJ:
            findMagnitudeJ = False
        else:
            findMagnitudeJ = True

        if findSpectralType or findTemperature or findMagnitudeK or findMagnitudeH or findMagnitudeJ:
            # Construct URL based on telluric coordinates and execute SIMBAD query to find the spectral type
            Simbad.add_votable_fields('flux(K)', 'flux(J)', 'flux(H)', 'sp')
            try:
                simbadStarTable = Simbad.query_region(coord.SkyCoord(ra=stdRA, dec=stdDEC, \
                    unit=(u.deg, u.deg), frame='fk5'), radius=0.1 * u.deg)
                # Viraja:  How are the RA and DEC formatted in XDpiped.csh??
            except:
                logger.error("Connection to SIMBAD was not established - either it is a network connection issue or")
                logger.error("SIMBAD is temporarily unavailable.")
                logger.info("Please manually provide the required parameters for the telluric in the configuration")
                logger.info("file.")
                logger.warning("Exiting script.\n")
                raise SystemExit

            if findSpectralType:
                # Get spectral type -- only the first 3 characters (strip off end of types like AIVn as they are not
                # in the 'reference_startable.txt')
                stdSpectralType = simbadStarTable['SP_TYPE'][0][0:3]
            else:
                logger.error("Cannot locate the spectral type of the telluric in the table generated by the")
                logger.error("SIMBAD query. Please update the parameter 'stdSpectralType' in the")
                logger.error("configuration file.")
                raise SystemExit
            
            if findMagnitudeK:
                stdMagnitudeK = str(simbadStarTable['FLUX_K'][0])
            else:
                logger.error("Cannot find a the K magnitude for the telluric in the table generated by the")
                logger.error("SIMBAD query.")
                logger.error("Please manually update the parameter 'stdMagnitudeK' in the configuration file.")
                logger.error("Exiting script.\n")
                raise SystemExit

            if findMagnitudeH:
                stdMagnitudeH = str(simbadStarTable['FLUX_H'][0])
            else:
                logger.error("Cannot find a the H magnitude for the telluric in the table generated by the")
                logger.error("SIMBAD query.")
                logger.error("Please manually update the parameter 'stdMagnitudeH' in the configuration file.")
                logger.error("Exiting script.\n")
                raise SystemExit

            if findMagnitudeJ:
                stdMagnitudeJ = str(simbadStarTable['FLUX_J'][0])
            else:
                logger.error("Cannot find a the J magnitude for the telluric in the table generated by the")
                logger.error("SIMBAD query.")
                logger.error("Please manually update the parameter 'stdMagnitudeJ' in the configuration file.")
                logger.error("Exiting script.\n")
                raise SystemExit

            if findTemperature:
                # Find temperature for the spectral type in 'reference_startable.txt'
                count = 0
                for line in referenceStarTable:
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

        tel_info = open('telinfofile.txt','w')
        tel_info.write('k K ' + stdMagnitudeK + ' ' + stdTemperature + '\n')
        tel_info.write('h H ' + stdMagnitudeH + ' ' + stdTemperature + '\n')
        tel_info.write('j J ' + stdMagnitudeJ + ' ' + stdTemperature + '\n')
        tel_info.write('j J ' + stdMagnitudeJ + ' ' + stdTemperature + '\n')
        tel_info.write('j J ' + stdMagnitudeJ + ' ' + stdTemperature + '\n')
        tel_info.write('j J ' + stdMagnitudeJ + ' ' + stdTemperature + '\n')
        tel_info.close()
        
        tel_info = open('telinfofile.txt','r').readlines()
        logger.info("Contents of telinfofile.txt:")
        for line in tel_info:
            logger.info("%s", line.split('\n'))

        logger.info("##############################################################################")
        logger.info("#                                                                            #")
        logger.info("#  COMPLETE - Getting telluric info completed for                            #")
        logger.info("#  %s", scipath                                                                )
        logger.info("#                                                                            #")
        logger.info("##############################################################################\n")

    os.chdir(path)  ## Return to the original directory

    return

#---------------------------------------------------------------------------------------------------------------------#

def nofits(filename):
    return filename.replace('.fits', '')
 
#---------------------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
