#!/usr/bin/env python

import os, glob, log, ConfigParser, gnirsHeaders
from astroquery.simbad import Simbad
from astropy import units as u
import astropy.coordinates as coord


# ----------------------------------------------------------------------------------------------------------------------
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
    logger.info('#################################################')

    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)

    manualMode = config.getboolean('defaults', 'manualMode')
    runtimedata = config.get('defaults', 'runtimeData')
    combinedsrc = config.get('runtimeFilenames', 'combinedsrc')
    extractRegularPrefix = config.get('runtimeFilenames', 'extractRegularPrefix')
    hLinePrefix = config.get('runtimeFilenames', 'hLinePrefix')
    dividedTelContinuumPrefix = config.get('runtimeFilenames', 'dividedTelContinuumPrefix')
    referenceStarTableFilename = config.get('fluxCalibration', 'referenceStarTableFilename')

    for scipath in config.options("ScienceDirectories"):

        if not config.getboolean("ScienceDirectories", scipath):  # only process directories that are marked True
            logger.info('Skipping getting telluric info in %s', scipath)
            continue

        #########################################################################
        #                                                                       #
        #                  BEGIN - OBSERVATION SPECIFIC SETUP                   #
        #                                                                       #
        #########################################################################

        stdpath = scipath + '/Standard/Intermediate'
        telpath = scipath + '/Telluric/Intermediate'
        scipath += '/Intermediate'

        # Record the right number of orders expected according to the GNIRS XD configuration.
        if 'LB_SXD' in scipath:
            orders = [3, 4, 5]
        elif 'LB_LXD' in scipath:
            orders = [3, 4, 5, 6, 7, 8]
        elif 'SB_SXD' in scipath:
            orders = [3, 4, 5, 6, 7, 8]
        else:
            logger.error("###################################################################################")
            logger.error("#                                                                                 #")
            logger.error("#   ERROR in get telluric info: unknown GNIRS XD configuration. Exiting script.   #")
            logger.error("#                                                                                 #")
            logger.error("###################################################################################")
            raise SystemExit

        if os.path.exists(stdpath):
            logger.info("Standard directory: %s", stdpath)
            logger.info("Skipping getting standard info as standard directory found in %s", scipath)
            continue
        elif os.path.exists(telpath):
            logger.info("Standard directory does not exist.")
            logger.info("Moving on to use the telluric as the standard to derive the approximate flux calibration.")
            stdpath = telpath
            logger.info("Telluric (used as standard) directory: %s\n", stdpath)
        else:
            logger.error("############################################################################################")
            logger.error("#                                                                                          #")
            logger.error("#   ERROR in get telluric info: standard and telluric data not available. Exiting script.  #")
            logger.error("#                                                                                          #")
            logger.error("############################################################################################")
            raise SystemExit

        # Check if required runtime data files and telluric corrected science source spectra available in their
        # respective paths
        logger.info("Checking if required runtime data files and continuum divided telluric source spectra available")
        logger.info("in %s and %s, respectively.", runtimedata, stdpath)

        if os.path.exists(runtimedata + '/' + referenceStarTableFilename):
            logger.info("Required reference star table of spectral types and temperatures available.")
            referenceStarTable = open(runtimedata + '/' + referenceStarTableFilename, "r").readlines()
        else:
            logger.warning("Required reference star table of spectral types and temperatures not available.")
            logger.warning("Please provide the star table in %s.", runtimedata)
            logger.warning("Exiting script.\n")
            raise SystemExit

        tel_dividedContinuum = sorted(glob.glob(stdpath + '/' + dividedTelContinuumPrefix + hLinePrefix + \
            extractRegularPrefix + nofits(combinedsrc) + '_order*_MEF.fits'))
        if len(tel_dividedContinuum) > 0:
            logger.info("Required continuum divided telluric source spectra available.")
            tel_header_info = gnirsHeaders.info(tel_dividedContinuum[0])
        else:
            logger.warning("Required continuum divided telluric spectra not available.")
            logger.warning("Please run gnirsTelluric.py to create the continuum divided spectra or provide them")
            logger.warning("manually in %s.", scipath)
            logger.warning("Exiting script.\n")
            raise SystemExit

        logger.info("Required runtime data files and continuum divided telluric spectra check complete.\n")
        
        #########################################################################
        #                                                                       #
        #                 COMPLETE - OBSERVATION SPECIFIC SETUP                 #
        #               BEGIN GET TELLURIC INFO FOR AN OBSERVATION              #
        #                                                                       #
        #########################################################################
        
        if manualMode:
            a = raw_input("About to enter getting telluric (used as standard star) info.")

        stdName = tel_header_info[os.path.basename(tel_dividedContinuum[0])]['OBJECT']

        logger.info("Checking if there exists a section (possibly) with parameters for the telluric %s in", stdName)
        logger.info("the configuration file.")

        parameters = ['stdRA', 'stdDEC', 'stdSpectralType', 'stdTemperature', 'stdMagnitude_order3',
            'stdMagnitude_order4', 'stdMagnitude_order5', 'stdMagnitude_order6', 'stdMagnitude_order7',
            'stdMagnitude_order8']
        if config.has_section(stdName):
            logger.info("Section for the telluric %s available in the configuration file.", stdName)
            for option in parameters:
                if config.has_option(stdName,option) and config.get(stdName,option):
                    logger.warning("Parameter %s for %s available and set in the configuration file.", option, stdName)
                    logger.warning("Not updating parameter %s from the query to SIMBAD.", option)
                elif config.has_option(stdName,option) and not config.get(stdName,option):
                    logger.warning("Parameter %s for %s available but not set in the configuration file.", option, stdName)
                    logger.info("Setting an empty parameter %s for %s in the configuration file.", option, stdName)
                    config.set(stdName,option,None)
                else:
                    logger.info("Parameter %s for %s not available in the configuration file.", option, stdName)
                    logger.info("Setting an empty parameter %s for %s in the configuration file.", option, stdName)
                    config.set(stdName,option,None)
        else:
            logger.info("Section for the telluric %s not available in the configuration file.", stdName)
            logger.info("Creating the telluric section with empty options in the configuration file.\n")
            config.add_section(stdName)
            [config.set(stdName,option,None) for option in parameters]

        for option in config.options(stdName):
            if option not in parameters:
                logger.info("Parameter %s for %s not required in the configuration file.", option, stdName)
                logger.info("Removing parameter %s for %s in the configuration file.\n", option, stdName)
                config.remove_option(stdName,option)

        with open('../../'+configfile, 'w') as f:
            logger.info('Updating the configuration file if any empty parameers were set.\n')
            config.write(f)
       
        stdRA = config.get(stdName,'stdRA')
        stdDEC = config.get(stdName,'stdDEC')
        stdSpectralType = config.get(stdName,'stdSpectralType')
        stdTemperature = config.get(stdName,'stdTemperature')
        stdMagnitude_order3 = config.get(stdName,'stdMagnitude_order3')
        stdMagnitude_order4 = config.get(stdName,'stdMagnitude_order4')
        stdMagnitude_order5 = config.get(stdName,'stdMagnitude_order5')
        stdMagnitude_order6 = config.get(stdName,'stdMagnitude_order6')
        stdMagnitude_order7 = config.get(stdName,'stdMagnitude_order7')
        stdMagnitude_order8 = config.get(stdName,'stdMagnitude_order8')

        # If user did not specify stdRA or stdDEC, get them from the science extension [1] of the continuum divided 
        # telluric spectrum header.
        if stdRA:
            pass
        else:
            stdRA = tel_header_info[os.path.basename(tel_dividedContinuum[0])]['RA']
        if stdDEC:
            pass
        else:
            stdDEC = tel_header_info[os.path.basename(tel_dividedContinuum[0])]['DEC']
        
        # If user did not specify stdMagnitude_order3, stdMagnitude_order4, stdMagnitude_order5, stdMagnitude_order6, 
        # stdMagnitude_order7, stdMagnitude_order8, or (stdTemperature or stdSpectralType), query SIMBAD for them.
        if stdTemperature:
            findSpectralType = False
            findTemperature = False
        else:
            findSpectralType = True
            findTemperature = True

        if stdMagnitude_order3:
            findMagnitude_order3 = False
        else:
            findMagnitude_order3 = True

        if stdMagnitude_order4:
            findMagnitude_order4 = False
        else:
            findMagnitude_order4 = True

        if stdMagnitude_order5:
            findMagnitude_order5 = False
        else:
            findMagnitude_order5 = True 

        if stdMagnitude_order6:
            findMagnitude_order6 = False
        else:
            findMagnitude_order6 = True

        if stdMagnitude_order7:
            findMagnitude_order7 = False
        else:
            findMagnitude_order7 = True

        if stdMagnitude_order8:
            findMagnitude_order8 = False
        else:
            findMagnitude_order8 = True                 

        if findSpectralType or findTemperature or findMagnitude_order3 or findMagnitude_order4 or findMagnitude_order5\
            or findMagnitude_order6 or findMagnitude_order7 or findMagnitude_order8:
            logger.info("Preparing to query SIMBAD for the standard star parameters.")
            # Construct URL based on telluric coordinates and execute SIMBAD query to find the spectral type
            Simbad.add_votable_fields('flux(K)', 'flux(J)', 'flux(H)', 'sp')
            try:
                logger.info("Querying SIMBAD...")
                simbadStarTable = Simbad.query_region(coord.SkyCoord(ra=stdRA, dec=stdDEC, unit=(u.deg, u.deg),
                    frame='fk5'), radius=0.1 * u.deg)
                # Viraja:  How are the RA and DEC formatted in XDpiped.csh???
            except:
                logger.error("Connection to SIMBAD was not established - either it is a network connection issue or")
                logger.error("SIMBAD is temporarily unavailable.")
                logger.info("Please manually provide the required parameters for the telluric in the configuration")
                logger.info("file.")
                logger.warning("Exiting script.\n")
                raise SystemExit

            logger.info("SIMBAD query successful!\n")

            if findSpectralType:
                try:
                    # Get spectral type -- only the first 3 characters (strip off end of types like AIVn as they are not
                    # in the referenceStarTable)
                    stdSpectralType = simbadStarTable['SP_TYPE'][0][0:3]
                except:
                    logger.error("Cannot locate the spectral type of the telluric in the table generated by the")
                    logger.error("SIMBAD query. Please update the parameter 'stdSpectralType' in the")
                    logger.error("configuration file.")
                    raise SystemExit
            
            if findTemperature:
                # Find temperature for the spectral type in referenceStarTable in runtime data path
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
                    logger.error("Please manually update the parameter 'stdTemperature' in the configuration file.")
                    logger.error("Exiting script.\n")
                    raise SystemExit
            
            if findMagnitude_order3:
                try:
                    stdMagnitude_order3 = str(simbadStarTable['FLUX_K'][0])
                except:
                    logger.error("Cannot find a K magnitude for order 3 of the telluric in the table generated by the")
                    logger.error("SIMBAD query.")
                    logger.error("Please manually update the parameter 'stdMagnitudeK' in the configuration file.")
                    logger.error("Exiting script.\n")
                    raise SystemExit

            if findMagnitude_order4:
                try:
                    stdMagnitude_order4 = str(simbadStarTable['FLUX_H'][0])
                except:
                    logger.error("Cannot find a H magnitude for order 4 of the telluric in the table generated by the")
                    logger.error("SIMBAD query.")
                    logger.error("Please manually update the parameter 'stdMagnitudeH' in the configuration file.")
                    logger.error("Exiting script.\n")
                    raise SystemExit

            if findMagnitude_order5:
                try:
                    stdMagnitude_order5 = str(simbadStarTable['FLUX_J'][0])
                except:
                    logger.error("Cannot find a J magnitude for order 5 of the telluric in the table generated by the")
                    logger.error("SIMBAD query.")
                    logger.error("Please manually update the parameter 'stdMagnitudeJ' in the configuration file.")
                    logger.error("Exiting script.\n")
                    raise SystemExit
            
            if findMagnitude_order6:
                try:
                    stdMagnitude_order6 = str(simbadStarTable['FLUX_J'][0])
                except:
                    logger.error("Cannot find a J magnitude for order 6 of the telluric in the table generated by the")
                    logger.error("SIMBAD query.")
                    logger.error("Please manually update the parameter 'stdMagnitudeJ' in the configuration file.")
                    logger.error("Exiting script.\n")
                    raise SystemExit
            
            if findMagnitude_order7:
                try:
                    stdMagnitude_order7 = str(simbadStarTable['FLUX_J'][0])
                except:
                    logger.error("Cannot find a J magnitude for order 7 of the telluric in the table generated by the")
                    logger.error("SIMBAD query.")
                    logger.error("Please manually update the parameter 'stdMagnitudeJ' in the configuration file.")
                    logger.error("Exiting script.\n")
                    raise SystemExit
            
            if findMagnitude_order8:
                try:
                    stdMagnitude_order8 = str(simbadStarTable['FLUX_J'][0])
                except:
                    logger.error("Cannot find a J magnitude for order 8 of the telluric in the table generated by the")
                    logger.error("SIMBAD query.")
                    logger.error("Please manually update the parameter 'stdMagnitudeJ' in the configuration file.")
                    logger.error("Exiting script.\n")
                    raise SystemExit

        config.set(stdName,'stdRA',stdRA)
        config.set(stdName,'stdDEC',stdDEC)
        config.set(stdName,'stdSpectralType',stdSpectralType)
        config.set(stdName,'stdTemperature',stdTemperature)
        config.set(stdName,'stdMagnitude_order3',stdMagnitude_order3)
        config.set(stdName,'stdMagnitude_order4',stdMagnitude_order4)
        config.set(stdName,'stdMagnitude_order5',stdMagnitude_order5)
        config.set(stdName,'stdMagnitude_order6',stdMagnitude_order6)
        config.set(stdName,'stdMagnitude_order7',stdMagnitude_order7)
        config.set(stdName,'stdMagnitude_order8',stdMagnitude_order8)
        
        with open('../../'+configfile, 'w') as f:
            logger.info('Updating the configuration file with telluric parameters obtained from SIMBAD.')
            config.write(f)
        
        logger.info("Updated the following section to the configuration file:")
        logger.info("[%s]", stdName)
        logger.info("stdRA = %s", stdRA)
        logger.info("stdDEC = %s", stdDEC)
        logger.info("stdSpectralType = %s", stdSpectralType)
        logger.info("stdTemperature = %s", stdTemperature)
        logger.info("stdMagnitude_order3 = %s", stdMagnitude_order3)
        logger.info("stdMagnitude_order4 = %s", stdMagnitude_order4)
        logger.info("stdMagnitude_order5 = %s", stdMagnitude_order5)
        logger.info("stdMagnitude_order6 = %s", stdMagnitude_order6)
        logger.info("stdMagnitude_order7 = %s", stdMagnitude_order7)
        logger.info("stdMagnitude_order8 = %s\n", stdMagnitude_order8)

        logger.info("##############################################################################")
        logger.info("#                                                                            #")
        logger.info("#  COMPLETE - Getting telluric info completed for                            #")
        logger.info("#  %s", scipath)
        logger.info("#                                                                            #")
        logger.info("##############################################################################")

    return


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
