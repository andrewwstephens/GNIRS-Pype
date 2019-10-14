#!/usr/bin/env python

from astropy import units as u
import astropy.coordinates as coords
from astroquery.simbad import Simbad
import ConfigParser
import obslog
import log
import os
import utils


# ----------------------------------------------------------------------------------------------------------------------
def start(configfile):
    """
    Make sure the config file includes all the required information for flux calibration:
       Spectral types, Temperatures, and Magnitudes for each flux calibrator
    The flux calibrator will either be in a 'Standard' directory or else we will use the Telluric standard.
    If the required parameters are not present try to query them from SIMBAD.

    Reads:
        Configuration file and or science extension [0] header of the order 0 continuum divided Telluric spectrum.
    
    Writes:
        Telluric parameters (if not already available) to the configuration file.
    """
    logger = log.getLogger('check_flux_calibrators')

    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)

    manualMode = config.getboolean('defaults', 'manualMode')
    runtimedata = config.get('defaults', 'runtimeData')
    combinedsrc = config.get('runtimeFilenames', 'combinedsrc')
    extractRegularPrefix = config.get('runtimeFilenames', 'extractRegularPrefix')
    hLinePrefix = config.get('runtimeFilenames', 'hLinePrefix')
    dividedTelContinuumPrefix = config.get('runtimeFilenames', 'dividedTelContinuumPrefix')
    stellar_temperature_data_file = config.get('fluxCalibration', 'StellarTemperatureData')

    utils.requires([runtimedata + '/' + stellar_temperature_data_file])
    stellar_temperature_data = open(runtimedata + '/' + stellar_temperature_data_file, 'r').readlines()

    for scipath in config.options("ScienceDirectories"):

        if not config.getboolean("ScienceDirectories", scipath):  # only process directories that are marked True
            logger.info('Skipping %s', scipath)
            continue

        logger.info(' -------------------------- ')
        logger.info('| Checking flux calibrator |')
        logger.info(' -------------------------- ')
        utils.pause(manualMode)

        logger.info('%s', scipath)
        stdpath = scipath + '/Standard/Intermediate'
        telpath = scipath + '/Telluric/Intermediate'
        scipath += '/Intermediate'

        if os.path.exists(stdpath):
            pass
        elif os.path.exists(telpath):
            stdpath = telpath
        else:
            logger.error('Could not find a Standard or Telluric directory')
            raise SystemExit

        # There are many ways to get the target name: FITS header, directory name, obslog...
        olog = obslog.readcsv(stdpath + '/obslog.csv')
        firstfile = olog.keys()[0]
        logger.debug('firsfile: %s', firstfile)
        target = olog[firstfile]['OBJECT']
        logger.info("Checking for flux calibration parameters for %s", target)

        param = {'RA': None, 'DEC': None, 'SpectralType': None, 'Temperature': None,
                 'Jmag': None, 'Hmag': None, 'Kmag': None}

        if config.has_section(target):
            logger.info('Found section for the %s in the config file.', target)
            for option in param.keys():
                if config.has_option(target, option) and config.get(target, option):
                    logger.warning("Parameter %s for %s available and set in the configuration file.", option, target)
                    logger.warning("Not updating parameter %s from the query to SIMBAD.", option)
                elif config.has_option(target, option) and not config.get(target, option):
                    logger.warning("Parameter %s for %s available but not set in the configuration file.", option, target)
                    logger.info("Setting an empty parameter %s for %s in the configuration file.", option, target)
                    config.set(target, option, None)
                else:
                    logger.info("Parameter %s for %s not available in the configuration file.", option, target)
                    logger.info("Setting an empty parameter %s for %s in the configuration file.", option, target)
                    config.set(target,option,None)
        else:
            logger.info('Creating a new config file section for %s', target)
            config.add_section(target)
            [config.set(target, option, None) for option in param]

        with open(configfile, 'w') as f:
            logger.info('Updating the configuration file...')
            config.write(f)

        for key in param.keys():
            param[key] = config.get(target, key)

        # If the RA and DEC are not in the config file then get them from the obslog:
        if param['RA'] is None or param['DEC'] is None:
            param['RA'] = olog[firstfile]['RA']
            param['DEC'] = olog[firstfile]['DEC']
            logger.debug('Coords: %s, %s', param['RA'], param['DEC'])

        # If no magnitudes are in the config file then query SIMBAD:
        if param['Temperature'] is None or param['Kmag'] is None:

            logger.info("Querying SIMBAD...")
            Simbad.add_votable_fields('flux(J)', 'flux(H)', 'flux(K)', 'sp')
            try:
                simbad_table = Simbad.query_region(
                    coords.SkyCoord(ra=param['RA'], dec=param['DEC'], unit=(u.deg, u.deg), frame='fk5'),
                    radius=0.1 * u.deg)
            except:
                logger.error('SIMBAD query failed')
                logger.info('Please manually provide the required parameters in the config file.')
                raise SystemExit

            logger.info('SIMBAD results: %s', simbad_table)

            if param['SpectralType'] is None:
                try:
                    param['SpectralType'] = simbad_table['SP_TYPE'][0]
                    logger.info('Spectral Type: %s', param['SpectralType'])
                except:
                    logger.error('No spectral type information found in the SIMBAD results.')
                    logger.error('Please update the parameter "SpectralType" in the config file.')
                    raise SystemExit

            if param['Temperature'] is None:  # Search the stellar temperatures database for the spectral type
                # Spectral types returned by simbad often include peculiarity codes, e.g. A0Vn or A0Vp.
                # For this reason search for table's class in the SIMBAD class (and not the other way around):
                for line in stellar_temperature_data:
                    field = line.split()
                    if field[0] in param['SpectralType']:
                        param['Temperature'] = float(field[1])
                        logger.info('Temperature: %d K', param['Temperature'])
                        break
                if param['Temperature'] is None:
                    logger.error('Could not find a temperature for spectral type %s', param['SpectralType'])
                    logger.error("Please manually update the 'Temperature' in the config file.")
                    raise SystemExit


            utils.pause(True, 'This is as far as you should go')

            


            if findMagnitude_order3:
                try:
                    stdMagnitude_order3 = str(simbad_table['FLUX_K'][0])
                except:
                    logger.error("Cannot find a K magnitude for order 3 of the telluric in the table generated by the")
                    logger.error("SIMBAD query.")
                    logger.error("Please manually update the parameter 'stdMagnitudeK' in the configuration file.")
                    logger.error("Exiting script.\n")
                    raise SystemExit

            if findMagnitude_order4:
                try:
                    stdMagnitude_order4 = str(simbad_table['FLUX_H'][0])
                except:
                    logger.error("Cannot find a H magnitude for order 4 of the telluric in the table generated by the")
                    logger.error("SIMBAD query.")
                    logger.error("Please manually update the parameter 'stdMagnitudeH' in the configuration file.")
                    logger.error("Exiting script.\n")
                    raise SystemExit

            if findMagnitude_order5:
                try:
                    stdMagnitude_order5 = str(simbad_table['FLUX_J'][0])
                except:
                    logger.error("Cannot find a J magnitude for order 5 of the telluric in the table generated by the")
                    logger.error("SIMBAD query.")
                    logger.error("Please manually update the parameter 'stdMagnitudeJ' in the configuration file.")
                    logger.error("Exiting script.\n")
                    raise SystemExit
            
            if findMagnitude_order6:
                try:
                    stdMagnitude_order6 = str(simbad_table['FLUX_J'][0])
                except:
                    logger.error("Cannot find a J magnitude for order 6 of the telluric in the table generated by the")
                    logger.error("SIMBAD query.")
                    logger.error("Please manually update the parameter 'stdMagnitudeJ' in the configuration file.")
                    logger.error("Exiting script.\n")
                    raise SystemExit
            
            if findMagnitude_order7:
                try:
                    stdMagnitude_order7 = str(simbad_table['FLUX_J'][0])
                except:
                    logger.error("Cannot find a J magnitude for order 7 of the telluric in the table generated by the")
                    logger.error("SIMBAD query.")
                    logger.error("Please manually update the parameter 'stdMagnitudeJ' in the configuration file.")
                    logger.error("Exiting script.\n")
                    raise SystemExit
            
            if findMagnitude_order8:
                try:
                    stdMagnitude_order8 = str(simbad_table['FLUX_J'][0])
                except:
                    logger.error("Cannot find a J magnitude for order 8 of the telluric in the table generated by the")
                    logger.error("SIMBAD query.")
                    logger.error("Please manually update the parameter 'stdMagnitudeJ' in the configuration file.")
                    logger.error("Exiting script.\n")
                    raise SystemExit

        config.set(target, 'stdRA', stdRA)
        config.set(target, 'stdDEC', stdDEC)
        config.set(target, 'stdSpectralType', stdSpectralType)
        config.set(target, 'stdTemperature', stdTemperature)
        config.set(target, 'stdMagnitude_order3', stdMagnitude_order3)
        config.set(target, 'stdMagnitude_order4', stdMagnitude_order4)
        config.set(target, 'stdMagnitude_order5', stdMagnitude_order5)
        config.set(target, 'stdMagnitude_order6', stdMagnitude_order6)
        config.set(target, 'stdMagnitude_order7', stdMagnitude_order7)
        config.set(target, 'stdMagnitude_order8', stdMagnitude_order8)
        
        with open('../../'+configfile, 'w') as f:
            logger.info('Updating the configuration file with telluric parameters obtained from SIMBAD.')
            config.write(f)
        
        logger.info("Updated the following section to the configuration file:")
        logger.info("[%s]", target)
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
        logger.info("#  COMPLETE - Getting telluric info completed for                            #")
        logger.info("#  %s", scipath)
        logger.info("##############################################################################")

    return


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs-pype.cfg')
