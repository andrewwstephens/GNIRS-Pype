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
    """
    logger = log.getLogger('check_flux_calibrators')

    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)

    manualMode = config.getboolean('defaults', 'manualMode')
    runtimedata = config.get('defaults', 'runtimeData')
    combinedsrc = config.get('runtimeFilenames', 'combinedsrc')
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
        logger.debug('firstfile: %s', firstfile)
        target = olog[firstfile]['OBJECT']
        logger.info("Checking for flux calibration parameters for %s", target)

        param = {'RA': None, 'DEC': None, 'SpectralType': None, 'Temperature': None,
                 'J': None, 'H': None, 'K': None}

        if config.has_section(target):
            logger.info('Found section for the %s in the config file.', target)
            for option in param.keys():
                if config.has_option(target, option) and config.get(target, option):
                    logger.warning("%s is already set in the config file.", option)
                elif config.has_option(target, option) and not config.get(target, option):
                    logger.warning("%s exists but is not set in the config file.", option)
                    logger.debug("Setting %s to None", option)
                    config.set(target, option, None)
                else:
                    logger.info("%s is not in the config file.", option)
                    logger.debug("Setting %s to None", option)
                    config.set(target, option, None)
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
            config.set(target, 'RA', param['RA'])
            config.set(target, 'DEC', param['DEC'])
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

            logger.debug('SIMBAD results: %s', simbad_table)

            if param['SpectralType'] is None:
                try:
                    param['SpectralType'] = simbad_table['SP_TYPE'][0]
                    config.set(target, 'SpectralType', param['SpectralType'])
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
                        config.set(target, 'Temperature', param['Temperature'])
                        logger.info('Temperature: %d K', param['Temperature'])
                        break
                if param['Temperature'] is None:
                    logger.error('Could not find a temperature for spectral type %s', param['SpectralType'])
                    logger.error("Please manually update the 'Temperature' in the config file.")
                    raise SystemExit

            simbad_label = {'J': 'FLUX_J', 'H': 'FLUX_H', 'K': 'FLUX_K'}
            for mag in ['J', 'H', 'K']:
                if param[mag] is None:
                    try:
                        param[mag] = simbad_table[simbad_label[mag]][0]
                        logger.info('%s = %s', mag, param[mag])
                        config.set(target, mag, param[mag])
                    except:
                        logger.error('Could not find %s in the SIMBAD response.', mag)
                        logger.error('Please manually update the config file.')
                        raise SystemExit

        logger.debug('Writing config file...')
        with open(configfile, 'w') as f:
            config.write(f)

        logger.info(" --------------------------------- ")
        logger.info("| Flux calibrators check complete |")
        logger.info(" --------------------------------- ")

    return


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs-pype.cfg')
