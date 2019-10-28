#!/usr/bin/env python

from astropy import units as u
from astroquery.simbad import Simbad
import ConfigParser
import math
import obslog
import log
import utils


# ----------------------------------------------------------------------------------------------------------------------
def start(configfile):

    logger = log.getLogger('get_redshift')
    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)
    manualMode = config.getboolean('defaults', 'manualMode')

    logger.info(' ------------------------- ')
    logger.info('| Getting Target redshift |')
    logger.info(' ------------------------- ')

    for scipath in config.options("ScienceDirectories"):

        if not config.getboolean("ScienceDirectories", scipath):  # only process directories that are marked True
            logger.info('Skipping %s', scipath)
            continue

        logger.info('%s', scipath)
        utils.pause(manualMode)

        olog = obslog.readcsv(scipath + '/Intermediate/obslog.csv')
        target = None
        for f in olog.keys():
            if olog[f]['OBSTYPE'] == 'OBJECT':
                target = olog[f]['OBJECT']
                break
        if target is None:
            logger.error('Could not get target from obslog.csv')
            raise SystemExit
        logger.info("Target: %s", target)

        if config.has_section(target):
            logger.info('Found config section for %s', target)
            if config.has_option(target, 'Redshift'):
                try:
                    redshift = float(config.get(target, 'Redshift'))
                    logger.info('Redshift: %s', redshift)
                    continue  # go to next target
                except:
                    logger.debug('Could not convert %s to float', config.get(target, 'Redshift'))
        else:
            logger.info('Creating a new config file section for %s', target)
            config.add_section(target)

        logger.info("Querying SIMBAD...")
        Simbad.add_votable_fields('rvz_radvel', 'rvz_type')
        try:
            results = Simbad.query_region(target, radius=0.01 * u.deg)
            logger.debug('Results: %s', results)
        except:
            logger.warning('SIMBAD query failed for the redshift of %s', target)
            results = None
            redshift = None

        if results:
            radvel = results['RVZ_RADVEL'][0]
            rvz_type = results['RVZ_TYPE'][0]
            logger.debug('RADVEL: %s', radvel)
            logger.debug('TYPE: %s', rvz_type)

            try:
                radvel = float(radvel)
            except:
                logger.error('Could not parse redshift from: %s', radvel)

            c = 299792.458  # km/s

            if rvz_type == 'z':    # redshift
                redshift = radvel
            elif rvz_type == 'v':  # velocity
                v = radvel
                redshift = math.sqrt((1.0 + v/c) / (1.0 - v/c)) - 1.0
            elif rvz_type == 'c':  # cz
                redshift = radvel / c
            else:
                logger.error('Unknown SIMBAD RVZ_TYPE: %s', rvz_type)
                redshift = None

        logger.info('Redshift: %s', redshift)

        logger.info('Updating configuration file...')
        config.set(target, 'Redshift', redshift)
        with open(configfile, 'w') as f:
            config.write(f)

    return


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs-pype.cfg')
