#!/usr/bin/env python
"""
Sort GNIRS data into directories
"""
import os
import shutil
import ConfigParser
import header
import log


# ----------------------------------------------------------------------------------------------------------------------
def start(configfile):
    """
    Sort and copy the GNIRS raw data to subdirectories
    :param configfile: string text configuration file

    Target_Date_Config_ObsID
    - Calibrations [for Tellurics this is a link to the best Science Calibrations directory]
    - Intermediate
    - Final
    - Telluric [for Science targets this is a link to best Telluric directory]
    """
    logger = log.getLogger('sort_data')

    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make config file options case-sensitive
    config.read(configfile)
    rawpath = config.get('DownloadData', 'RawDataDir')

    info = header.info(rawpath)

    path = os.getcwd()
    logger.info("Path to raw data: %s", rawpath)
    logger.info('Sorted data will be copied to %s', path)
    logger.info("Creating new directories and copying files...")

    caldirs = []   # List of output calibration directories
    scidirs = []   # List of output science directories
    teldirs = []   # List of output Telluric standard directories

    for filename in sorted(info.keys()):

        logger.debug('%s: OBSCLASS:%s OBSTYPE:%s OBJECT:%s SLIT:%s GCALLAMP:%s',
                     filename, info[filename]['OBSCLASS'], info[filename]['OBSTYPE'], info[filename]['OBJECT'],
                     info[filename]['SLIT'], info[filename]['GCALLAMP'])

        if info[filename]['OBSTYPE'] == 'OBJECT' and 'acq' not in info[filename]['OBSCLASS']:      # SCIENCE & TELLURICS

            newpath = path + '/' + \
                info[filename]['OBJECT'] + '_' + \
                info[filename]['DATE-OBS'] + '_' + \
                info[filename]['CONFIG'] + '_' + \
                info[filename]['OBSID']

            if info[filename]['OBSCLASS'] == 'science':
                if newpath not in scidirs:
                    scidirs.append(newpath)

            elif info[filename]['OBSCLASS'] == 'partnerCal':
                if newpath not in teldirs:
                    teldirs.append(newpath)

            newpath += '/Intermediate'
            logger.debug('newpath: %s', newpath)
            if not os.path.exists(newpath):
                logger.debug('Creating %s', newpath)
                os.makedirs(newpath)
            logger.debug('Copying %s to %s', filename, newpath)
            shutil.copy2(rawpath + '/' + filename, newpath)

        elif info[filename]['OBSTYPE'] in ['FLAT', 'ARC', 'DARK']:                                        # CALIBRATIONS
            # Flats & Arcs may be taken in a different observation with a different Obs-ID.
            # Try to find the matching science observation(s) here by associating calibrations with
            # science targets at the same coordaintes on the same night with the same config:
            matches = {}
            for f in info.keys():
                if (info[f]['OBSCLASS'] == 'science' and
                    info[f]['DATE-OBS'] == info[filename]['DATE-OBS'] and
                    info[f]['CONFIG'] == info[filename]['CONFIG'] and
                    info[f]['COORDS'] == info[filename]['COORDS'] and
                    info[f]['OBJECT'] not in matches) or \
                        ('Pinholes' in info[filename]['SLIT'] and  # pinholes don't match the full config
                         info[f]['OBSCLASS'] == 'science' and
                         info[f]['DATE-OBS'] == info[filename]['DATE-OBS'] and
                         info[f]['CAMERA'] == info[filename]['CAMERA'] and
                         info[f]['OBJECT'] not in matches.keys()):
                    matches[info[f]['OBJECT']] = {'CONFIG': info[f]['CONFIG'], 'OBSID': info[f]['OBSID']}
                    logger.debug('match: %s: %s', info[f]['OBJECT'], matches[info[f]['OBJECT']])
            logger.debug('This calibration matches: %s', matches.keys())
            for m in matches.keys():
                newpath = path + '/' + \
                          m + '_' + \
                          info[filename]['DATE-OBS'] + '_' + \
                          matches[m]['CONFIG'] + '_' + \
                          matches[m]['OBSID'] + '/' + \
                          'Calibrations'
                logger.debug('newpath: %s', newpath)
                if newpath not in caldirs:
                    caldirs.append(newpath)

                if not os.path.exists(newpath):
                    logger.debug('Creating %s', newpath)
                    os.makedirs(newpath)
                logger.debug('Copying %s to %s', filename, newpath)
                shutil.copy2(rawpath + '/' + filename, newpath)

        else:
            logger.warning('Unhandlded file: %s', filename)

    logger.info('Creating Final directories...')
    for d in scidirs + teldirs:
        newpath = d + '/Final'
        if not os.path.exists(newpath):
            logger.debug('Creating %s', newpath)
            os.makedirs(newpath)

    logger.info('Adding directories to the config file...')

    for d in scidirs:
        logger.info('%s', d)
        config.set('ScienceDirectories', d, True)

    for d in teldirs:
        logger.info('%s', d)
        config.set('TelluricDirectories', d, True)

    for d in caldirs:
        logger.info('%s', d)
        config.set('CalibrationDirectories', d, True)

    with open(configfile, 'w') as f:
        config.write(f)

    return


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs-pype.cfg')
