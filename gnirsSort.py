#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ConfigParser
import datetime
import log
import os
import shutil
import gnirsHeaders


def start(configfile):
    """
    Sort and copy the GNIRS raw data to subdirectories: target / date / config.
    """
    logger = log.getLogger('Sort')

    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)
    rawpath = config.get('getDataConfig', 'rawPath')

    info = gnirsHeaders.start(rawpath)

    path = os.getcwd()
    logger.info("Path to raw data: %s", rawpath)
    logger.info('Sorted data will be copied to %s', path)
    logger.info("Creating new directories and copying files...")

    scidirs = []     # List of output science directories
    teldirs = []     # List of output Telluric standard directories
    acq = {}         # Acqusition parameters OBJECT:{POFFSET, QOFFSET}

    for filename in sorted(info.keys()):

        logger.debug('%s: OBSCLASS:%s OBSTYPE:%s OBJECT:%s SLIT:%s GCALLAMP:%s',
                     filename, info[filename]['OBSCLASS'], info[filename]['OBSTYPE'], info[filename]['OBJECT'],
                     info[filename]['SLIT'], info[filename]['GCALLAMP'])

        if info[filename]['OBSTYPE'] == 'OBJECT' and info[filename]['OBSCLASS'] == 'acq':

            acq[info[filename]['OBJECT']] = {'POFFSET': info[filename]['POFFSET'], 'QOFFSET': info[filename]['QOFFSET']}
            logger.debug('Acquisition parmameters: %s', acq[info[filename]['OBJECT']])

        elif info[filename]['OBSTYPE'] == 'OBJECT' and info[filename]['OBSCLASS'] == 'science':

            newpath = path + '/' + \
                info[filename]['OBJECT'] + '/' + \
                info[filename]['DATE-OBS'] + '/' + \
                info[filename]['CONFIG'] + '/' + \
                'Science/' + \
                info[filename]['OBSID']
            logger.debug('newpath: %s', newpath)
            if newpath not in scidirs:
                scidirs.append(newpath)
            if not os.path.exists(newpath):
                logger.debug('Creating %s', newpath)
                os.makedirs(newpath)
            logger.debug('Copying %s to %s', filename, newpath)
            shutil.copy2(rawpath + '/' + filename, newpath)

            if inslit(slit=info[filename]['SLIT'], decker=info[filename]['DECKER'],
                      p=info[filename]['POFFSET'] - acq[info[filename]['OBJECT']]['POFFSET'],
                      q=info[filename]['QOFFSET'] - acq[info[filename]['OBJECT']]['QOFFSET']):
                append(filename, newpath + '/src.list')
            else:
                append(filename, newpath + '/sky.list')

        elif info[filename]['OBSTYPE'] == 'OBJECT' and info[filename]['OBSCLASS'] == 'partnerCal':

            # Associate this Telluric with all the science targets with the same config within 1.5 hours:
            matches = []
            for f in info.keys():
                if info[f]['OBSCLASS'] == 'science' and \
                        info[f]['CONFIG'] == info[filename]['CONFIG'] and \
                        abs(info[f]['AVETIME'] - info[filename]['AVETIME']) < datetime.timedelta(hours=1.5) and \
                        info[f]['OBJECT'] not in matches:
                    matches.append(info[f]['OBJECT'])
            logger.debug('This Telluric matches: %s', matches)
            for m in matches:
                newpath = path + '/' + m + '/' + \
                    info[filename]['DATE-OBS'] + '/' + \
                    info[filename]['CONFIG'] + '/' + \
                    'Tellurics/' + \
                    info[filename]['OBSID']
                logger.debug('newpath: %s', newpath)
                if newpath not in teldirs:
                    teldirs.append(newpath)
                if not os.path.exists(newpath):
                    logger.debug('Creating %s', newpath)
                    os.makedirs(newpath)
                logger.debug('Copying %s to %s', filename, newpath)
                shutil.copy2(rawpath + '/' + filename, newpath)
                append(filename, newpath + '/src.list')

        elif info[filename]['OBSTYPE'] in ['FLAT', 'ARC', 'DARK']:

            # Associate these calibrations with all the science targets on the same night with the same config
            matches = {}
            for f in info.keys():
                if (info[f]['OBSCLASS'] == 'science' and
                    info[f]['CONFIG'] == info[filename]['CONFIG'] and
                    info[f]['DATE-OBS'] == info[filename]['DATE-OBS'] and
                    info[f]['OBJECT'] not in matches) or \
                        ('Pinholes' in info[filename]['SLIT'] and  # pinholes don't match the full config
                         info[f]['OBSCLASS'] == 'science' and
                         info[f]['CAMERA'] == info[filename]['CAMERA'] and
                         info[f]['DATE-OBS'] == info[filename]['DATE-OBS'] and
                         info[f]['OBJECT'] not in matches.keys()):
                    matches[info[f]['OBJECT']] = info[f]['CONFIG']
            logger.debug('This calibration matches: %s', matches.keys())

            for m in matches.keys():
                newpath = path + '/' + m + '/' + \
                    info[filename]['DATE-OBS'] + '/' + \
                    matches[m] + '/' + \
                    'Calibrations'
                logger.debug('newpath: %s', newpath)

                if not os.path.exists(newpath):
                    logger.debug('Creating %s', newpath)
                    os.makedirs(newpath)
                logger.debug('Copying %s to %s', filename, newpath)
                shutil.copy2(rawpath + '/' + filename, newpath)

                if info[filename]['OBSTYPE'] == 'ARC':
                    append(filename, newpath + '/arcs.list')

                elif info[filename]['OBSTYPE'] == 'DARK':
                    append(filename, newpath + '/darks.list')

                elif info[filename]['OBSTYPE'] == 'FLAT' and \
                        info[filename]['GCALLAMP'] == 'QH' and \
                        'Pinholes' in info[filename]['SLIT']:
                    append(filename, newpath + '/pinholes.list')

                elif info[filename]['OBSTYPE'] == 'FLAT' and \
                        info[filename]['GCALLAMP'] == 'QH'and \
                        'Pinholes' not in info[filename]['SLIT']:
                    append(filename, newpath + '/QHflats.list')

                elif info[filename]['OBSTYPE'] == 'FLAT' and \
                        info[filename]['GCALLAMP'] == 'IRhigh':
                    append(filename, newpath + '/IRflats.list')

                else:
                    logger.error('Unknown configuration for %s', filename)
        else:
            logger.warning('Unhandlded file: %s', filename)

    logger.info('Updating config file...')

    for d in scidirs:
        logger.debug('Setting %s = True', d)
        config.set('ScienceDirectories', d, True)

    for d in teldirs:
        logger.debug('Setting %s = True', d)
        config.set('TelluricDirectories', d, True)

    with open(configfile, 'w') as f:
        config.write(f)

    return

# ----------------------------------------------------------------------------------------------------------------------


def append(string, filename):
    logger = log.getLogger('append')
    with open(filename, 'a+') as f:
        if string not in f.read().split('\n'):
            logger.debug('Adding %s to %s', string, filename)
            f.write(string + '\n')

# ----------------------------------------------------------------------------------------------------------------------


def inslit(slit, decker, p, q):
    logger = log.getLogger('inslit')
    logger.debug('Offsets:  %.1f  %.1f arcsec', p, q)
    if 'SC' in decker and 'XD' in decker:
        slitlength = 7.0
    elif 'LC' in decker and 'XD' in decker:
        slitlength = 5.0
    else:
        logger.error('cannot calculate slit length')
        raise SystemExit
    slitwidth = float(slit[:slit.find('arcsec')])
    logger.debug('Slit size: %.1f x %.1f arcsec', slitwidth, slitlength)
    if abs(p) < slitwidth / 2.0 and abs(q) < slitlength / 2.0:
        return True
    else:
        return False

# ----------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
