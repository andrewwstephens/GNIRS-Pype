#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Sort GNIRS data into directories

This function should be split into pieces to facilitate testing and validation:
1. Copy files
2. Create file lists
3. Match Tellurics to science and match cals to Tellurics and create links
"""
import ConfigParser
import datetime
import log
import numpy
import obslog
import os
import re
import shutil
import gnirsHeaders


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

    1. Copy the FITS files where they need to be.
    2. Go through the list of science and Telluric ObsIDs and call 'writecsv' to put in each folder.
    3. Go through the list of science ObsIDs and link to the best Telluric
    4. Go through the list of Telluric ObsIDs and link the best Calibrations

    """
    logger = log.getLogger('gnirsSort.start')

    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make config file options case-sensitive
    config.read(configfile)
    rawpath = config.get('getData', 'rawPath')

    info = gnirsHeaders.info(rawpath)

    path = os.getcwd()
    logger.info("Path to raw data: %s", rawpath)
    logger.info('Sorted data will be copied to %s', path)
    logger.info("Creating new directories and copying files...")

    alldirs = []   # List of all output directories
    caldirs = []   # List of output calibration directories
    scidirs = []   # List of output science directories
    teldirs = []   # List of output Telluric standard directories
    obstype = []   # list of observation types (science, telluric, calibration)
    obstime = []   # List of average observation times

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
                newpath += '/Intermediate'
                if newpath not in scidirs:
                    alldirs.append(newpath)
                    scidirs.append(newpath)
                    obstype.append('science')

            elif info[filename]['OBSCLASS'] == 'partnerCal':
                if newpath not in teldirs:
                    alldirs.append(newpath)
                    teldirs.append(newpath)
                    obstype.append('telluric')

            logger.debug('newpath: %s', newpath)
            if not os.path.exists(newpath):
                logger.debug('Creating %s', newpath)
                os.makedirs(newpath)
            logger.debug('Copying %s to %s', filename, newpath)
            shutil.copy2(rawpath + '/' + filename, newpath)
            append(filename, newpath + '/all.list')

        elif info[filename]['OBSTYPE'] in ['FLAT', 'ARC', 'DARK']:                                        # CALIBRATIONS
            # Flats & Arcs may be taken in a different observation (with a different Obs-ID).
            # Try to find the matching science observation(s) here by associating calibrations with
            # all the science targets at the same coordaintes on the same night with the same config
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
                    alldirs.append(newpath)
                    caldirs.append(newpath)
                    obstype.append('calibration')

                if not os.path.exists(newpath):
                    logger.debug('Creating %s', newpath)
                    os.makedirs(newpath)
                logger.debug('Copying %s to %s', filename, newpath)
                shutil.copy2(rawpath + '/' + filename, newpath)
                append(filename, newpath + '/all.list')

                if info[filename]['OBSTYPE'] == 'ARC':
                    append(filename, newpath + '/arcs.list')

                elif info[filename]['OBSTYPE'] == 'DARK':
                    append(filename, newpath + '/darks.list')

                elif info[filename]['OBSTYPE'] == 'FLAT' and \
                        info[filename]['GCALLAMP'] == 'QH' and \
                        'Pinholes' in info[filename]['SLIT']:
                    append(filename, newpath + '/pinholes.list')

                elif info[filename]['OBSTYPE'] == 'FLAT' and \
                        info[filename]['GCALLAMP'] == 'QH' and \
                        'Pinholes' not in info[filename]['SLIT']:
                    append(filename, newpath + '/QHflats.list')

                elif info[filename]['OBSTYPE'] == 'FLAT' and \
                        info[filename]['GCALLAMP'] == 'IRhigh':
                    append(filename, newpath + '/IRflats.list')

                else:
                    logger.error('Unknown configuration for %s', filename)
        else:
            logger.warning('Unhandlded file: %s', filename)

    for d in scidirs + teldirs:  # Add an obslog.csv to each folder
        chunk = re.split('[/_]', d)
        logger.debug('chunk: %s', chunk)
        if 'Intermediate' in chunk:      # Tellurics have a shorter path than science
            obslog.writecsv(obsid=chunk[-2], date=chunk[-8], rawpath=rawpath, output=d+'/obslog.csv')
        else:
            obslog.writecsv(obsid=chunk[-1], date=chunk[-7], rawpath=rawpath, output=d+'/obslog.csv')

    # This is complicated by the fact that the calibrations might not have the same ObsID as the science.
    # Also, the "AVETIME" is the average over all files with the same OBSID, so highly weighted towards the science.
    # What we really want is the average time of the files in each directory (sci, tel, cal).

    # Okay, this is probably overkill.  We could have recorded the 'AVETIME' above when we were copying files, however,
    # there the average time included science for the cals (and cals for the science).  Here we go back through and
    # re-calculate the average time using only the files in each directory (so without mixing science and cals).
    for d in alldirs:
        info = gnirsHeaders.info(d)
        firstfile = sorted(info.keys())[0]  # use the first file to avoid daytime pinholes
        logger.debug('firstfile = %s', firstfile)
        obstime.append(info[firstfile]['AVETIME'])

    obstime = numpy.array(obstime)
    obstype = numpy.array(obstype)

    logger.info('Searching for best Telluric standards...')
    for d in scidirs:
        logger.info('Science: %s', d)
        scitime = obstime[alldirs.index(d)]
        logger.debug('scitime = %s', scitime)
        dt = obstime[obstype == 'telluric'] - scitime
        logger.debug('dt = %s', dt)
        imin = numpy.argmin(dt)
        logger.debug('imin = %s', imin)
        logger.info('Best Telluric: %s (dt: %s)', teldirs[imin], str(dt[imin]))
        if abs(dt[imin] > datetime.timedelta(hours=1.5)):
            logger.warning('Time difference > 1.5 hours')
        dest = d[:d.rfind('/')] + '/Telluric'
        logger.debug('Creating link to %s', dest)
        os.symlink(teldirs[imin], dest)

    logger.info('Searching for best Telluric calibrations...')
    for d in teldirs:
        logger.info('Telluric: %s', d)
        teltime = obstime[alldirs.index(d)]
        logger.debug('teltime = %s', teltime)
        dt = obstime[obstype == 'calibration'] - teltime
        logger.debug('dt = %s', dt)
        imin = numpy.argmin(dt)
        logger.debug('imin = %s', imin)
        logger.info('Best Calibration: %s (dt: %s)', caldirs[imin], str(dt[imin]))
        if abs(dt[imin] > datetime.timedelta(hours=1)):
            logger.warning('Time difference > 1 hour')
        dest = d + '/Calibration'
        logger.debug('Creating link to %s', dest)
        os.symlink(caldirs[imin], dest)
        print('done')

    # TODO: Create A.list and B.list

    logger.info('Updating config file...')

    for d in scidirs:
        logger.debug('Setting %s = True', d)
        config.set('ScienceDirectories', d, True)

    for d in teldirs:
        logger.debug('Setting %s = True', d)
        config.set('TelluricDirectories', d, True)

    for d in caldirs:
        logger.debug('Setting %s = True', d)
        config.set('CalibrationDirectories', d, True)

    with open(configfile, 'w') as f:
        config.write(f)

    return


# ----------------------------------------------------------------------------------------------------------------------
def append(string, filename):
    logger = log.getLogger('gnirsSort.append')
    with open(filename, 'a+') as f:
        if string not in f.read().split('\n'):
            logger.debug('Adding %s to %s', string, filename)
            f.write(string + '\n')


# ----------------------------------------------------------------------------------------------------------------------
def inslit(slit, decker, p, q):
    logger = log.getLogger('gnirsSort.inslit')
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
