#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Associate and link GNIRS calibrations
"""
import ConfigParser
import datetime
import log
import numpy
import os
import gnirsHeaders


# ----------------------------------------------------------------------------------------------------------------------
def start(configfile):
    """
    :param configfile:
    :return:

    1. Go through the list of science ObsIDs and link to the best Telluric
    2. Go through the list of Telluric ObsIDs and link the best Calibrations
    """
    logger = log.getLogger('link_cals')

    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make config file options case-sensitive
    config.read(configfile)

    overwrite = config.get('defaults', 'overwrite')
    scipath = config.items("ScienceDirectories")
    telpath = config.items("TelluricDirectories")
    calpath = config.items("CalibrationDirectories")

    # Survey all the directories and record the average time when each was observed.
    # This is complicated because calibrations might not have the same ObsID as the science.
    # Also, the "AVETIME" is the average over all files with the same OBSID, so highly weighted towards the science.
    # What we really want is the average time of the files in each directory (sci, tel, cal).
    # Re-calculate the average time using only the files in each directory (so without mixing science and cals).

    alldirs = []   # List of all output directories
    caldirs = []   # List of calibration directories
    scidirs = []   # List of science directories
    teldirs = []   # List of Telluric standard directories
    obstype = []   # list of observation types (science, telluric, calibration)
    obstime = []   # List of average observation times in each directory

    for path, process in scipath:
        logger.debug('%s = %s', path, process)
        if process:
            header = gnirsHeaders.info(path + '/Intermediate')
            first_file = sorted(header.keys())[0]
            logger.debug('first_file = %s', first_file)
            alldirs.append(path)
            scidirs.append(path)
            obstype.append('science')
            obstime.append(header[first_file]['AVETIME'])

    for path, process in telpath:
        logger.debug('%s = %s', path, process)
        if process:
            header = gnirsHeaders.info(path + '/Intermediate')
            first_file = sorted(header.keys())[0]
            logger.debug('first_file = %s', first_file)
            alldirs.append(path)
            teldirs.append(path)
            obstype.append('telluric')
            obstime.append(header[first_file]['AVETIME'])

    for path, process in calpath:
        logger.debug('%s = %s', path, process)
        if process:
            header = gnirsHeaders.info(path)
            first_file = sorted(header.keys())[0]  # use the first file to avoid daytime pinholes
            logger.debug('first_file = %s', first_file)
            alldirs.append(path)
            caldirs.append(path)
            obstype.append('calibration')
            obstime.append(header[first_file]['AVETIME'])

    obstime = numpy.array(obstime)
    obstype = numpy.array(obstype)

    logger.debug('alldirs: %s', alldirs)

    logger.info('------------------------------------------')
    logger.info('Calibration Summary')  # Print a summary table for humans to decide if the best choices were made
    for i in range(len(alldirs)):
        logger.info('%s %s %s', alldirs[i], obstype[i], obstime[i])
    logger.info('------------------------------------------')

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
        dest = d + '/Telluric'
        if os.path.exists(dest):
            if overwrite:
                logger.warning('Removing old link: %s', dest)
                os.remove(dest)
            else:
                logger.warning('Cannot create link %s', dest)
                logger.warning('Link already exists and overwrite = False')
                continue
        logger.info('Linking: %s -> %s', teldirs[imin], dest)
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
        dest = d + '/Calibrations'
        if os.path.exists(dest):
            if overwrite:
                logger.warning('Removing old link: %s', dest)
                os.remove(dest)
            else:
                logger.warning('Cannot create link %s', dest)
                logger.warning('Link already exists and overwrite = False')
                continue
        logger.info('Linking %s -> %s', caldirs[imin], dest)
        os.symlink(caldirs[imin], dest)

    return


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
