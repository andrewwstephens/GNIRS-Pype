#!/usr/bin/env python
"""
Associate and link GNIRS calibrations
"""
import ConfigParser
import datetime
import header
import log
import numpy
import os


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
    # What we really want is the average time of the files in each directory (sci, tel, cal), so
    # recalculate the average times using only the files in each directory (so without mixing science and cals).

    alldirs = []   # List of all output directories
    caldirs = []   # List of calibration directories
    scidirs = []   # List of science directories
    teldirs = []   # List of Telluric standard directories
    obstype = []   # list of observation types (science, telluric, calibration)
    obstime = []   # List of average observation times in each directory
    obsconf = []   # List of observation configurations

    for path, process in scipath:
        logger.debug('%s = %s', path, process)
        if process:
            info = header.info(path + '/Intermediate')
            first_file = sorted(info.keys())[0]
            logger.debug('first_file = %s', first_file)
            alldirs.append(path)
            scidirs.append(path)
            obsconf.append(info[first_file]['CONFIG'])
            obstype.append('science')
            obstime.append(info[first_file]['AVETIME'])

    for path, process in telpath:
        logger.debug('%s = %s', path, process)
        if process:
            info = header.info(path + '/Intermediate')
            first_file = sorted(info.keys())[0]
            logger.debug('first_file = %s', first_file)
            alldirs.append(path)
            teldirs.append(path)
            obsconf.append(info[first_file]['CONFIG'])
            obstype.append('telluric')
            obstime.append(info[first_file]['AVETIME'])

    for path, process in calpath:
        logger.debug('%s = %s', path, process)
        if process:
            info = header.info(path)
            first_file = sorted(info.keys())[0]  # use the first file to avoid daytime pinholes
            logger.debug('first_file = %s', first_file)
            alldirs.append(path)
            caldirs.append(path)
            obsconf.append(info[first_file]['CONFIG'])
            obstype.append('calibration')
            obstime.append(info[first_file]['AVETIME'])

    obstime = numpy.array(obstime)
    obstype = numpy.array(obstype)
    obsconf = numpy.array(obsconf)

    logger.debug('alldirs: %s', alldirs)

    logger.info('------------------------------------------------------------------')
    logger.info('Observation Summary')  # Print a summary table for humans to decide if the best choices were made
    for i in range(len(alldirs)):
        logger.info('%s %s %s %s', alldirs[i], obstype[i], obsconf[i], obstime[i])
    logger.info('------------------------------------------------------------------')

    logger.info('Searching for best Telluric standard for each science observation...')
    for d in scidirs:
        logger.info('Science: %s', d)
        indx = alldirs.index(d)
        sciconf = obsconf[indx]
        scitime = obstime[indx]
        match = numpy.where((obstype == 'telluric') & (obsconf == sciconf))
        logger.debug('match: %s', match)
        dt = obstime[match] - scitime
        logger.debug('dt = %s', dt)
        imin = numpy.argmin(dt)  # the index of the min dt
        logger.debug('imin: %s', imin)

        indx = match[imin][0]   # the index of the best observation in the list of all observations
        logger.debug('indx: %s', indx)

        logger.info('Best Telluric: %s (dt = %s)', alldirs[indx], dt[imin])
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
        logger.info('Linking: %s -> %s', alldirs[indx], dest)
        os.symlink(alldirs[indx], dest)

    logger.info('Searching for best calibrations for each Telluric standard...')
    for d in teldirs:
        logger.info('Telluric: %s', d)
        indx = alldirs.index(d)
        telconf = obsconf[indx]
        teltime = obstime[indx]
        match = numpy.where((obstype == 'calibration') & (obsconf == telconf))
        dt = obstime[match] - teltime
        logger.debug('dt = %s', dt)

        imin = numpy.argmin(dt)  # the index of the min dt
        logger.debug('imin: %s', imin)

        indx = match[numpy.argmin(dt)][0]  # the index of the best observation in the list of all observations
        logger.debug('indx = %s', indx)

        logger.info('Best Calibration: %s (dt = %s)', alldirs[indx], dt[imin])
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
        logger.info('Linking %s -> %s', alldirs[indx], dest)
        os.symlink(alldirs[indx], dest)

    return


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs-pype.cfg')
