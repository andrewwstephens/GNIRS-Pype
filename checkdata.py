#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ConfigParser
import datetime
import log
import os
import gnirsHeaders


def start(configfile):
    """
    For each Science Directory:
    - check that src.list and all.list exist and that they are at least the minimum size
    - check that all the files in the list exist, have the same config, coords, exptime, and target
    - find the best calibration directory and do the same list checks
    - find the closest Telluric directory and do the same list checks
    """

    logger = log.getLogger('checkdata')
    logger.info('Checking that each observation has the required calibrations.')

    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)

    caldirs = config.options('CalibrationDirectories')
    scidirs = config.options('ScienceDirectories')
    teldirs = config.options('TelluricDirectories')
    logger.debug('Calibration directories: %s', caldirs)
    logger.debug('Science directories: %s', scidirs)
    logger.debug('Telluric directories: %s', teldirs)

    if len(scidirs) < 1:
        logger.error('No science directories are listed in %s', configfile)
    if len(caldirs) < 1:
        logger.error('No calibration directories are listed in %s', configfile)
    if len(teldirs) < 1:
        logger.error('No telluric directories are listed in %s', configfile)

    for sdir in scidirs:
        if config.getboolean('ScienceDirectories', sdir):  # Only check directories marked True

            logger.info('Checking science directory %s...', sdir)
            sci_info = gnirsHeaders.start(sdir)
            checklist('all.list', path=sdir, headerdict=sci_info)
            checklist('src.list', path=sdir, headerdict=sci_info)
            if os.path.exists(sdir + '/sky.list'):  # If there is a sky.list then check it too
                checklist('sky.list', path=sdir, headerdict=sci_info)

            sfile = next(iter(sci_info))  # The sci exposures should all be the same, so just use the first

            logger.info('Searching for matching calibrations...')
            for cdir in caldirs:
                logger.debug('...%s', cdir)
                cal_info = gnirsHeaders.start(cdir)
                # The arcs and flats will match but the pinholes will not.  Find an IR flats to compare with:

                # TODO: We should probably compare ALL the calibration configs (IR, QH, arcs) to make sure they match
                
                irflats = [k for k in cal_info.keys() if
                           cal_info[k]['OBSTYPE'] == 'FLAT' and cal_info[k]['GCALLAMP'] == 'IRhigh']
                logger.debug('IR flats: %s', irflats)
                f = irflats[0]
                if cal_info[f]['DATE-OBS'] == sci_info[sfile]['DATE-OBS'] and \
                        cal_info[f]['CONFIG'] == sci_info[sfile]['CONFIG'] and \
                        cal_info[f]['COORDS'] == sci_info[sfile]['COORDS']:
                    logger.info('Calibration directory %s matches', cdir)
                    # Check that we have all the cals we need:
                    checklist('arcs.list', path=cdir, headerdict=cal_info)
                    checklist('IRflats.list', path=cdir, headerdict=cal_info)
                    checklist('QHflats.list', path=cdir, headerdict=cal_info)
                    checklist('arcs.list', path=cdir, headerdict=cal_info)
                    break  # no need to check any other calibration directories

            logger.info('Searching for matching Telluric standards...')
            dt = {}
            for tdir in teldirs:
                logger.debug('...%s', tdir)
                tel_info = gnirsHeaders.start(tdir)
                tfile = next(iter(tel_info))  # use the first Telluric file
                if tel_info[tfile]['CONFIG'] == sci_info[sfile]['CONFIG'] and \
                        tel_info[tfile]['DATE-OBS'] == sci_info[sfile]['DATE-OBS']:
                    dt[tdir] = abs(tel_info[tfile]['AVETIME'] - sci_info[sfile]['AVETIME'])
                    logger.debug('This Telluric directory matches; dt = %s', dt[tdir])

            if len(dt) == 0:
                logger.error('No matching Telluric standard found.')
            else:
                logger.info('Found %d Tellurics with the same config on the same night', len(dt))
                best = dt.keys()[dt.values().index(min(dt.values()))]
                logger.info('The best Telluric is: %s', best)
                if dt[best] > datetime.timedelta(hours=1.5):
                    logger.warning('Telluric was taken %s from the science', dt[best])
                if best != tdir:  # re-read the FITS headers if necessary
                    tel_info = gnirsHeaders.start(best)
                checklist('all.list', path=tdir, headerdict=tel_info)
                checklist('src.list', path=tdir, headerdict=tel_info)

    return

# ----------------------------------------------------------------------------------------------------------------------


def checklist(filelist, path, headerdict):
    logger = log.getLogger('checklist')

    if os.path.exists(path + '/' + filelist):
        logger.info('Found %s', filelist)
        with open(path + '/' + filelist, "r") as f:
            files = f.read().strip().split('\n')
        logger.debug('files: %s', files)

        configs = []
        coords = []
        exptimes = []
        objects = []
        obstypes = []

        for f in files:

            if f not in headerdict:  # Check that all the files in the list exist
                logger.error('%s is in %s but was not found', f, filelist)
                continue

            configs.append(headerdict[f]['CONFIG'])
            coords.append(headerdict[f]['COORDS'])
            exptimes.append(headerdict[f]['EXPTIME'])
            objects.append(headerdict[f]['OBJECT'])
            obstypes.append(headerdict[f]['OBSTYPE'])

        if headerdict[f]['OBSTYPE'] == 'OBJECT':  # Science or Telluric
            minfiles = 2
        elif headerdict[f]['OBSTYPE'] == 'ARC':
            minfiles = 1
        elif headerdict[f]['OBSTYPE'] == 'FLAT':  # IR, QH, or pinholes
            minfiles = 1
        else:
            minfiles = 1
        logger.debug('Minimum number of files: %d', minfiles)
        if len(files) < minfiles:
            logger.error('%s only has %d files', filelist, len(files))

        if len(list(set(configs))) > 1:  # Check that all the configs are the same
            logger.error('Multiple configurations: %s', configs)
        else:
            logger.debug('Configs match')

        if len(list(set(coords))) > 1:  # Check that all the coordinates are the same
            logger.error('Multiple coordinates: %s', coords)
        else:
            logger.debug('Coordinates match')

        if len(list(set(exptimes))) > 1:  # Check that all the exposure times are the same
            logger.error('Multiple exposure times: %s', exptimes)
        else:
            logger.debug('Exposure times match')

        if len(list(set(objects))) > 1:  # Check that all the target names are the same
            logger.error('Multiple target names: %s', objects)
        else:
            logger.debug('Target names match')

        if len(list(set(obstypes))) > 1:  # Check that all the obstypes are the same
            logger.error('Multiple observation types: %s', obstypes)
        else:
            logger.debug('Observation types match')

    else:
        logger.error('Could not find %s', filelist)
    return

# ----------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='INFO')
    start('gnirs.cfg')
