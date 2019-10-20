#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
FITS header handler
"""
from astropy.io import fits
import datetime
import dateutil.parser
import glob
import log
import os
import re


# ----------------------------------------------------------------------------------------------------------------------
def info(files_or_directory):
    """
    Return a dictionary of relevant FITS header information.
    Input may be a single FITS file, a list of FITS files, or the path to a directory of FITS files.
    """
    logger = log.getLogger('header.info')

    if isinstance(files_or_directory, list):
        files = files_or_directory
    elif isinstance(files_or_directory, basestring):
        if os.path.isfile(files_or_directory):
            files = [files_or_directory]
        elif os.path.isdir(files_or_directory):
            files = sorted(glob.glob(files_or_directory + '/N*.fits'))
        else:
            logger.error('Having trouble understanding input string')
            raise SystemExit
    else:
        logger.error('Having trouble understanding input')
        raise SystemExit

    logger.info('Reading header information...')
#    logger.debug('files: %s', files)

    data = {}
    for filename in files:
        # logger.debug('%s', filename)
        f = filename[filename.rfind('/')+1:]  # filename base to be used as output dictionary key
        data[f] = {}

        if not os.path.exists(filename):
            logger.error('Cannot find %s', filename)
            raise SystemExit

        header = fits.open(filename)[0].header

        if header['INSTRUME'] != 'GNIRS':
            logger.debug('Skipping %s which is not from GNIRS', filename)
            continue

        if 'SXD' not in header['PRISM'] and 'LXD' not in header['PRISM']:
            logger.debug('Skipping %s which is not cross-dispersed', filename)
            continue

        for key in ['PRISM', 'OBSTYPE', 'OBSID', 'OBSCLASS', 'OBJECT', 'RA', 'DEC', 'DATE-OBS', 'TIME-OBS', 'EXPTIME',
                    'CAMERA', 'DECKER', 'GRATING', 'SLIT', 'GRATWAVE', 'POFFSET', 'QOFFSET', 'GCALLAMP']:
            try:
                data[f][key] = header[key].strip() if isinstance(header[key], str) else header[key]
            except:
                logger.debug('%s[%s] = None', f, key)
                data[f][key] = None

        data[f]['DATETIME'] = dateutil.parser.parse(data[f]['DATE-OBS'] + ' ' + data[f]['TIME-OBS'])
        data[f]['COORDS'] = '%.2f %.2f' % (data[f]['RA'], data[f]['DEC'])
        data[f]['OBJECT'] = re.sub('[^a-zA-Z0-9]', '', data[f]['OBJECT'])  # replace non-alphanumeric characters
        data[f]['DATE-OBS'] = data[f]['DATE-OBS'].replace('-', '')
        # info[f]['OBSID'] = info[f]['OBSID'].replace('-', '_')
        data[f]['CONFIG'] = \
            re.sub('(ong|hort|lue|ed)', '', data[f]['CAMERA'][:data[f]['CAMERA'].find('_')]) + '_' + \
            data[f]['PRISM'][data[f]['PRISM'].find('+')+1:data[f]['PRISM'].find('_')] + '_' + \
            data[f]['GRATING'][:data[f]['GRATING'].find('/')] + '_' + \
            data[f]['SLIT'][:data[f]['SLIT'].find('_')] + '_' + \
            '%6.4fum' % data[f]['GRATWAVE']

    obsids = [data.values()[i]['OBSID'] for i in range(len(data.values()))]

    for o in list(set(obsids)):  # Add the average time for each Obs-ID
        time = []
        for f in data.keys():
            if data[f]['OBSID'] == o:
                time.append(data[f]['DATETIME'])
        avetime = min(time) + datetime.timedelta(seconds=(max(time)-min(time)).total_seconds()/2.)
        for f in data.keys():
            if data[f]['OBSID'] == o:
                data[f]['AVETIME'] = avetime

#    logger.debug('%s', data)
    return data


if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    info("rawData")  # Directory
    # info('N20110516S0151.fits')  # single file
    # info(['N20110516S0155.fits', 'N20110516S0156.fits', 'N20110516S0157.fits', 'N20110516S0159.fits'])  # List

