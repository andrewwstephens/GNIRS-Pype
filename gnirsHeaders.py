#!/usr/bin/env python

import datetime
import dateutil.parser
import glob
import log
from astropy.io import fits

def start(data_directory):
    """
    Create a dictionary of all the relevant header info! :-)
    """
    logger = log.getLogger('Headers')
    info = {}
    for filename in sorted(glob.glob(data_directory + '/N*.fits')):
        f = filename[filename.rfind('/')+1:]
        info[f] = {}
        header = fits.open(filename)[0].header
        for key in ['INSTRUME', 'PRISM', 'OBSTYPE', 'OBSID', 'OBSCLASS', 'OBJECT', 'DATE-OBS', 'TIME-OBS', 'CAMERA',
                    'DECKER', 'GRATING', 'SLIT', 'GRATWAVE', 'POFFSET', 'QOFFSET', 'GCALLAMP']:
            try:
                info[f][key] = header[key].strip() if isinstance(header[key],str) else header[key]
            except:
                info[f][key] = None
                
        info[f]['DATETIME'] = dateutil.parser.parse(info[f]['DATE-OBS'] + ' ' + info[f]['TIME-OBS'])

        info[f]['CONFIG'] = info[f]['CAMERA'][:info[f]['CAMERA'].find('_')]+'_' + \
                 info[f]['PRISM'][info[f]['PRISM'].find('+')+1:info[f]['PRISM'].find('_')] + '_' + \
                 info[f]['GRATING'][:info[f]['GRATING'].find('/')] + '_' + \
                 info[f]['SLIT'][:info[f]['SLIT'].find('_')] + '_' + \
                 '%6.4fum' % info[f]['GRATWAVE']

    obsids = [info.values()[i]['OBSID'] for i in range(len(info.values()))]

    for o in list(set(obsids)): # Add the average time for each Obs-ID
        time = []
        for f in info.keys():
            if info[f]['OBSID'] == o:
                time.append(info[f]['DATETIME'])
        avetime = min(time) + datetime.timedelta(seconds=(max(time)-min(time)).total_seconds()/2.)
        for f in info.keys():
            if info[f]['OBSID'] == o:
                info[f]['AVETIME'] = avetime

    logger.debug('%s', info)
    return info

if __name__ == '__main__':
    start("rawData")
