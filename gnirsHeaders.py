#!/usr/bin/env python

import glob, datetime, dateutil.parser
from astropy.io import fits

def start(data_directory):
    """
    Create a dictionary of all the relevant header info! :-)
    """
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

    obsids = [info.values()[i]['OBSID'] for i in range(len(info.values()))]

    for o in list(set(obsids)):
        time = []
        for f in info.keys():
            if info[f]['OBSID'] == o:
                time.append(info[f]['DATETIME'])
        info[f]['AVETIME'] = min(time) + datetime.timedelta(seconds=(max(time)-min(time)).total_seconds()/2.)
    
    return info

if __name__ == '__main__':
    start(data_directory)
