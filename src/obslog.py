#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Obslog handler

This library provides methods to:
- download the raw obslog.txt from the archive
- parse the raw obslog.txt
- extract the relevant obslog data for a particular Obs-ID
- output a csv file for the requested Obs-ID with the original header info plus (P,Q) corrected for acq offsets
- read an obslog.csv and return a dictionary of the relevant FITS header keywords (all the FITS header keywords?)
We may want to provide some stand-alone obslog handling functionality.
"""
import csv
import gnirsHeaders
import log
import os
import re
import urllib


# ----------------------------------------------------------------------------------------------------------------------
def download(progid, date):
    """
    Download an observing log from the Gemini archive.
    :param progid: string program ID
    :param date: string date, format YYYYMMDD
    """
    logger = log.getLogger('obslog.download')

    if not re.match(r'^G[NS]-[0-9]{4}[AB]-([CQ]|DD|FT|LP|SV)-[0-9]{0,3}$', progid):
        logger.error('This does not look like a program ID: %s', progid)
        raise SystemExit

    obslog = date + '_' + progid + '_obslog.txt'
    url = 'https://archive.gemini.edu/file/' + obslog
    logger.debug('URL: %s', url)
    logger.info('Downloading %s...', obslog)
    urllib.urlretrieve(url, obslog)
    return


# ----------------------------------------------------------------------------------------------------------------------
def writecsv(obsid, date, output='obslog.csv', rawpath=None):
    """
    Extract the data for the requested Obs-ID from the obslog text file
    Find the last acquisition image before the requested obsid
    Extract the P,Q of that acquisition image
    Add new keys for the absolute P and Q offsets (header offsets - acquisition offsets)
    Write the header info to a csv file
    :param obsid: string observation ID
    :param date: string date, format YYYYMMDD
    :param rawpath: string path to raw data and obslog text file
    :param output: string output file name
    """
    logger = log.getLogger('obslog.writecsv')
    progid = obsid[:obsid.rfind('-')]
    logger.debug('Program ID: %s', progid)
    obslog = date + '_' + progid + '_obslog.txt'
    if rawpath:
        obslog = rawpath + '/' + obslog
    data = readtxt(obslog)

    output_data = {}
    first_spectrum = None
    files = sorted(data.keys(), reverse=True)  # Note the reverse sort here
    for f in files:  # Go through the whole list in case there were interruptions for re-aqcuisitions
        if data[f]['Observation ID'] == obsid:
            output_data[f] = data[f]
            logger.debug('Including %s', f)
            first_spectrum = f
    logger.debug('First spectrum: %s', first_spectrum)

    last_acq = None
    for i in range(files.index(first_spectrum)+1, len(files)):  # again, files is reverse sorted
        if data[files[i]]['ACQ'] == 'Y':
            last_acq = files[i]
            logger.info('Last acqusition file: %s', last_acq)
            break

    # Get the header info for the requested images plus the last acquisition image:
    if rawpath:
        fitsfiles = [rawpath + '/' + f for f in ([last_acq] + sorted(output_data.keys()))]
    else:
        fitsfiles = [last_acq] + sorted(output_data.keys())
    headerinfo = gnirsHeaders.info(fitsfiles)

    for f in output_data.keys():  # Add new keys for the absolute P and Q offsets:
        headerinfo[f]['P'] = headerinfo[f]['POFFSET'] - headerinfo[last_acq]['POFFSET']
        headerinfo[f]['Q'] = headerinfo[f]['QOFFSET'] - headerinfo[last_acq]['QOFFSET']
    logger.debug('Updated Info: %s', headerinfo)

    def mergedict(a, b):
        a.update(b)
        return a

    logger.info('Writing %s...', output)  # Write the info for the requested Obs-ID into a csv file:
    with open(output, mode='w') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=['FITSFILE'] + headerinfo[f].keys())
        writer.writeheader()
        for k, d in sorted(headerinfo.items()):
            writer.writerow(mergedict({'FITSFILE': k}, d))

    return


# ----------------------------------------------------------------------------------------------------------------------
def readcsv(csvfile):
    """
    Read an obslog csv file
    :param csvfile: text csv obslog filename
    :return: dictionary of obslog info
    """
    logger = log.getLogger('obslog.readcsv')

    if not os.path.exists(csvfile):
        logger.error('Cannot access %s', csvfile)
        raise SystemExit

    data = {}
    with open(csvfile, mode='r') as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            data[row['FITSFILE']] = row
    logger.debug('Data: %s', data)
    return data


# ----------------------------------------------------------------------------------------------------------------------
def readtxt(obslog):
    """
    Parse the raw obslog .txt file from the archive.
    :param obslog: string obslog file name
    The archived obslog file is in fixed-witdth format, but the column width changes between each obslog according to
    the width of the values in the table.  There are also free-form comments from the observer interspersed in the file
    with no indication of which lines are comments.  This parser uses the location of the column headers to set the
    width of the fixed columns and sets patterm match guidelines for the contents of each column.  If any of the parsed
    columns do not match the expected pattern the whole line is assumed to be a comment.
    NOTE: Many old obslogs (before 2012-Feb-10) are corrupt and will need to be regenerated by Gemini.
    """

    logger = log.getLogger('obslog.readtxt')

    if not os.path.exists(obslog):
        logger.error('Cannot access %s', obslog)
        raise SystemExit

    logger.info('Reading %s', obslog)

    with open(obslog) as f:  # Since we will have to go through the data twice, read the whole file at once.
        data = f.readlines()

    header = ['Observation ID', 'Data Labels', 'File Numbers', 'Dataset UT', 'Target Name', 'Filters',  'Slit',
              'Grating/Wavelength', 'Camera/Prism', 'ExpTime/LNR/Coadds', 'ACQ']

    pattern = dict()  # Enforce formatting rules to avoid parsing comments as data:
    pattern['Observation ID'] = re.compile(r'^G[NS]-[0-9]{4}[AB]-([CQ]|DD|FT|LP|SV)-[0-9]{0,3}-[0-9]+$')
    pattern['Data Labels'] = re.compile(r'[0-9]+-*[0-9]*')                      # 1, 2-3, 45-67, 890-1234
    pattern['File Numbers'] = re.compile(r'[0-9]+-*[0-9]*')                     # 1, 2-3, 45-67, 890-1234
    pattern['Dataset UT'] = re.compile(r'^[0-9]{2}:[0-9]{2}:[0-9]{2}$')         # 09:58:15
    pattern['Target Name'] = re.compile(r'[a-zA-Z0-9_-]+')                      # Match any string
    pattern['Filters'] = re.compile(r'[A-Z0-9\-]+')                             # H, XD, H2, X, J, H
    pattern['Slit'] = re.compile(r'[a-zA-Z0-9]+')                               # 0.675, ACQ, LgPin
    pattern['Grating/Wavelength'] = re.compile(r'[0-9]{2,3}/[0-9]\.[0-9]{2}')   # 32/1.65, 111/1.68
    pattern['Camera/Prism'] = re.compile(r'[A-Z]{2}/[A-Z]{3}')                  # LB/MIR, SB/SXD
    pattern['ExpTime/LNR/Coadds'] = re.compile(r'[0-9]+\.[0-9]/[0-9]+/[0-9]+')  # 0.2/1/25, 300.0/32/1
    pattern['ACQ'] = re.compile(r'^Y*$')                                        # Y or ''

    indx = {}
    for line in data:
        if 'Electronic Observing Log' in line:
            date = line.split()[-1][7:]
            logger.debug('Log date: %s', date)
        if line[0:14] == 'Observation ID':   # This defines the start of the header row
            for h in header:
                indx[h] = line.find(h)        # Find where each column starts
            break                             # No need to go farther

    width = {}                        # Find the width of each row
    for i in range(len(header) - 1):  # This requires that 'header' be an ordered array (not a dictionary)
        width[header[i]] = indx[header[i + 1]] - indx[header[i]]
    width[header[i+1]] = 1            # The ACQ field is either 'Y' or blank

    val = {}
    match = {}
    info = {}
    for line in data:
        logger.debug('\n%s', line)
        files = []
        for h in header:
            val[h] = line[indx[h]: indx[h] + width[h]].strip()
            match[h] = re.match(pattern[h], val[h])
            logger.debug('%s: "%s"   %s' % (h, val[h], match[h]))

        # Maybe throw a warning if only match 1 fails; indicating a likely bad pattern specification?

        if None in match.values():
            logger.debug('Failed to match all patterns  ->  This is a comment')
            continue

        if '-' in val['File Numbers']:
            start, stop = val['File Numbers'].split('-')
            for i in range(int(start), int(stop)+1):
                files.append(i)
        else:
            files.append(int(val['File Numbers']))

        for filenum in files:
            f = 'N%sS%04d.fits' % (date, filenum)
            logger.debug('File: %s', f)
            info[f] = {}
            for h in [header[0]] + header[3:]:  # Skip 'Data Labels' and "File Numbers'
                info[f][h] = val[h]

    logger.debug('info: %s', info)
    return info


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    # readtxt('rawData/20110516_GN-2011A-Q-126_obslog.txt')
    # writecsv('GN-2011A-Q-126-6', '20110516', rawpath='/home/astephens/github/GNIRS_PIPELINE/rawData')  # Science
    writecsv('GN-2011A-Q-126-13', '20110516', rawpath='/home/astephens/github/GNIRS_PIPELINE/rawData')  # Telluric
    # readcsv('obslog.csv')
    # download(progid='GN-2018A-FT-112', date='20180729')

