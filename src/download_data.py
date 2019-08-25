#!/usr/bin/env python

import ConfigParser
import glob
import log
from ndmapperDownloader import download_query_gemini
import os


def start(configfile):
    """
    Get raw GNIRS data.
    :param configfile: name of GNIRS-pype configuration file [string]
    :return:
    """
    logger = log.getLogger('download_data')

    logger.info('####################################')
    logger.info('#                                  #')
    logger.info('#      Download GNIRS data         #')
    logger.info('#                                  #')
    logger.info('####################################')

    config = ConfigParser.SafeConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)

    rawdatadir = config.get('getData', 'RawDataDir')
    program = config.get('getData', 'Program')
    logger.debug('rawdata: %s', rawdatadir)
    logger.debug('program: %s', program)

    if program == 'None':
        logger.debug('No data to download')
        return

    if not os.path.exists(rawdatadir):
        try:
            os.mkdir(rawdatadir)
        except Exception, e:
            logger.error('Failed to create %s:\n%s', rawdatadir, str(e))

    rawfiles = glob.glob(rawdatadir + '*.fits')
    logger.debug('rawfiles: %s', rawfiles)
    if len(rawfiles) > 0:
        logger.warning('There are already FITS files in %s', rawdatadir)
        logger.warning('Continue with download of all the raw data for %s ? (y/[n])', program)
        answer = raw_input('')
        if 'n' in answer.lower():
            logger.info('Skipping the data download.')
            logger.info('Set Program = None in your config file to avoid this in the future.')
            return

    url = 'https://archive.gemini.edu/download/' + program + '/notengineering/NotFail/present/canonical'
    logger.info('Downloading data from Gemini public archive (this may take a few minutes)...')
    logger.debug('URL: %s', url)
    download_query_gemini(url, rawdatadir)

    return


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
