#!/usr/bin/env python

from xml.dom.minidom import parseString
from astropy.io import fits

import urllib, os, sys, ConfigParser, log
import numpy as np

# Import NDMapper gemini data download, by James E.H. Turner.
from downloadFromGeminiPublicArchive import download_query_gemini


def start(configfile):
    """
    Get the required raw GNIRS data.

    Args:
    -  configfile: gnirs.cfg configuration file.
    """
    
    log.configure('gnirs.log', filelevel='INFO', screenlevel='INFO')
    logger = log.getLogger('getData')

    # Store current working directory for later use.
    path = os.getcwd()

    # Set up the logger file.
    #log = os.getcwd()+'/gnirs.log'

    logger.info('####################################')
    logger.info('#                                  #')
    logger.info('#       Get raw GNIRS data         #')
    logger.info('#                                  #')
    logger.info('####################################\n')

    logger.info("Parameters read from %s\n", configfile)
    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)
    # Read general config.
    manualMode = bool(config.get('defaults','manualMode'))
    overwrite = bool(config.get('defaults','overwrite'))
    # Read getData specific config.
    rawPath = config.get('getDataConfig','rawPath')
    program = config.get('getDataConfig','program')
    
    # Check for invalid command line input. Cannot both copy from Gemini and sort local files. Exit if -q <path to raw 
    # frame files> and -c True are specified at command line (cannot copy from Gemini North internal network AND use 
    # local raw data).
    if rawPath and program:
        logger.info("\n##############################################################################################")
        logger.info("###############################################################################################")
        logger.info("                                                                                               ")
        logger.info("           Error in getData: both a local path and a Gemini program ID (to download            ")
        logger.info("                             from Gemini Public Archive) were provided.                        ")
        logger.info("                                                                                               ")
        logger.info("###############################################################################################")
        logger.info("##############################################################################################\n")      
        raise SystemError

    # Download data from gemini public archive to ./rawData/.
    if program:
        url = 'https://archive.gemini.edu/download/'+ program + '/notengineering/NotFail/present/canonical'
        if not os.path.exists('./rawData'):
            os.mkdir('./rawData')
        logger.info('\nDownloading data from Gemini public archive to ./rawData. This will take a few minutes.')
        logger.info('\nURL used for the download: \n' + str(url))
        if proprietaryCookie:
            download_query_gemini(url, './rawData', proprietaryCookie)
        else:
            download_query_gemini(url, './rawData')
        rawPath = os.getcwd()+'/rawData'

    return rawPath

if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
