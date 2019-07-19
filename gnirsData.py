#!/usr/bin/env python

from xml.dom.minidom import parseString
import urllib, logging, os, glob, sys, shutil, math, pkg_resources, argparse, ConfigParser, time, datetime, re, datetime, dateutil.parser
from astropy.io import fits
import numpy as np
from gnirsUtils import datefmt

# Import NDMapper gemini data download, by James E.H. Turner.
from downloadFromGeminiPublicArchive import download_query_gemini

#RECIPES_PATH = pkg_resources.resource_filename('nifty', 'recipes/')
RECIPES_PATH = 'recipes/'
#RUNTIME_DATA_PATH = pkg_resources.resource_filename('nifty', 'runtimeData/')
RUNTIME_DATA_PATH = 'runtimeData/'

def start():
    """
    Get the required GNIRS data.

    Args:
        program (string): OT observation id (used only within Gemini network). Specified with -p at command line. "GN-2013B-Q-109".
    """
    '''
    # Store current working directory for later use.
    path = os.getcwd()
    '''
    # Format logging options.
    FORMAT = '%(asctime)s %(message)s'
    DATEFMT = datefmt()

    # Set up the logging file.
    log = os.getcwd()+'/gnirs.log'

    logging.info('\n####################################')
    logging.info('#                                  #')
    logging.info('#   Get the required GNIRS data    #')
    logging.info('#                                  #')
    logging.info('####################################\n')

    # Load reduction parameters from runtimeData/defaultConfig.cfg.
    config = ConfigParser.RawConfigParser()
    config.read(RECIPES_PATH + 'defaultConfig.cfg')
    # Read general config.
    overwrite = config.getboolean('defaults','overwrite')
    manualMode = config.getboolean('defaults','manualMode')
    # Read data specific config.
    rawPath = config.get('dataConfig','rawPath')
    program = config.get('dataConfig','program')
    
    # Check for invalid command line input. Cannot both copy from Gemini and sort local files. Exit if -q <path to raw frame files> and -c True are
    # specified at command line (cannot copy from Gemini North internal network AND use local raw data).
    if rawPath and program:
        logging.info("\nError in sort: both a local path and a Gemini program ID (to download from Gemini Public Archive) were provided.\n")
        raise SystemError

    # Download data from gemini public archive to ./rawData/.
    if program:
        url = 'https://archive.gemini.edu/download/'+ program + '/notengineering/NotFail/present/canonical'
        if not os.path.exists('./rawData'):
            os.mkdir('./rawData')
        logging.info('\nDownloading data from Gemini public archive to ./rawData. This will take a few minutes.')
        logging.info('\nURL used for the download: \n' + str(url))
        if proprietaryCookie:
            download_query_gemini(url, './rawData', proprietaryCookie)
        else:
            download_query_gemini(url, './rawData')
        rawPath = os.getcwd()+'/rawData'

    return rawPath

if __name__ == '__main__':
    # Set up logging if called as a standalone script.
    # Format logging options.
    FORMAT = '%(asctime)s %(message)s'
    DATEFMT = datefmt()

    # Set up the logging file.
    logging.basicConfig(filename='gnirsSort.log',format=FORMAT,datefmt=DATEFMT,level=logging.DEBUG)
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    # This lets us logging.info(to stdout AND a logfile. Cool, huh?
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Start gnirsSort from the beginning!
    start()
