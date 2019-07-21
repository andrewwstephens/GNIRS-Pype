#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MIT License

# Copyright (c) 2015, 2017 Marie Lemoine-Busserolle

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

################################################################################
#                Import some useful Python utilities/modules                   #
################################################################################

import os, sys, shutil, pkg_resources, argparse, datetime, ConfigParser
import log

import gnirsData
import gnirsHeaders
import gnirsSort
import gnirsBaselineCalibration
#import gnirsReduce
#import gnirsTelluric
#import gnirsFluxCalibrate
#import gnirsCombineOrdersXD
# Import gnirs utilities module.
from gnirsUtils import datefmt

# TODO(nat): fix this to import the version from setup.py.
__version__ = "1.0.1"


def start(args):
    """
    This script is a full-featured GNIRS data reduction pipeline. It can call up to five "Steps".

    This script does two things:
        - gets data reduction parameters; either from an interactive input session (not adapted for GNIRS as of July 2019) or an input file, and
        - launches appropriate scripts to do the work. It can call up to five scripts directly:
                1) gnirsData.py
                2) gnirsHeaders.py
                3) gnirsSort.py
                4) gnirsBaselineCalibration.py
                5) gnirsReduce.py
    """

    log.configure('gnirs.log', filelevel='INFO', screenlevel=args.loglevel)
    logger = log.getLogger('Pipeline')

    # Save starting path for later use and change one directory up.
    path = os.getcwd()

#    RECIPES_PATH = pkg_resources.resource_filename('GitHub', 'recipes/')  ## this is not adapted for GNIRS as of July 2019
    RECIPES_PATH = 'recipes/'
#    RUNTIME_DATA_PATH = pkg_resources.resource_filename('GitHub', 'runtimeData/')  ## this is not adapted for GNIRS as of July 2019
    RUNTIME_DATA_PATH = 'runtimeData/'

    logger.info("####################################")
    logger.info("#                                  #")
    logger.info("#             NIFTY                #")
    logger.info("#  GNIRS Data Reduction Pipeline   #")
    logger.info("#         Version %5s            #", __version__)
    logger.info("#         July 25th, 2017          #")
    logger.info("#     Marie Lemoine-Busserolle     #")
    logger.info("# Gemini Observatory, Hilo, Hawaii #")
    logger.info("#                                  #")
    logger.info("####################################")

    # Make sure to change this if you change the default logfile.
    logger.info('The log file is gnirs.log.')

    # Read or write a configuration file, interactively, from defaults or from a provided file.
    # Second argument is the name of the current script. This could be used to get script-dependent configuration.
#    GetConfig(args, "gnirsPipeline")  ## this functionality is not adapted for GNIRS as of July 2019

    logger.info("Parameters read from %s", args.cfg)
    config = ConfigParser.RawConfigParser()
    config.read(args.cfg)

    for section in config.sections():
        logger.info('[%s]', section)
        for option in config.options(section):
            logger.info(option + " = " + config.get(section,option))
    logger.info("")

    # Andy wonders why there are TWO config files: config.cfg and defaultConfig.cfg

    # Load pipeline configuration from recipes/defaultConfig.cfg
    manualMode = config.getboolean('defaults','manualMode')
    calibrationReduction = config.getboolean('gnirsPipelineConfig','calibrationReduction')
    telluricReduction = bool(config.get('gnirsPipelineConfig','telluricReduction'))
    scienceReduction = bool(config.get('gnirsPipelineConfig','scienceReduction'))
    telluricCorrection = bool(config.get('gnirsPipelineConfig','telluricCorrection'))
    fluxCalibration = bool(config.get('gnirsPipelineConfig','fluxCalibration'))
    combineOrdersXD = bool(config.get('gnirsPipelineConfig','combineOrdersXD'))

    ###########################################################################
    ##                         SETUP COMPLETE                                ##
    ##                      BEGIN DATA REDUCTION                             ##
    ##                                                                       ##
    ##   Six Main Steps:                                                    ##  
    ##       1) Get Data - gnirsData.py                                      ##
    ##       2) Read Headers - gnirsHeaders.py                               ##
    ##       2) Sort the Raw Data - gnirsSort.py                             ##
    ##       3) Reduce baseline calibrations - gnirsBaselineCalibration.py   ##
    ##       4) Reduce telluric observations - gnirsReduce.py                ##
    ##       5) Reduce science observations - gnirsReduce.py                 ##
    ##                                                                       ##
    ###########################################################################

    ###########################################################################
    ##                      STEP 1: Get the data.                            ##
    ###########################################################################
    
    if config.getboolean('gnirsPipelineConfig','getdata'):
        if config.getboolean('defaults','manualMode'):
            a = raw_input('About to enter gnirsData.')
        data_directory = gnirsData.start(config)

    ###########################################################################
    ##                      STEP 2: Read the file headers.                   ##
    ###########################################################################
    
    if config.getboolean('gnirsPipelineConfig','headers'):
        if config.getboolean('defaults','manualMode'):
            a = raw_input('About to enter gnirsHeaders.')
        headerinfo = gnirsHeaders.start(data_directory)
    
    ###########################################################################
    #                       STEP 3: Sort the raw data.                        #
    ###########################################################################

    if config.getboolean('gnirsPipelineConfig','sort'):
        if config.getboolean('defaults','manualMode'):
            a = raw_input('About to enter gnirsSort.')
        gnirsSort.start(headerinfo, config)

    # By now, we should have paths to the three types of raw data. Print them out.
    #printDirectoryLists()

    ###########################################################################
    ##                STEP 4: Reduce baseline calibrations.                  ##
    ###########################################################################
    
    if calibrationReduction:
        if manualMode:
            a = raw_input('About to enter gnirsBaselineCalibration.')
        gnirsBaselineCalibration.start()


    """
    ###########################################################################
    ##                STEP 5: Reduce telluric observations.                  ##
    ###########################################################################

    if telluricReduction:
        if manualMode:
            a = raw_input('About to enter gnirsReduce to reduce Tellurics.')
        gnirsReduce.start('Telluric')

    ###########################################################################
    ##                 STEP 6: Reduce science observations.                  ##
    ###########################################################################

    if scienceReduction:
        if manualMode:
            a = raw_input('About to enter gnirsReduce to reduce science.')
        gnirsReduce.start('Science')
    if telluricCorrection:
        if manualMode:
            a = raw_input('About to enter gnirsTelluric to make and create telluric corrected spectra.')
        gnirsTelluric.run()

    if fluxCalibration:
        if manualMode:
            a = raw_input('About to enter gnirsFluxCalibrate to make and create flux calibrated and telluric corrected '
                          'spectra.')
        gnirsFluxCalibrate.run()

    if combineOrders:
        if manualMode:
            a = raw_input('About to enter gnirsMerge to combine spectra from different orders into a single spectrum.')
        combineOrders.run()

    ###########################################################################
    ##                    Data Reduction Complete!                           ##
    ##                  Good luck with your science!                         ##
    ###########################################################################

    logger.info('###########################################################')
    logger.info('#                                                         #')
    logger.info('#               DATA REDUCTION COMPLETE                   #')
    logger.info('#             Good luck with your science!                #')
    logger.info('# Check out http://nifty4gemini.readthedocs.io/en/latest/ #')
    logger.info('#           For docs, tutorials and examples.             #')
    logger.info('#                                                         #')
    logger.info('###########################################################')

    """
    return


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='This is the GNIRS Data Reduction Pipeline...')

    parser.add_argument('cfg', nargs='?', default='gnirs.cfg',
                        help='Configuration file')

    parser.add_argument('--loglevel', action='store', default='INFO',
                        choices=['DEBUG','INFO', 'WARNING', 'ERROR'],
                        help='Screen logger level')

    args = parser.parse_args()
    start(args)
