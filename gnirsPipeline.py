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

# STDLIB

import logging, os, sys, shutil, pkg_resources, argparse, datetime, ConfigParser

# LOCAL

# Import major Nifty Steps.
import gnirsData as gnirsData
import gnirsHeaders as gnirsHeaders
import gnirsSort as gnirsSort
#import gnirsBaselineCalibration as gnirsBaselineCalibration
#import gnirsReduce as gnirsReduce
#import gnirsTelluric as gnirsTelluric
#import gnirsFluxCalibrate as gnirsFluxCalibrate
#import gnirsCombineOrdersXD as gnirsCombineOrdersXD
# Import gnirs utilities module.
from gnirsUtils import datefmt
# Import custom pipeline setup Class.
#from objectoriented.GetConfig import GetConfig  ## this is not adapted for GNIRS as of July 2019
# Conveniently import some utility functions so we don't have to type the full name.

#                                +
#
#
#
#              +
#         +         +         +
#
#                     +      +
#
#
#      +       +   + + + + +    + + + +  + + + + +   +    +
#     + +     +       +        +            +         + +
#    +   +   +       +        + +          +           +
#   +     + +       +        +            +           +
#  +       +   + + + + +    +            +           +
#
#
#                                      +
#                                   +     +
#                                       +
#                                      +
#

# Welcome to Nifty!

# The current version:
# TODO(nat): fix this to import the version from setup.py.
__version__ = "1.0.1"

# The time when Nifty was started is:
startTime = str(datetime.datetime.now())
print(startTime)

## IMPORTANT NOTE:  The command line options are not available for GNIRS as of July 2019.

def start():
    """
    This script is a full-featured GNIRS data reduction pipeline. It can call up to five "Steps".

    This script does two things. It:
        - gets data reduction parameters; either from an interactive input session (not adapted for GNIRS as of July 2019) or an input file, and
        - launches appropriate scripts to do the work. It can call up to five scripts directly:
                1) gnirsData.py
                2) gnirsHeaders.py
                3) gnirsSort.py
                4) gnirsBaselineCalibration.py
                5) gnirsReduce.py
    """
    # Save starting path for later use and change one directory up.
    path = os.getcwd()
    print("IT WORKED!")
    # Get paths to built-in Nifty data. Special code in setup.py makes sure recipes/ and runtimeData/ will be installed when someone installs Nifty, 
    # and accessible in this way.
#    RECIPES_PATH = pkg_resources.resource_filename('GitHub', 'recipes/')  ## this is not adapted for GNIRS as of July 2019
    RECIPES_PATH = 'recipes/'
#    RUNTIME_DATA_PATH = pkg_resources.resource_filename('GitHub', 'runtimeData/')  ## this is not adapted for GNIRS as of July 2019
    RUNTIME_DATA_PATH = 'runtimeData/'

    # Format logging options.
    FORMAT = '%(asctime)s %(message)s'
    DATEFMT = datefmt()

    # Set up the main logging file.
    logging.basicConfig(filename='gnirs.log',format=FORMAT,datefmt=DATEFMT,level=logging.DEBUG)
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    # This lets us logging.info(to stdout AND a logfile. Cool, huh?
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logging.info("\n####################################")
    logging.info("#                                  #")
    logging.info("#             NIFTY                #")
    logging.info("#  GNIRS Data Reduction Pipeline   #")
    logging.info("#         Version "+ __version__+ "            #")
    logging.info("#         July 25th, 2017          #")
    logging.info("#     Marie Lemoine-Busserolle     #")
    logging.info("# Gemini Observatory, Hilo, Hawaii #")
    logging.info("#                                  #")
    logging.info("####################################\n")

    # Make sure to change this if you change the default logfile.
    logging.info('The log file is gnirs.log.')

    # Read or write a configuration file, interactively, from defaults or from a provided file.
    # Second argument is the name of the current script. This could be used to get script-dependent configuration.
#    GetConfig(args, "gnirsPipeline")  ## this functionality is not adapted for GNIRS as of July 2019

    logging.info("\nParameters for this data reduction as read from " + str(RECIPES_PATH) + "defaultConfig.cfg:\n")
    config = ConfigParser.RawConfigParser()
    config.read(RECIPES_PATH + 'defaultConfig.cfg')

    for i in config.sections():
        for k in config.options(i):
            logging.info(k + " " + config.get(i,k))
    logging.info("")
    
    # Load pipeline configuration from /recipes/defaultConfig.cfg that is used by this script.
    manualMode = config.getboolean('defaults','manualMode')

    # Load pipeline specific configuration.
    data = bool(config.get('gnirsPipelineConfig','data'))
    headers = bool(config.get('gnirsPipelineConfig','headers'))
    sort = bool(config.get('gnirsPipelineConfig','sort'))
    calibrationReduction = bool(config.get('gnirsPipelineConfig','calibrationReduction'))
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
    
    if data:
        if manualMode:
            a = raw_input('About to enter gnirsData.')
        data_directory = gnirsData.start()

    ###########################################################################
    ##                      STEP 2: Read the file headers.                   ##
    ###########################################################################
    
    if headers:
        if manualMode:
            a = raw_input('About to enter gnirsHeaders.')
        headerinfo = gnirsHeaders.start(data_directory)
    
    ###########################################################################
    ##                      STEP 3: Sort the raw data.                       ##
    ###########################################################################
    
    if sort:
        if manualMode:
            a = raw_input('About to enter gnirsSort.')
        gnirsSort.start(data_directory, headerinfo)
    # By now, we should have paths to the three types of raw data. Print them out.
    #printDirectoryLists()

    ###########################################################################
    ##                STEP 4: Reduce baseline calibrations.                  ##
    ###########################################################################
    '''
    if calibrationReduction:
        if manualMode:
            a = raw_input('About to enter gnirsBaselineCalibration.')
        gnirsBaselineCalibration.start()

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

    logging.info('\n###########################################################')
    logging.info('#                                                         #')
    logging.info('#               DATA REDUCTION COMPLETE                   #')
    logging.info('#             Good luck with your science!                #')
    logging.info('# Check out http://nifty4gemini.readthedocs.io/en/latest/ #')
    logging.info('#           For docs, tutorials and examples.             #')
    logging.info('#                                                         #')
    logging.info('###########################################################\n')

    return
    '''
if __name__ == '__main__':
    # This block could let us call start gnirsPipeline.py from the command line. It is disabled for now.
    # start(args)
#    pass
    start()
