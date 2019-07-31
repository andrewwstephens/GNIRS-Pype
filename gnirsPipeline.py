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

import gnirsGetData
import gnirsSort
import gnirsCheckData
import gnirsBaselineCalibration
import gnirsReduce
import gnirsCombineSpectra2D
import gnirsExtractSpectra1D
import gnirsTelluric
import gnirsFluxCalibrate
import gnirsCombineOrdersXD
import gnirsCalculateSpectrumSNR
import writeDataSheet

# Import custom pipeline setup Class.
#from objectoriented.GetConfig import GetConfig  ## this is not adapted for GNIRS as of July 2019
# Conveniently import some utility functions so we don't have to type the full name.

# The current version:
# TODO(nat): fix this to import the version from setup.py.
__version__ = "1.0.1"

## IMPORTANT NOTE:  The command line options are not available for GNIRS as of July 2019.


def start(args):
    """
    This script is a full-featured GNIRS data reduction pipeline. It can call up to five "Steps".

    This script does two things. It:
        - gets data reduction parameters; either from an interactive input session (not adapted for GNIRS as of July 
          2019) or an input file, and
        - launches appropriate scripts to do the work. It can call up to five scripts directly:             ---------- TBD??
                1) gnirsGetData.py
                2) gnirsSort.py
                3) gnirsCheckCalibrations.py
                4) gnirsBaselineCalibration.py
                5) gnirsReduce.py (science and/or telluric)
    """
    log.configure('gnirs.log', filelevel='INFO', screenlevel=args.loglevel)
    logger = log.getLogger('gnirsPipeline')

    # Save starting path for later use and change one directory up.
    path = os.getcwd()

    # Get paths to built-in Nifty data. Special code in setup.py makes sure recipes/ and runtimeData/ will be installed 
    # when someone installs Nifty, and accessible in this way.

    logger.info("\n####################################")
    logger.info("#                                  #")
    logger.info("#             NIFTY                #")
    logger.info("#  GNIRS Data Reduction Pipeline   #")
    logger.info("#         Version %5s            #", __version__)
    logger.info("#         July 25th, 2017          #")
    logger.info("#     Marie Lemoine-Busserolle     #")
    logger.info("# Gemini Observatory, Hilo, Hawaii #")
    logger.info("#                                  #")
    logger.info("####################################\n")

    # Make sure to change this if you change the default logfile.
    logger.info('The log file is gnirs.log.')

    # Read or write a configuration file, interactively, from defaults or from a provided file.
    # Second argument is the name of the current script. This could be used to get script-dependent configuration.
#    GetConfig(args, "gnirsPipeline")  ## this functionality is not adapted for GNIRS as of July 2019

    logger.info("Parameters read from %s\n", args.config)
    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(args.config)

    for section in config.sections():
        logger.info('[%s]', section)
        for option in config.options(section):
            logger.info(option + " = " + config.get(section,option))
    logger.info("")

    # Andy wonders why there are TWO config files: gnirs.cfg and defaultConfig.cfg              
    # Viraja:  If the script gives the flexibility of how the user wishes to reduce the data, which I think gnirs.cfg \
    # does offer, then we do not need a gnirsDefault.cfg
    
    # Load general configuration from gnirs.cfg.
    manualMode = config.getboolean('defaults','manualMode')
    # Load pipeline specific configuration.
    getData = config.getboolean('gnirsPipeline','getData')
    sort = config.getboolean('gnirsPipeline','sort')
    checkData = config.getboolean('gnirsPipeline','checkData')
    calibrationReduction = config.getboolean('gnirsPipeline','calibrationReduction')
    telluricReduction = config.getboolean('gnirsPipeline','telluricReduction')
    scienceReduction = config.getboolean('gnirsPipeline','scienceReduction')
    combineSpectra2D = config.getboolean('gnirsPipeline','combineSpectra2D')
    extractSpectra1D = config.getboolean('gnirsPipeline','extractSpectra1D')
    telluricCorrection = config.getboolean('gnirsPipeline','telluricCorrection')
    fluxCalibration = config.getboolean('gnirsPipeline','fluxCalibration')
    combineOrdersXD = config.getboolean('gnirsPipeline','combineOrdersXD')
    calculateSpectrumSNR = config.getboolean('gnirsPipeline','calculateSpectrumSNR')
    writeDataSheet = config.getboolean('gnirsPipeline','writeDataSheet')

    ######################################################################################
    ##                         SETUP COMPLETE                                           ##
    ##                      BEGIN DATA REDUCTION                                        ##
    ##                                                                                  ##
    ##   Twelve Steps:                                                                  ##  
    ##       1) Get raw data - gnirsData.py                                             ##
    ##       2) Sort raw data - gnirsSort.py                                            ##
    ##       3) Check raw data - gnirsCheckData.py                                      ##
    ##       4) Reduce baseline calibrations - gnirsBaselineCalibration.py              ##
    ##       5) Reduce observations (science and/or telluric) - gnirsReduce.py          ##
    ##       6) Combine 2D spectra - gnirsCombine2Dspectra.py                           ##
    ##       7) Extract 1D spectra - gnirsExtract1Dspectra.py                           ##
    ##       8) Perform telluric correction - gnirsTelluric.py                          ##
    ##       9) Perform flux calibration - gnirsFluxCalibrate.py                        ##
    ##      10) Combine XD Orders - gnirsCombineOrdersXD.py                             ##
    ##      11) Calculate the SNR spectrum - gnirsCalculateSpectrumSNR.py               ##
    ##      12) Write data sheet - gnirsWriteDataSheet.py                               ##
    ##                                                                                  ##
    ######################################################################################

    ###########################################################################
    ##                      STEP 1: Get raw data                             ##
    ###########################################################################
    
    if getData:
        if manualMode:
            a = raw_input('About to enter gnirsGetData to locate the raw data directory or download data from the Gemini public archive.')
        gnirsGetData.start(args.config)
    
    ###########################################################################
    ##                        STEP 2: Sort raw data                          ##
    ###########################################################################

    if sort:
        if manualMode:
            a = raw_input('About to enter gnirsSort to sort observation types, and copy files and lists into appropriate directories for reduction.')
        gnirsSort.start(args.config)

    ###########################################################################
    ##                         STEP 3: Check raw Data                        ##
    ###########################################################################

    if checkData:
        if manualMode:
            a = raw_input('About to enter gnirsCheckData.')
        gnirsCheckData.start(args.config)

    ###########################################################################
    ##                STEP 4: Reduce baseline calibrations                   ##
    ###########################################################################
    
    if calibrationReduction:
        if manualMode:
            a = raw_input('About to enter gnirsBaselineCalibration to reduce baseline calibrations.')
        gnirsBaselineCalibration.start(args.config)

    ###########################################################################
    ##         STEP 5: Reduce observations (telluric and/or science)         ##
    ###########################################################################

    if scienceReduction:
        if manualMode:
            a = raw_input('About to enter gnirsReduce to reduce science observations.')
        gnirsReduce.start('Science', args.config)

    if telluricReduction:
        if manualMode:
            a = raw_input('About to enter gnirsReduce to reduce telluric observations.')
        gnirsReduce.start('Telluric', args.config)

    ###########################################################################
    ##                      STEP 6: Combine 2D spectra                       ##
    ###########################################################################

    if combineSpectra2D:
        if manualMode:
            a = raw_input('About to enter gnirsCombineSpectra2D to combine different spectral orders into a single 2D spectrum.')
        gnirsCombineSpectra2D.start(args.config)

    ###########################################################################
    ##                      STEP 7: Extract 1D spectra                       ##
    ###########################################################################

    if extractSpectra1D:
        if manualMode:
            a = raw_input('About to enter gnirsExtractSpectra1D to extract different spectral orders into a individual 1D spectra.')
        gnirsExtractSpectra1D.start(args.config)

    ###########################################################################
    ##                   STEP 8: Perform telluric correction                 ##
    ###########################################################################

    if telluricCorrection:
        if manualMode:
            a = raw_input('About to enter gnirsTelluric to remove telluric absorption features from the extracted science spectra.')
        gnirsTelluricCorrection.start(args.config)

    ###########################################################################
    ##                    STEP 9: Perform flux calibration                  ##
    ###########################################################################

    if fluxCalibration:
        if manualMode:
            a = raw_input('About to enter gnirsFluxCalibrate to create flux calibrated science spectra.')
        gnirsFluxCalibrate.start(args.config)

    ###########################################################################
    ##                       STEP 10: Combine XD Orders                      ##
    ###########################################################################

    if combineOrdersXD:
        if manualMode:
            a = raw_input('About to enter gnirsCombineOrdersXD to combine diferent spectral orders.')
        combineOrdersXD.start(args.config)
    
    ###########################################################################
    ##                    STEP 11: Calculate SNR spectrum                    ##
    ###########################################################################

    if calculateSpectrumSNR:
        if manualMode:
            a = raw_input('About to enter gnirsCalculateSpectrumSNR to calculate the SNR spectrum.')
        calculateSpectrumSNR.start(args.config)
    
    ###########################################################################
    ##                        STEP 12: Write Data Sheet                      ##
    ###########################################################################

    if writeDataSheet:
        if manualMode:
            a = raw_input('About to enter gnirsWriteDataSheet to write a data sheet of results from GNIRS data reduction.')
        writeDataSheet.start(args.config)

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
    logger.info('###########################################################\n')

    return

#---------------------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='This is the GNIRS Data Reduction Pipeline...')

    parser.add_argument('config', nargs='?', default='gnirs.cfg',
                        help='Configuration file')

    parser.add_argument('--loglevel', action='store', default='INFO',
                        choices=['DEBUG','INFO', 'WARNING', 'ERROR'],
                        help='Screen logger level')

    args = parser.parse_args()
    start(args)
