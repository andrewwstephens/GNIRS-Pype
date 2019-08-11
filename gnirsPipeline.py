#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import ConfigParser
import log
import os
import gnirsGetData
import sort_data
import make_lists
import link_cals
import gnirsCheckData
import gnirsBaselineCalibration
import gnirsReduce
import gnirsCombineSpectra2D
import gnirsExtractSpectra1D
import gnirsTelluric
#import gnirsFluxCalibrate
#import gnirsCombineOrdersXD
#import gnirsCalculateSpectrumSNR
#import writeDataSheet

__version__ = "2019.08.11"


# ----------------------------------------------------------------------------------------------------------------------
def start(args):
    """
    This script is a full-featured GNIRS data reduction pipeline.

    This script:
        - gets data reduction parameters; either from an interactive input session (not adapted for GNIRS as of July 
          2019) or an input file, and
        - launches appropriate scripts to do the work:
                1) gnirsGetData.py
                2) gnirsSort.py
                3) gnirsCheckCalibrations.py
                4) gnirsBaselineCalibration.py
                5) gnirsReduce.py (science and/or telluric)
    """
    log.configure('gnirs.log', filelevel='INFO', screenlevel=args.loglevel)
    logger = log.getLogger('gnirsPipeline')

    path = os.getcwd()  # Save starting path for later use

    # Get paths to built-in Nifty data. Special code in setup.py makes sure recipes/ and runtimeData/ will be installed 
    # when someone installs Nifty, and accessible in this way.

    logger.info(" -------------------------------------- ")
    logger.info("|                                      |")
    logger.info("|              GNIRS-Pype              |")
    logger.info("| GNIRS Python Data Reduction Pipeline |")
    logger.info("|          Gemini Observatory          |")
    logger.info("|          Version %s          |", __version__)
    logger.info("|                                      |")
    logger.info(" -------------------------------------- ")

    # Make sure to change this if you change the default logfile.
    logger.info('The log file is gnirs.log.')

    logger.info("Reading parameters from %s\n", args.config)
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

    manualMode = config.getboolean('defaults','manualMode')

    ####################################################################################
    #                         SETUP COMPLETE                                           #
    #                      BEGIN DATA REDUCTION                                        #
    #                                                                                  #
    #   Twelve Steps:                                                                  #
    #       1) Get raw data - gnirsData.py                                             #
    #       2) Sort raw data - sort_data, make_lists, link_cals                        #
    #       3) Check raw data - gnirsCheckData.py                                      #
    #       4) Reduce baseline calibrations - gnirsBaselineCalibration.py              #
    #       5) Reduce observations (science and/or telluric) - gnirsReduce.py          #
    #       6) Combine 2D spectra - gnirsCombine2Dspectra.py                           #
    #       7) Extract 1D spectra - gnirsExtract1Dspectra.py                           #
    #       8) Perform telluric correction - gnirsTelluric.py                          #
    #       9) Perform flux calibration - gnirsFluxCalibrate.py                        #
    #      10) Combine XD Orders - gnirsCombineOrdersXD.py                             #
    #      11) Calculate the SNR spectrum - gnirsCalculateSpectrumSNR.py               #
    #      12) Write data sheet - gnirsWriteDataSheet.py                               #
    #                                                                                  #
    ####################################################################################

    #########################################################################
    #                      STEP 1: Get raw data                             #
    #########################################################################
    
    if config.getboolean('gnirsPipeline', 'getData'):
        if manualMode:
            a = raw_input('About to enter gnirsGetData to locate the raw data directory or download data from the Gemini public archive.')
        gnirsGetData.start(args.config)

    # ------------------------------------------------------------------------------------------------------------------
    # STEP 2: Sort the data into directories, make the file lists, and link the best calibrations:
    if config.getboolean('gnirsPipeline', 'sort'):
        sort_data.start(args.config)
        make_lists.start(args.config)
        link_cals.start(args.config)

    #########################################################################
    #                         STEP 3: Check raw Data                        #
    #########################################################################

    if config.getboolean('gnirsPipeline', 'checkData'):
        if manualMode:
            a = raw_input('About to enter gnirsCheckData.')
        gnirsCheckData.start(args.config)

    #########################################################################
    #                STEP 4: Reduce baseline calibrations                   #
    #########################################################################
    
    if config.getboolean('gnirsPipeline', 'calibrationReduction'):
        if manualMode:
            a = raw_input('About to enter gnirsBaselineCalibration to reduce baseline calibrations.')
        gnirsBaselineCalibration.start(args.config)

    #########################################################################
    #         STEP 5: Reduce observations (telluric and/or science)         #
    #########################################################################

    if config.getboolean('gnirsPipeline', 'scienceReduction'):
        if manualMode:
            a = raw_input('About to enter gnirsReduce to reduce science observations.')
        gnirsReduce.start('Science', args.config)

    if config.getboolean('gnirsPipeline', 'telluricReduction'):
        if manualMode:
            a = raw_input('About to enter gnirsReduce to reduce telluric observations.')
        gnirsReduce.start('Telluric', args.config)

    #########################################################################
    #                      STEP 6: Combine 2D spectra                       #
    #########################################################################

    if config.getboolean('gnirsPipeline', 'combineSpectra2D'):
        if manualMode:
            a = raw_input('About to enter gnirsCombineSpectra2D to combine different spectral orders into a single 2D spectrum.')
        gnirsCombineSpectra2D.start(args.config)

    #########################################################################
    #                      STEP 7: Extract 1D spectra                       #
    #########################################################################

    if config.getboolean('gnirsPipeline', 'extractSpectra1D'):
        if manualMode:
            a = raw_input('About to enter gnirsExtractSpectra1D to extract different spectral orders into a individual 1D spectra.')
        gnirsExtractSpectra1D.start(args.config)

    #########################################################################
    #                   STEP 8: Perform telluric correction                 #
    #########################################################################

    if config.getboolean('gnirsPipeline', 'telluricCorrection'):
        if manualMode:
            a = raw_input('About to enter gnirsTelluric to remove telluric absorption features from the extracted science spectra.')
        gnirsTelluric.start(args.config)

    #########################################################################
    #                    STEP 9: Perform flux calibration                   #
    #########################################################################

    if config.getboolean('gnirsPipeline', 'fluxCalibration'):
        if manualMode:
            a = raw_input('About to enter gnirsFluxCalibrate to create flux calibrated science spectra.')
        gnirsFluxCalibrate.start(args.config)

    #########################################################################
    #                       STEP 10: Combine XD Orders                      #
    #########################################################################

    if config.getboolean('gnirsPipeline', 'combineOrdersXD'):
        if manualMode:
            a = raw_input('About to enter gnirsCombineOrdersXD to combine diferent spectral orders.')
        combineOrdersXD.start(args.config)
    
    #########################################################################
    #                    STEP 11: Calculate SNR spectrum                    #
    #########################################################################

    if config.getboolean('gnirsPipeline', 'calculateSpectrumSNR'):
        if manualMode:
            a = raw_input('About to enter gnirsCalculateSpectrumSNR to calculate the SNR spectrum.')
        calculateSpectrumSNR.start(args.config)
    
    #########################################################################
    #                        STEP 12: Write Data Sheet                      #
    #########################################################################

    if config.getboolean('gnirsPipeline', 'writeDataSheet'):
        if manualMode:
            a = raw_input('About to enter gnirsWriteDataSheet to write a data sheet of results from GNIRS data reduction.')
        writeDataSheet.start(args.config)

    #########################################################################
    #                    Data Reduction Complete!                           #
    #                  Good luck with your science!                         #
    #########################################################################

    logger.info('###########################################################')
    logger.info('#                                                         #')
    logger.info('#               DATA REDUCTION COMPLETE                   #')
    logger.info('#             Good luck with your science!                #')
    logger.info('# Check out http://nifty4gemini.readthedocs.io/en/latest/ #')
    logger.info('#           For docs, tutorials and examples.             #')
    logger.info('#                                                         #')
    logger.info('###########################################################')

    return


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='This is the GNIRS Python Data Reduction Pipeline...')

    parser.add_argument('config', nargs='?', default='gnirs.cfg',
                        help='Configuration file')

    parser.add_argument('--loglevel', action='store', default='INFO',
                        choices=['DEBUG','INFO', 'WARNING', 'ERROR'],
                        help='Screen logger level')

    args = parser.parse_args()
    start(args)
