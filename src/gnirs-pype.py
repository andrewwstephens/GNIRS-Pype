#!/usr/bin/env python

import argparse
import ConfigParser
import log
import download_data
import sort_data
import make_lists
import link_cals
import check_data
import baseline_cals
import reduce
import gnirsCombineSpectra2D
import extract_spectra
import telluric_correct
import gnirsGetTelluricInfo
import gnirsFluxCalibrate
import combine_orders
import noise_spectrum
import pdf_summary

__version__ = "2019.09.29"


# ----------------------------------------------------------------------------------------------------------------------
def start(args):
    """
    GNIRS data reduction pipeline.
    """
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel=args.loglevel)
    logger = log.getLogger('gnirs-pype')

    logger.info(" -------------------------------------- ")
    logger.info("|              GNIRS-Pype              |")
    logger.info("| GNIRS Python Data Reduction Pipeline |")
    logger.info("|          Gemini Observatory          |")
    logger.info("|          Version %s          |", __version__)
    logger.info(" -------------------------------------- ")

    logger.info('The log file is gnirs-pype.log.')

    logger.info("Reading parameters from %s\n", args.config)
    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(args.config)
    for section in config.sections():
        logger.info('[%s]', section)
        for option in config.options(section):
            logger.info(option + " = " + config.get(section,option))
    logger.info('')

    # STEP 1: Download the raw data:
    if config.getboolean('gnirsPipeline', 'DownloadData'):
        download_data.start(args.config)

    # STEP 2: Sort the data into directories, make the file lists, and link the best calibrations:
    if config.getboolean('gnirsPipeline', 'sort'):
        sort_data.start(args.config)
        make_lists.start(args.config)
        link_cals.start(args.config)

    # STEP 3: Check raw Data
    if config.getboolean('gnirsPipeline', 'checkData'):
        check_data.start(args.config)

    # STEP 4: Reduce baseline calibrations
    if config.getboolean('gnirsPipeline', 'calibrationReduction'):
        baseline_cals.start(args.config)

    # STEP 5: Reduce observations (telluric and/or science)
    if config.getboolean('gnirsPipeline', 'scienceReduction'):
        reduce.start('Science', args.config)
    if config.getboolean('gnirsPipeline', 'telluricReduction'):
        reduce.start('Telluric', args.config)

    # STEP 6: Combine 2D spectra
    if config.getboolean('gnirsPipeline', 'combineSpectra2D'):
        gnirsCombineSpectra2D.start(args.config)

    # STEP 7: Extract spectra
    if config.getboolean('gnirsPipeline', 'extractSpectra'):
        extract_spectra.start(args.config)

    # STEP 8: Perform telluric correction
    if config.getboolean('gnirsPipeline', 'telluricCorrection'):
        telluric_correct.start(args.config)

    # STEP 9: Perform flux calibration
    if config.getboolean('gnirsPipeline', 'fluxCalibration'):
        gnirsGetTelluricInfo.start(args.config)
        gnirsFluxCalibrate.start(args.config)

    # STEP 10: Combine XD Orders
    if config.getboolean('gnirsPipeline', 'CombineOrders'):
        combine_orders.start(args.config)
    
    # STEP 11: Calculate noise and S/N spectrum
    if config.getboolean('gnirsPipeline', 'CalculateSNR'):
        noise_spectrum.start(args.config)

    # STEP 12: Write a PDF summary
    if config.getboolean('gnirsPipeline', 'PDFSummary'):
        pdf_summary.start(args.config)

    logger.info(' -------------------------------------------------')
    logger.info('|            DATA REDUCTION COMPLETE              |')
    logger.info('|          Good luck with your science!           |')
    logger.info('| See http://gnirs-pype.readthedocs.io/en/latest/ |')
    logger.info('|       for docs, tutorials and examples.         |')
    logger.info(' -------------------------------------------------')

    return


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='This is the GNIRS Python Data Reduction Pipeline.')

    parser.add_argument('config', nargs='?', default='gnirs-pype.cfg',
                        help='Configuration file')

    parser.add_argument('--loglevel', action='store', default='INFO',
                        choices=['DEBUG','INFO', 'WARNING', 'ERROR'],
                        help='Screen logger level')

    args = parser.parse_args()
    start(args)
