#!/usr/bin/env python

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

import log, os, sys, glob, ConfigParser, datetime, dateutil.parser
from astropy.io import fits
from pyraf import iraf, iraffunctions
import numpy as np

## IMPORTANT NOTE:  The command line options are not available for GNIRS as of July 2019.


def start(configfile):
    """
    This module contains all the functions needed to perform the full reduction of SCIENCE or TELLURIC data.

    Parameters are loaded from gnirs.cfg configuration file. This script will automatically detect if it is being run
    on telluric data or science data. There are 5 steps.

    INPUT FILES:
        - Configuration file
        - Science or Telluric frames
        - mdfshiftrefimage
        - masterflat
        - /database files from the appropriate calibrations directory

    OUTPUT FILES:
        - If telluric:  cleaned (optional), prepared, radiation-event corrected, reduced, spatial distortion corrected, 
          and transformed images
        - If science:  cleaned (optional), prepared, radiation-event corrected, reduced, spatial distortion corrected, 
          and transformed images

    Args:
        - kind (string): Either 'Science' or 'Telluric'
        - configfile: gnirs.cfg configuration file.
                - Paths to the Science (str), reduction truth value (boolean)
                  E.g. 'target/date/config/{Sci,Tel}_ObsID/{Calibrations,Intermediate}', True
                - Paths to the Tellurics (str), reduction truth value (boolean)
                  E.g. 'target/date/config/{Sci,Tel}_ObsID/{Calibrations,Intermediate}', True
                - manualMode (boolean): Enable optional manualModeging pauses? Default: False
                - overwrite (boolean): Overwrite old files? Default: False
                # And gnirsReduce specific settings
    """
    logger = log.getLogger('gnirsReduce.start')

    # TODO(nat): Right now the pipeline will crash if you decide to skip, say, doing a bad pixel correction. This is 
    # because each step adds a prefix to the frame name, and most following steps depend on that prefix being there.
    # One way to fix this is if a step is to be skipped, iraf.copy() is called instead to copy the frame and add the 
    # needed prefix. Messy but it might work for now.

    ###########################################################################
    ##                                                                       ##
    ##                  BEGIN - GENERAL EXTRACT 1D SETUP                     ##
    ##                                                                       ##
    ###########################################################################

    # Store current working directory for later use.
    path = os.getcwd()

    logger.info('####################################################')
    logger.info('#                                                  #')
    logger.info('#        Start Extracting GNIRS 1D Spectra         #')
    logger.info('#                                                  #')
    logger.info('####################################################')

    # Set up/prepare IRAF.
    iraf.gemini()
    iraf.gemtools()
    iraf.gnirs()

    # Reset to default parameters the used IRAF tasks.
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs,iraf.imcopy)

    # From http://bishop.astro.pomona.edu/Penprase/webdocuments/iraf/beg/beg-image.html:
    # Before doing anything involving image display the environment variable stdimage must be set to the correct frame 
    # buffer size for the display servers (as described in the dev$graphcap file under the section "STDIMAGE devices") 
    # or to the correct image display device. The task GDEVICES is helpful for determining this information for the 
    # display servers.
    iraf.set(stdimage='imt1024')

    # Prepare the IRAF package for GNIRS.
    # NSHEADERS lists the header parameters used by the various tasks in the GNIRS package (excluding headers values 
    # which have values fixed by IRAF or FITS conventions).
    iraf.nsheaders("gnirs",logfile=logger.root.handlers[0].baseFilename)

    # Set clobber to 'yes' for the script. This still does not make the gemini tasks overwrite files, so: YOU WILL 
    # LIKELY HAVE TO REMOVE FILES IF YOU RE_RUN THE SCRIPT.
    us_clobber=iraf.envget("clobber")
    iraf.reset(clobber='yes')
    
    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)
    # Read general config.
    manualMode = config.getboolean('defaults','manualMode')
    overwrite = config.getboolean('defaults','overwrite')

    # config required for extracting 1D spectra
    observationSections = ['ScienceDirectories','TelluricDirectories']
    nsextractInter = config.getboolean('interactive','nsextractInter')
    
    # extract1Spectra1D specific config
    useApall = config.getboolean('extract1Spectra1D','useApall')
    extractionApertureRadius = config.getfloat('extract1Spectra1D','extractionApertureRadius')
    extractStepwise = config.getboolean('extract1Spectra1D','extractStepwise')
    extractionStepsize = config.getfloat('extract1Spectra1D','extractionStepsize')

    ###########################################################################
    ##                                                                       ##
    ##                 COMPLETE - GENERAL EXTRACT 1D SETUP                   ##
    ##                                                                       ##
    ###########################################################################

    # gnirsExtractSpectra1D will first check if the reduction truth value of the science and telluric directories is 
    # True -- if it is, it will then check if the required spectra to be extracted are available in the directories
    # (and proceed only if it finds them there); else, it will warn the user and request to provide the spectra for 
    # extracting. If the reduction truth value of the science and telluric directories is False, the script will skip
    # extracting 1D spectra in those directories.

    # Loop through all the observation (telluric and science) directories to extract 1D spectra in each one.
    for section in observationSections:
        for obspath in config.options(section):
            if config.getboolean(section,obspath):

                ###########################################################################
                ##                                                                       ##
                ##                  BEGIN - OBSERVATION SPECIFIC SETUP                   ##
                ##                                                                       ##
                ###########################################################################

                os.chdir(obspath)
                # Change the iraf directory to the current directory.
                iraffunctions.chdir(obspath)

                # Print the current directory of observations being reduced.
                logger.info("Currently working on extracting 1D spectra in %s\n", obspath)
                
                allobsfilename = 'all.list'
                allobslist = open(allobsfilename, "r").readlines()                                 
                allobslist = [filename.strip() for filename in allobslist]
                rawHeader = fits.open(allobslist[0])[0].header
                srclistfilename = 'src.list'
                srclist = open(srclistfilename, "r").readlines()
                srclist = [filename.strip() for filename in srclist]
                skylistfilename = 'sky.list'
                skylist = open(skylistfilename, "r").readlines()
                skylist = [filename.strip() for filename in skylist]

                # Check if calibrations done in the calibrations directory path
                logger.info("Checking if required calibrations available in %s", calpath)
                
                QHflatslist = open(calpath+'/QHflats.list', "r").readlines()
                mdfshiftimage = calpath+'/n'+QHflatslist[0].strip()
                logger.debug("%s", mdfshiftimage) 
                if os.path.exists(mdfshiftimage):  ## Check for the reference image to calculate MDF
                    logger.info("Reference image to calculate MDF information available.")
                    calflag = True  ## calflag defined here
                else:
                    logger.warning("Reference image to calculate MDF information not available.")
                    calflag = False  ## calflag defined here
                
                masterflat = calpath+'/masterflat.fits'
                if os.path.exists(masterflat):  ## Check for the masterflat
                    logger.info("Masterflat to apply flat field correction available.")
                    calflag = calflag and True
                else:
                    logger.warning("Masterflat to apply flat field correction not available.")
                    calflag = calflag and False
                
                databasepath = calpath+'/database'  
                if os.path.exists(databasepath):  ## Check for /database directory (possibly) containing calibration files
                    logger.info("Reference /database directory (possibly) containing calibration files available.")
                    calflag = calflag and True
                else:
                    logger.warning("Reference /database directory (possibly) containing calibration files not available.")
                    calflag = calflag and False
                
                sdistfiles = glob.glob(databasepath+'/*_sdist')
                if not sdistfiles:  ## Check for spatial distortion correction calibration files
                    logger.info("Reference files to apply spatial distortion correction not available in the /database ")
                    logger.info("directory.")
                    calflag = calflag and False
                else:
                    logger.warning("Reference files to apply spatial distortion correction available in the /database ")
                    logger.warning("directory.")
                    calflag = calflag and True
                    sdistfileslength = len(sdistfiles)
                
                wavecallampfiles = glob.glob(databasepath+'/*_lamp')
                if not wavecallampfiles:  ## Check for spectral transformation calibration files
                    logger.info("Reference files to apply spectral transformation not available in the /database ")
                    logger.info("directory.")
                    calflag = calflag and False
                else:
                    logger.warning("Reference files to apply spectral transformation available in the /database ")
                    logger.warning("directory.")
                    calflag = calflag and True
                    wavecallampfileslength = len(wavecallampfiles)
                
                logger.info("Calibrations check complete.")
                if calflag:
                    logger.info("All reference calibration files available in %s", calpath)
                else:
                    logger.warning("One or more reference calibration files not available in %s ", calpath + ". Please ")
                    logger.warning("run gnirsBaselineCalibration.py to create the required calibration files or update ")
                    logger.warning("them manually. Exiting script.")
                    raise SystemExit

                combinedarc = 'arc_comb.fits'  ## The combined arc filename used in the calibrations directory path
                
                ###########################################################################
                ##                                                                       ##
                ##                 COMPLETE - OBSERVATION SPECIFIC SETUP                 ##
                ##            BEGIN EXTRACTING 1D SPECTRA FOR AN OBSERVATION             ##
                ##                                                                       ##
                ###########################################################################
                    
                if manualMode:
                    a = raw_input("About to enter step 7: extract 1D spectra.")

                if kind == 'Telluric':
                    allcomb = 'all_comb.fits'
                    if useApall:
                        apertureTracingColumns = 20
                        extract1Dspectra(allcomb, useApall, apertureTracingColumns, databasepath, \
                            extractionApertureRadius, overwrite)
                    else:
                        apertureTracingColumns = 10
                        extract1Dspectra(allcomb, useApall, apertureTracingColumns, databasepath, \
                            extractionApertureRadius, overwrite)
                if kind == 'Science':
                    # For science frames, only the source combined spectrum is extracted unless the user has 
                    # specified that the SNR spectrum is to be calculated in which case the sky combined spectrum 
                    # will also be extracted.
                    srccomb = 'src_comb.fits'
                    if useApall:
                        # This performs a weighted extraction
                        apertureTracingColumns = 20
                        extract1Dspectra(srccomb, useApall, apertureTracingColumns, databasepath, \
                            extractionApertureRadius, overwrite)
                        if calcSNRspectrum:
                            logger.info("Extracting the combined sky spectrum reduced without sky subtraction.\n")
                            skycomb = 'sky_comb.fits'
                            extract1Dspectra(skycomb, useApall, apertureTracingColumns, databasepath, \
                                extractionApertureRadius, overwrite)
                    else:
                        apertureTracingColumns = 10
                        extract1Dspectra(srccomb, useApall, apertureTracingColumns, databasepath, \
                            extractionApertureRadius, overwrite)
                        if calcSNRspectrum:
                            logger.info("Extracting the combined sky spectrum reduced without sky subtraction.\n")
                            skycomb = 'sky_comb.fits'
                            extract1Dspectra(skycomb, useApall, apertureTracingColumns, databasepath, \
                                extractionApertureRadius, overwrite)

            logger.info("##############################################################################")
            logger.info("#                                                                            #")
            logger.info("#  COMPLETE - Extracting 1D spectra completed for                            #")
            logger.info("#  %s", obspath                                                                )
            logger.info("#                                                                            #")
            logger.info("##############################################################################\n")

    # Return to directory script was begun from.
    os.chdir(path)

    return

##################################################################################################################
#                                                     ROUTINES                                                   #
##################################################################################################################

def extract1Dspectra(combinedimage, useApall, apertureTracingColumns, databasepath, extractionApertureRadius, overwrite):
    """
    Extracting 1D spectra from the combined 2D spectra using nsextract.
    """
    logger = log.getLogger('gnirsReduce.extractSpectra1D')

    if os.path.exists('v'+combinedimage):
        if overwrite:
            logger.warning("Removing old v%s", combinedimage)
            os.remove('v'+combinedimage)
            iraf.nsextract(inimages=combinedimage, outspectra='', outprefix='v', dispaxis=1, database=databasepath, \
                line=700, nsum=apertureTracingColumns, ylevel='INDEF', upper=str(extractionApertureRadius), \
                lower='-'+str(extractionApertureRadius), background='none', fl_vardq='yes', fl_addvar='no', \
                fl_skylines='yes', fl_inter='no', fl_apall=useApall, fl_trace='no', \
                aptable='gnirs$data/apertures.fits', fl_usetabap='no', fl_flipped='yes', fl_project='yes', \
                fl_findneg='no', bgsample='*', trace='', tr_nsum=10, tr_step=10, tr_nlost=3, tr_function='legendre', \
                tr_order=5, tr_sample='*', tr_naver=1, tr_niter=0, tr_lowrej=3.0, tr_highrej=3.0, tr_grow=0.0, \
                weights='variance', logfile=logger.root.handlers[0].baseFilename, verbose='yes', mode='al')
        else:
            logger.warning("Old %s exists and -overwrite not set - skipping nsextract for observations.", combinedimage)
    else:
        iraf.nsextract(inimages=combinedimage, outspectra='', outprefix='v', dispaxis=1, database=databasepath, \
            line=700, nsum=apertureTracingColumns, ylevel='INDEF', upper=str(extractionApertureRadius), \
            lower='-'+str(extractionApertureRadius), background='none', fl_vardq='yes', fl_addvar='no', \
            fl_skylines='yes', fl_inter='no', fl_apall=useApall, fl_trace='no', aptable='gnirs$data/apertures.fits', \
            fl_usetabap='no', fl_flipped='yes', fl_project='yes', fl_findneg='no', bgsample='*', trace='', tr_nsum=10,\
            tr_step=10, tr_nlost=3, tr_function='legendre', tr_order=5, tr_sample='*', tr_naver=1, tr_niter=0, \
            tr_lowrej=3.0, tr_highrej=3.0, tr_grow=0.0, weights='variance', \
            logfile=logger.root.handlers[0].baseFilename, verbose='yes', mode='al')

#---------------------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
