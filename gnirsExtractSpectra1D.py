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

import log, os, sys, pkg_resources, glob, ConfigParser, cleanir, datetime, dateutil.parser
from astropy.io import fits
from pyraf import iraf, iraffunctions
import numpy as np

## IMPORTANT NOTE:  The command line options are not available for GNIRS as of July 2019.


def start(kind, configfile):
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
    ##                  BEGIN - GENERAL REDUCTION SETUP                      ##
    ##                                                                       ##
    ###########################################################################

    # Store current working directory for later use.
    path = os.getcwd()

    logger.info('####################################################')
    logger.info('#                                                  #')
    logger.info('#  Start the GNIRS Science and Telluric Reduction  #')
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
    
    if kind == 'Science':  # science reduction specific config
        observationSection = 'ScienceDirectories'
        observationDirectories = config.options('ScienceDirectories')
        start = config.getint('scienceReduction','Start')
        stop = config.getint('scienceReduction','Stop')
        cleanir = config.getboolean('scienceReduction','cleanir')
        radiationCorrectionMethod = config.get('scienceReduction','radiationCorrectionMethod')
        radiationThreshold = config.getfloat('scienceReduction','radiationThreshold')
        nscombineInter = config.getboolean('interactive','nscombineInter')
        calculateSpectrumSNR = config.getboolean('gnirsPipeline','calculateSpectrumSNR')
    elif kind == 'Telluric':  # telluric reduction specific config
        observationSection = 'TelluricDirectories'
        observationDirectories = config.options('TelluricDirectories')
        start = config.getint('telluricReduction','Start')
        stop = config.getint('telluricReduction','Stop')
        cleanir = config.getboolean('telluricReduction','cleanir')
        radiationCorrectionMethod = config.get('telluricReduction','radiationCorrectionMethod')
        radiationThreshold = config.getfloat('telluricReduction','radiationThreshold')
        nscombineInter = config.getboolean('interactive','nscombineInter')
        calculateSpectrumSNR = config.getboolean('gnirsPipeline','calculateSpectrumSNR')
    else:
        logger.error("###########################################################################")
        logger.error("###########################################################################")
        logger.error("#                                                                         #")
        logger.error("#     ERROR in reduce: invalid kind of reduction. Please enter either     #")
        logger.error("#                      <Science> or <Telluric>. Exiting script.           #")
        logger.error("#                                                                         #")
        logger.error("###########################################################################")
        logger.error("###########################################################################")
        raise SystemExit  
        ## TODO(Viraja): can ask for a raw input on the command line to avoid a system exit

    ###########################################################################
    ##                                                                       ##
    ##                 COMPLETE - GENERAL REDUCTION SETUP                    ##
    ##                                                                       ##
    ###########################################################################

    # gnirsReduce will reduce observations in each science or telluric directory only if the reduction truth value for
    # that directory is True; else, it will skip the reductions (of science or telluric frames) in that dorectory.

    # Loop through all the observation (telluric or science) directories to perform a reduction on each one.
    for obspath in observationDirectories:
        if config.getboolean(observationSection,obspath):

            ###########################################################################
            ##                                                                       ##
            ##                  BEGIN - OBSERVATION SPECIFIC SETUP                   ##
            ##                                                                       ##
            ###########################################################################

            os.chdir(obspath)
            # Change the iraf directory to the current directory.
            iraffunctions.chdir(obspath)

            # Print the current directory of observations being reduced.
            logger.info("Currently working on reductions in %s\n", obspath)

            tempObspath = obspath.split(os.sep)            
            calpath = "/".join(tempObspath[:-1])+'/Calibrations'
            # Print the directory from which reduced calibratios are being used.
            logger.info("The path to the calibrations is %s\n", calpath)
            
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
            ##                BEGIN DATA REDUCTION FOR AN OBSERVATION                ##
            ##                                                                       ##
            ###########################################################################
            
            # Check start and stop values for reduction steps. Ask user for a correction if input is not valid.
            valindex = start
            while valindex > stop  or valindex < 1 or stop > 7:
                logger.warning("#####################################################################")
                logger.warning("#####################################################################")
                logger.warning("#                                                                   #")
                logger.warning("#     WARNING in reduce: invalid start/stop values of observation   #")
                logger.warning("#                        reduction steps.                           #")
                logger.warning("#                                                                   #")
                logger.warning("#####################################################################")
                logger.warning("#####################################################################\n")

                valindex = int(raw_input("Please enter a valid start value (1 to 7, default 1): "))
                stop = int(raw_input("Please enter a valid stop value (1 to 7, default 7): "))

            while valindex <= stop :

                #############################################################################
                ##  STEP 1: Clean raw observations.                                        ##
                ##  Output: Cleaned science or telluric frames.                            ##
                #############################################################################

                if valindex == 1:
                    if manualMode:
                        a = raw_input("About to enter step 1: clean raw observations.")
                    
                    if cleanir:
 #                       cleanir(allobslist)
                        pass  ## this step is not modifoed to work with this pipeline as of July 2019  
                    else:
                        logger.info("######################################################################")
                        logger.info("######################################################################")
                        logger.info("#                                                                    #")
                        logger.info("#       WARNING in reduce: raw observations not cleaned.             #")
                        logger.info("#                                                                    #")
                        logger.info("######################################################################")
                        logger.info("######################################################################\n")
                    
                    logger.info("##################################################################")
                    logger.info("#                                                                #")
                    logger.info("#       STEP 1: Clean raw observations - COMPLETED               #")
                    logger.info("#                                                                #")
                    logger.info("##################################################################\n")
                    
                ###########################################################################
                ##    STEP 2: Prepare observations (science or telluric frames)          ##
                ###########################################################################

                if valindex == 2:
                    if manualMode:
                        a = raw_input("About to enter step 2: locate the spectrum.")
                    
                    # Use the most appropriate bad pixel mask. For data taken Before the summer 2012 lens replacement,
                    # use 'gnirs$data/gnirsn_2011apr07_bpm.fits'; after summer 2012, use 
                    # 'gnirs$data/gnirsn_2012dec05_bpm.fits'. Use keyword 'ARRAYID' in the raw file header to check
                    # which camera was used for observations, and accordingly set parameter 'bpm' for the bpmfile to be
                    # used from the gnirs$database.
                    arrayid = rawHeader['ARRAYID'].strip()
                    if arrayid == 'SN7638228.1':
                        bpmfile = 'gnirs$data/gnirsn_2011apr07_bpm.fits'
                    elif arrayid == 'SN7638228.1.2':
                        bpmfile = 'gnirs$data/gnirsn_2012dec05_bpm.fits'
                    else:
                        logger.error("######################################################################")
                        logger.error("######################################################################")
                        logger.error("#                                                                    #")
                        logger.error("#        ERROR in reduce: unknown array ID. Exiting script.          #")
                        logger.error("#                                                                    #")
                        logger.error("######################################################################")
                        logger.error("######################################################################\n")
                        raise SystemExit
                    
                    prepareObservations(mdfshiftimage, bpmfile, overwrite)

                    logger.info("##############################################################################")
                    logger.info("#                                                                            #")
                    logger.info("#       STEP 2: Locate the Spectrum (prepare observations) - COMPLETED       #")
                    logger.info("#                                                                            #")
                    logger.info("##############################################################################\n")
                
                ###########################################################################
                ##  STEP 3: Prepare bad pixel masks and do radiation-event correction    ##
                ###########################################################################

                elif valindex == 3:
                    if manualMode:
                        a = raw_input("About to enter step 3: radiation event correction.")

                    preparedHeader = fits.open('n'+allobslist[0])[0].header
                
                    date = dateutil.parser.parse(preparedHeader['DATE-OBS'].strip() + ' ' + preparedHeader['TIME-OBS'])
                    rdnoise = preparedHeader['RDNOISE']
                    gain = preparedHeader['GAIN']

                    if radiationEventCorrectionMethod == 'fixpix':
                        radiationEventCorrectionFixpix(allobslist, srclist, skylist, radiationThreshold, rdnoise, \
                            gain, overwrite)
                    elif radiationEventCorrectionMethod == 'dqplane':
                        radiationEventCorrectionDQplane(allobslist, srclist, skylist, radiationThreshold, rdnoise, \
                            gain, overwrite)
                    else:
                        logger.warning("########################################################################")
                        logger.warning("########################################################################")
                        logger.warning("#                                                                      #")
                        logger.warning("#   WARNING in reduce: invalid/no radiation event correction method    #")
                        logger.warning("#                      method found. Checking the observation date.    #")
                        logger.warning("#                                                                      #")
                        logger.warning("########################################################################")
                        logger.warning("########################################################################")
                        
                        if date < datetime.datetime(year=2012, month=7, day=1, hour=0, minute=0, second=0):
                            logger.info("Observation date is " + str(date) + ". Performing the radiation event ")
                            logger.info("correction using the 'dqplane' method.\n")
                            radiationEventCorrectionMethod = 'dqplane'
                            radiationEventCorrectionDQplane(allobslist, srclist, skylist, radiationThreshold, rdnoise,\
                                gain, overwrite)
                        else:
                            logger.info("Observation date is " + str(date) + ". GNIRS data from 2012B not affected ")
                            logger.info("by radiation events. Skipping radiation event correction.\n")
                            pass

                    logger.info("##############################################################################")
                    logger.info("#                                                                            #")
                    logger.info("#  STEP 3: Radiation event correction - COMPLETED                            #")
                    logger.info("#                                                                            #")
                    logger.info("##############################################################################\n")
                
                ##############################################################################
                ##  STEP 4: Flat field and subtract sky background (science or telluric)    ##
                ##############################################################################

                elif valindex == 4:
                    if manualMode:
                        a = raw_input("About to enter step 4: flat fielding and sky subtraction.")
                    
                    skySubtraction = 'yes'
                    reduce_outputPrefix = 'r'  ## output prefix assigned to the reduced observations
                    reduceObservations(skySubtraction, reduce_outputPrefix, masterflat, radiationThreshold, overwrite)
                    if kind == 'Science':
                        # If the parameter 'calcSNRspectrum' is set to 'yes' in the configuration file, the script will 
                        # reduce the science frames without sky subtraction (by setting the parameter 'fl_sky' in 
                        # nsreduce to 'no') in which case the parameter 'outprefix' in nsreduce will be 'k'; else, 
                        # the science frames will be reduced with sky subtraction and the output prefix will be 'r'.
                        if calcSNRspectrum:
                            logger.info("Reducing science frames without sky subtraction to generate an error spectrum ")   
                            logger.info("later by the pipeline.\n")
                            skySubtraction = 'no'
                            reduce_outputPrefix = 'k'
                            reduceObservations(skySubtraction, reduce_outputPrefix, masterflat, radiationThreshold, \
                                overwrite)
                
                    logger.info("##############################################################################")
                    logger.info("#                                                                            #")
                    logger.info("#  STEP 4: Flat fielding and sky subtraction - COMPLETED                     #")
                    logger.info("#                                                                            #")
                    logger.info("##############################################################################\n")
                
                #################################################################################
                ##  STEP 5: Apply spatial distortion correction and spectral transformation    ##
                #################################################################################

                elif valindex == 5:
                    if manualMode:
                        a = raw_input("About to enter step 5: spatial distortion correction and spectral \
                            transformation.")

                    reduce_outputPrefix = 'r'
                    SdistCorrection_SpectralTransform(databasepath, reduce_outputPrefix, allobslist, sdistfileslength,\
                        wavecallampfileslength, combinedarc, overwrite)
                    if kind == 'Science':
                        if calcSNRspectrum:
                            logger.info("Applying spatial distortion correction and spectral transformation to the "+ \
                                "science frames reduced without sky subtraction.\n")
                            reduce_outputPrefix = 'k'
                            SdistCorrection_SpectralTransform(databasepath, reduce_outputPrefix, allobslist, \
                                sdistfileslength, wavecallampfileslength, combinedarc, overwrite)

                    logger.info("##################################################################################")
                    logger.info("#                                                                                #")
                    logger.info("# STEP 5: spatial distortion correction and spectral transformation - COMPLETED  #")
                    logger.info("#                                                                                #")
                    logger.info("##################################################################################\n")
                
                #################################################################################
                ##  STEP 6: Combine 2D spectra (individual orders)                             ##                                       -------------------- which reduced sky frames should be combined for calculating the SNR spectrum later?
                #################################################################################

                elif valindex == 6:
                    if manualMode:
                        a = raw_input("About to enter step 6: combine 2D spectra.")

                    if kind == 'Telluric':
                        allcombinlist = 'all.list'
                        reduce_outputPrefix = 'r'
                        allcomb = 'all_comb.fits'
                        crossCorrelation = 'yes'
                        combine2Dspectra(allcombinlist, reduce_outputPrefix, allcomb, crossCorrelation, overwrite)
                    if kind == 'Science':
                        # For science frames, only the source observations are combined unless the user has specified 
                        # that the SNR spectrum is to be calculated in which case the sky observations will also be
                        # be combined.
                        srccombinlist = 'src.list'
                        reduce_outputPrefix = 'r'
                        srccomb = 'src_comb.fits'
                        crossCorrelation = 'no'
                        combine2Dspectra(srccombinlist, reduce_outputPrefix, srccomb, crossCorrelation, overwrite)
                        if calcSNRspectrum:
                            logger.info("Combining the sky observations reduced without sky subtraction.\n")
                            skycombinlist = 'sky.list'
                            reduce_outputPrefix = 'k'
                            skycomb = 'sky_comb.fits'
                            combine2Dspectra(skycombinlist, reduce_outputPrefix, skycomb, crossCorrelation, overwrite)

                    logger.info("###########################################################################")
                    logger.info("#                                                                         #")
                    logger.info("# STEP 6: combine 2D spectra - COMPLETED                                  #")
                    logger.info("#                                                                         #")
                    logger.info("###########################################################################\n")
                
                #################################################################################
                ##  STEP 7: Extract 1D spectra (individual orders)                             ##
                #################################################################################

                elif valindex == 7:
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

                    logger.info("#############################################################################")
                    logger.info("#                                                                           #")
                    logger.info("# STEP 7: extract 1D spectra - COMPLETED                                    #")
                    logger.info("#                                                                           #")
                    logger.info("#############################################################################\n")
                
                else:
                    logger.error("########################################################################")
                    logger.error("########################################################################")
                    logger.error("#                                                                      #")
                    logger.error("#  ERROR in reduce: %d", valindex + " is not valid. Exiting script.     #")
                    logger.error("#                                                                      #")
                    logger.error("########################################################################")
                    logger.error("########################################################################\n")
                    raise SystemExit
               
                valindex += 1

            logger.info("##############################################################################")
            logger.info("#                                                                            #")
            logger.info("#  COMPLETE - Reductions completed for                                       #")
            logger.info("#  %s", obspath                                                                )
            logger.info("#                                                                            #")
            logger.info("##############################################################################\n")

    # Return to directory script was begun from.
    os.chdir(path)

    return

##################################################################################################################
#                                                     ROUTINES                                                   #
##################################################################################################################

def prepareObservations(mdfshiftimage, bpmfile, overwrite):
    """
    Prepare raw observations (science or telluric) using nsprepare. Output prefix "n" is added to raw filenames.

    Processing with NSPREPARE (this task is used only for GNIRS data but other instruments have their own preparation 
    tasks with similar actions) will rename the data extension and add variance and data quality extensions. By default 
    (see NSHEADERS) the extension names are SCI for science data, VAR for variance, and DQ for data quality (0 = good). 
    NSPREPARE will add an MDF file (extension MDF) describing the GNIRS image pattern and how all images 
    cross-coorelate with the first pinhole flat frame supplied as the "mdfshiftimage" reference with the "shiftimage"
    parameter in NSPREPARE.
    """
    logger = log.getLogger('gnirsReduce.prepareObservations')

    # Update frames with mdf offset value and generate variance and data quality extensions. This code structure checks 
    # if iraf output files already exist. If output files exist and overwrite is specified, iraf output is overwritten.
    oldfiles = glob.glob('./nN*.fits')
    if not oldfiles:
        iraf.nsprepare(inimages='@all.list', rawpath='', outimages='', outprefix='n', bpm=bpmfile, \
            logfile=logger.root.handlers[0].baseFilename, fl_vardq='yes', fl_cravg='no', crradius=0.0, \
            fl_dark_mdf='no', fl_correct='no', fl_saturated='yes', fl_nonlinear='yes', fl_checkwcs='yes', \
            fl_forcewcs='yes', arraytable='gnirs$data/array.fits', configtable='gnirs$data/config.fits', \
            specsec='[*,*]', offsetsec='none', pixscale='0.15', shiftimage=mdfshiftimage, shiftx=0., shifty=0., \
            obstype='FLAT', fl_inter='no', verbose='yes', mode='al')
    else:
        if overwrite:
            logger.warning('Removing old nN*.fits files.')
            [os.remove(filename) for filename in oldfiles]
            iraf.nsprepare(inimages='@all.list', rawpath='', outimages='', outprefix='n', bpm=bpmfile, \
                logfile=logger.root.handlers[0].baseFilename, fl_vardq='yes', fl_cravg='no', crradius=0.0, \
                fl_dark_mdf='no', fl_correct='no', fl_saturated='yes', fl_nonlinear='yes', fl_checkwcs='yes', \
                fl_forcewcs='yes', arraytable='gnirs$data/array.fits', configtable='gnirs$data/config.fits', \
                specsec='[*,*]', offsetsec='none', pixscale='0.15', shiftimage=mdfshiftimage, shiftx=0., shifty=0., \
                obstype='FLAT', fl_inter='no', verbose='yes', mode='al')
        else:
            logger.warning("Output files exist and -overwrite not set - skipping nsprepare for all observations.")

#---------------------------------------------------------------------------------------------------------------------#

def radiationEventCorrectionFixpix(allobslist, srclist, skylist, radiationThreshold, rdnoise, gain, overwrite):
    """
    Prepare bad pixel masks for all observations (science or telluric) using the bad pixel interpolation method.

    Here, the task proto.fixpix is set to interpolate along lines to avoid bad effects when a radiation event is right\
    on a skyline. Column interpolation may be better for events landing on the peak of target spectra, but effects \
    from bad interpolation over skylines seem to do worse.
    """
    logger = log.getLogger('gnirsReduce.radiationEventCorrectionFixpix')

    minsrc = 'min_src.fits'
    minsky = 'min_sky.fits'
    oldfiles_minimages = glob.glob('./min*.fits')
    oldfiles_masks = glob.glob('./mask*')
    oldfiles_gmasks = glob.glob('./gmask*')
    oldfiles_radiationclean = glob.glob('./lnN*.fits')
    oldfiles_combinedmasks = glob.glob('./bgmask*')
    oldfiles = oldfiles_minimages + oldfiles_masks + oldfiles_gmasks + oldfiles_radiationclean + oldfiles_combinedmasks
    if not oldfiles:
        # Combine source and sky observations separately using GEMCOMBINE to create a "minimum" image using "minmax"
        # combining algorithm.
        logger.info("Creating a relatively clean 'minimum' image of the source observations.\n")  ##                                      -------------------- probably all telluric frames should be used to create a minimum image as none of them are on-sky observations
        iraf.gemcombine(input='n//@src.list', output=minsrc, title="", combine="median", reject="minmax", \
            offsets="none", masktype="goodvalue", maskvalue=0., scale="none", zero="none", weight="none", \
            statsec="[*,*]", expname="EXPTIME", lthreshold='INDEF', hthreshold='INDEF', nlow=0, nhigh=len(srclist)-1, \
            nkeep=1, mclip='yes', lsigma=3., hsigma=3., key_ron="RDNOISE", key_gain="GAIN", ron=0.0, gain=1.0, \
            snoise="0.0", sigscale=0.1, pclip=-0.5, grow=0.0, bpmfile="", nrejfile="", sci_ext="SCI", var_ext="VAR", \
            dq_ext="DQ", fl_vardq='yes', logfile=logger.root.handlers[0].baseFilename, fl_dqprop='no', verbose='yes')
        logger.info("Creating initial masks for the source observations.\n")
        for src in srclist:
            iraf.imexpr("(a-b)>"+str(radiationThreshold)+"*sqrt("+str(rdnoise)+"**2+2*b/"+str(gain)+") ? 1 : 0", \
                output='mask'+src[:-5]+'.pl', a='n'+src+'[SCI]', b=minsrc+'[SCI]', dims="auto", intype="int", \
                outtype="auto", refim="auto", bwidth=0, btype="nearest", bpixval=0., rangecheck='yes', verbose='yes', \
                exprdb="none")
        
        logger.info("Creating a relatively clean 'minimum' image of the sky observations.\n")
        iraf.gemcombine(input='n//@sky.list', output=minsky, title="", combine="median", reject="minmax", \
            offsets="none", masktype="goodvalue", maskvalue=0., scale="none", zero="none", weight="none", \
            statsec="[*,*]", expname="EXPTIME", lthreshold='INDEF', hthreshold='INDEF', nlow=0, nhigh=len(skylist)-1, \
            nkeep=1, mclip='yes', lsigma=3., hsigma=3., key_ron="RDNOISE", key_gain="GAIN", ron=0.0, gain=1.0, \
            snoise="0.0", sigscale=0.1, pclip=-0.5, grow=0.0, bpmfile="", nrejfile="", sci_ext="SCI", var_ext="VAR", \
            dq_ext="DQ", fl_vardq='yes', logfile=logger.root.handlers[0].baseFilename, fl_dqprop='no', verbose='yes')
        logger.info("Creating initial masks for the sky observations.\n")
        for sky in skylist:
            iraf.imexpr("(a-b)>"+str(radiationThreshold)+"*sqrt("+str(rdnoise)+"**2+2*b/"+str(gain)+") ? 1 : 0", \
                output='mask'+sky[:-5]+'.pl', a='n'+sky+'[SCI]', b=minsky+'[SCI]', dims="auto", intype="int", \
                outtype="auto", refim="auto", bwidth=0, btype="nearest", bpixval=0., rangecheck='yes', verbose='yes', \
                exprdb="none")
        
        with open('./masks.list', 'a+') as f:
            for filename in sorted(glob.glob('./mask*.pl')):
                filename = filename[filename.rfind('/')+1:]
                if filename not in f.read().split('\n'):
                    f.write(filename + '\n')
                    
        iraf.crgrow(input='@masks.list', output='g//@masks.list', radius=1.5, inval="INDEF", outval="INDEF")

        logger.info("Adding bad pixels from the [DQ] plane generated by nsprepare to the masks so that they can be \
            fixed as well. Interpolating over the radiation events.")
        for src in srclist:  ##                                                                                                           -------------------- consequently, use the same minimum image to interpolate over the bad pixels in all telluric frames
            iraf.copy(input='n'+src, output='ln'+src, verbose='no')
            iraf.imexpr(expr="a||b", output='bgmask'+src[:-5]+'.pl', a='gmask'+src[:-5]+'.pl', b=minsrc+'[DQ]', \
                dims="auto", intype="int", outtype="auto", refim="auto", bwidth=0, btype="nearest", bpixval=0., \
                rangecheck='yes', verbose='yes', exprdb="none")
            iraf.proto.fixpix(images='ln'+src+'[SCI,1]', masks='bgmask'+src[:-5]+'.pl', linterp=1, cinterp="INDEF", \
                verbose='no', pixels='no')

        for sky in skylist:
            iraf.copy(input='n'+sky, output='ln'+sky, verbose='no') 
            iraf.imexpr(expr="a||b", output='bgmask'+sky[:-5]+'.pl', a='gmask'+sky[:-5]+'.pl', b=minsky+'[DQ]', \
                dims="auto", intype="int", outtype="auto", refim="auto", bwidth=0, btype="nearest", bpixval=0., \
                rangecheck='yes', verbose='yes', exprdb="none")
            iraf.proto.fixpix(images='ln'+src+'[SCI,1]', masks='bgmask'+sky[:-5]+'.pl', linterp=1, cinterp="INDEF", \
                verbose='yes', pixels='no')
    else:
        if overwrite:
            logger.warning('Removing old %s and %s', minsrc, minsky)
            logger.warning('Removing old mask*, gmask*.pl, lnN*.fits, bgmask*.pl files.')
            [os.remove(filename) for filename in oldfiles]
            logger.info("Creating a relatively clean 'minimum' image of the source observations.\n")
            iraf.gemcombine(input='n//@src.list', output=minsrc, title="", combine="median", reject="minmax", \
                offsets="none", masktype="goodvalue", maskvalue=0., scale="none", zero="none", weight="none", \
                statsec="[*,*]", expname="EXPTIME", lthreshold='INDEF', hthreshold='INDEF', nlow=0, \
                nhigh=len(srclist)-1, nkeep=1, mclip='yes', lsigma=3., hsigma=3., key_ron= "RDNOISE", key_gain="GAIN",\
                ron=0.0, gain=1.0, snoise="0.0", sigscale=0.1, pclip=-0.5, grow=0.0, bpmfile="", nrejfile="", \
                sci_ext="SCI", var_ext="VAR", dq_ext="DQ", fl_vardq='yes', \
                logfile=logger.root.handlers[0].baseFilename, fl_dqprop='no', verbose='yes')
            logger.info("Creating initial masks for the source observations.\n")
            for src in srclist:
                iraf.imexpr(expr="(a-b)>"+str(radiationThreshold)+"*sqrt("+str(rdnoise)+"**2+2*b/"+str(gain)+") ? 1 : 0",\
                    output='mask'+src[:-5]+'.pl', a='n'+src+'[SCI]', b=minsrc+'[SCI]', dims="auto", intype="int", \
                    outtype="auto", refim="auto", bwidth=0, btype="nearest", bpixval=0., rangecheck='yes', \
                    verbose='yes', exprdb="none")
            
            logger.info("Creating a relatively clean 'minimum' image of the sky observations.\n")
            iraf.gemcombine(input='n//@sky.list', output=minsky, title="", combine="median", reject="minmax", \
                offsets="none", masktype="goodvalue", maskvalue=0., scale="none", zero="none", weight="none", \
                statsec="[*,*]", expname="EXPTIME", lthreshold='INDEF', hthreshold='INDEF', nlow=0, \
                nhigh=len(skylist)-1, nkeep=1, mclip='yes', lsigma=3., hsigma=3., key_ron="RDNOISE", key_gain="GAIN",\
                ron=0.0, gain=1.0, snoise="0.0", sigscale=0.1, pclip=-0.5, grow=0.0, bpmfile="", nrejfile="", \
                sci_ext="SCI", var_ext="VAR", dq_ext="DQ", fl_vardq='yes', logfile=logger.root.handlers[0].baseFilename,\
                fl_dqprop='no', verbose='yes')
            logger.info("Creating initial masks for the sky observations.\n")
            for sky in skylist:
                iraf.imexpr(expr="(a-b)>"+str(radiationThreshold)+"*sqrt("+str(rdnoise)+"**2+2*b/"+str(gain)+") ? 1 : 0",\
                    output='mask'+sky[:-5]+'.pl', a='n'+sky+'[SCI]', b=minsrc+'[SCI]', dims="auto", intype="int", \
                    outtype="auto", refim="auto", bwidth=0, btype="nearest", bpixval=0., rangecheck='yes', \
                    verbose='yes', exprdb="none")

            with open('./masks.list', 'a+') as f:
                for filename in sorted(glob.glob('./mask*.pl')):
                    filename = filename[filename.rfind('/')+1:]
                    if filename not in f.read().split('\n'):
                        f.write(filename + '\n')

            iraf.crgrow(input='@masks.list', output='g//@masks.list', radius=1.5, inval="INDEF", outval="INDEF")

            logger.info("Adding bad pixels from the [DQ] plane generated by nsprepare to the masks so that they can be\
                fixed as well. Interpolating over the radiation events.")
            for src in srclist:
                iraf.copy(input='n'+src, output='ln'+src, verbose='no') 
                iraf.imexpr(expr="a||b", output='bgmask'+src[:-5]+'.pl', a='gmask'+src[:-5]+'.pl', b=minsrc+'[DQ]', \
                    dims="auto", intype="int", outtype="auto", refim="auto", bwidth=0, btype="nearest", bpixval=0., \
                    rangecheck='yes', verbose='yes', exprdb="none")
                iraf.proto.fixpix(images='ln'+src+'[SCI,1]', masks='bgmask'+src[:-5]+'.pl', linterp=1, \
                    cinterp="INDEF", verbose='yes', pixels='no')

            for sky in skylist:
                iraf.copy(input='n'+sky, output='ln'+sky, verbose='no') 
                iraf.imexpr(expr="a||b", output='bgmask'+sky[:-5]+'.pl', a='gmask'+sky[:-5]+'.pl', b=minsky+'[DQ]', \
                    dims="auto", intype="int", outtype="auto", refim="auto", bwidth=0, btype="nearest", bpixval=0., \
                    rangecheck='yes', verbose='yes', exprdb="none")
                iraf.proto.fixpix(images='ln'+src+'[SCI,1]', masks='bgmask'+sky[:-5]+'.pl', linterp=1, \
                    cinterp="INDEF", verbose='yes', pixels='no')
        else:
            logger.warning("Output files exist and -overwrite not set - skipping radiation event correction for all \
                observations.")
            '''
#---------------------------------------------------------------------------------------------------------------------#

def radiationEventCorrectionDQplane(allobslist, srclist, skylist, radiationThreshold, rdnoise, gain, overwrite):
    """
    Prepare bad pixel masks for all observations (science or telluric) using the data quality plane of the prepared 
    observations.

    TODO(Viraja):  Describe 'dqplane' method.
    """
    logger = log.getLogger('gnirsReduce.radiationEventCorrectionDQplane')

    minsrc = 'min_src.fits'
    minsky = 'min_sky.fits'
    oldfiles_minimages = glob.glob('./min*.fits')
    oldfiles_masks = glob.glob('./mask*')
    oldfiles_gmasks = glob.glob('./gmask*')
    oldfiles_radiationclean = glob.glob('./lnN*.fits')
    oldfiles_combinedmasks = glob.glob('./bgmask*')
    oldfiles = oldfiles_minimages + oldfiles_masks + oldfiles_gmasks + oldfiles_radiationclean + oldfiles_combinedmasks
    if not oldfiles:
        # Combine source and sky observations separately using GEMCOMBINE to create a "minimum" image using "minmax"
        # combining algorithm.
        logger.info("Creating a relatively clean 'minimum' image of the source observations.\n")
        iraf.gemcombine(input='n//@src.list', output=minsrc, title="", combine="median", reject="minmax", \
            offsets="none", masktype="goodvalue", maskvalue=0., scale="none", zero="none", weight="none", \
            statsec="[*,*]", expname="EXPTIME", lthreshold='INDEF', hthreshold='INDEF', nlow=0, nhigh=len(srclist)-1, \
            nkeep=1, mclip='yes', lsigma=3., hsigma=3., key_ron="RDNOISE", key_gain="GAIN", ron=0.0, gain=1.0, \
            snoise="0.0", sigscale=0.1, pclip=-0.5, grow=0.0, bpmfile="", nrejfile="", sci_ext="SCI", var_ext="VAR", \
            dq_ext="DQ", fl_vardq='yes', logfile=logger.root.handlers[0].baseFilename, fl_dqprop='no', verbose='yes')
        logger.info("Creating a relatively clean 'minimum' image of the sky observations.\n")
        iraf.gemcombine(input='n//@sky.list', output=minsky, title="", combine="median", reject="minmax", \
            offsets="none", masktype="goodvalue", maskvalue=0., scale="none", zero="none", weight="none", \
            statsec="[*,*]", expname="EXPTIME", lthreshold='INDEF', hthreshold='INDEF', nlow=0, nhigh=len(skylist)-1, \
            nkeep=1, mclip='yes', lsigma=3., hsigma=3., key_ron="RDNOISE", key_gain="GAIN", ron=0.0, gain=1.0, \
            snoise="0.0", sigscale=0.1, pclip=-0.5, grow=0.0, bpmfile="", nrejfile="", sci_ext="SCI", var_ext="VAR", \
            dq_ext="DQ", fl_vardq='yes', logfile=logger.root.handlers[0].baseFilename, fl_dqprop='no', verbose='yes')
        
        logger.info("Copying the [DQ] planes generated by nsprepare as the initial masks for all observations.\n")  ##                                      -------------------- in this case, we do not need to combine the images to create a minimum image in the previous step, right?
        for obs in allobslist:
            iraf.imcopy(input='n'+obs+'[DQ,1]', output='mask'+obs[:-5]+'.pl', verbose='yes')
        
        with open('./masks.list', 'a+') as f:
            for filename in sorted(glob.glob('./mask*.pl')):
                filename = filename[filename.rfind('/')+1:]
                if filename not in f.read().split('\n'):
                    f.write(filename + '\n')
    
        iraf.crgrow(input='@masks.list', output='g//@masks.list', radius=1.5, inval="INDEF", outval="INDEF")

        logger.info("Adding 'minimum' images of the source and the sky observations to the masks so that they can be \
            replaced as well. Replacing only the pixels identified as bad pixels by the final masks with corresponding\
            pixels from the 'minimum' images.")
        for src in srclist:  ##                                                                                                                 -------------------- or should we still have a minimum image and add the masks from the [DQ] planes to it? 
            iraf.copy(input='n'+src, output='ln'+src, verbose='no') 
            iraf.imexpr(expr="a||b", output='bgmask'+src[:-5]+'.pl', a='gmask'+src[:-5]+'.pl', b=minsrc+'[DQ]', \
                dims="auto", intype="int", outtype="auto", refim="auto", bwidth=0, btype="nearest", bpixval=0., \
                rangecheck='yes', verbose='yes', exprdb="none")
            iraf.imexpr(expr="c>1 ? b : a", output='ln'+src+'[SCI,overwrite]', a='n'+src+'[SCI]', b=minsrc+'[SCI]', \
                c='bgmask'+src[:-5]+'.pl', dims="auto", intype="int", outtype="auto", refim="auto", bwidth=0, \
                btype="nearest", bpixval=0., rangecheck='yes', verbose='yes', exprdb="none")

        for sky in skylist:
            iraf.copy(input='n'+sky, output='ln'+sky, verbose='no') 
            iraf.imexpr(expr="a||b", output='bgmask'+sky[:-5]+'.pl', a='gmask'+sky[:-5]+'.pl', b=minsky+'[DQ]', \
                dims="auto", intype="int", outtype="auto", refim="auto", bwidth=0, btype="nearest", bpixval=0., \
                rangecheck='yes', verbose='yes', exprdb="none")
            iraf.imexpr(expr="c>1 ? b : a", output='ln'+sky+'[SCI,overwrite]', a='n'+sky+'[SCI]', b=minsky+'[SCI]', \
                c='bgmask'+sky[:-5]+'.pl', dims="auto", intype="int", outtype="auto", refim="auto", bwidth=0, \
                btype="nearest", bpixval=0., rangecheck='yes', verbose='yes', exprdb="none")
    else:
        if overwrite:
            logger.warning('Removing old %s and %s', minsrc, minsky)
            logger.warning('Removing old mask*, gmask*.pl, lnN*.fits, bgmask*.pl files.')
            [os.remove(filename) for filename in oldfiles]
            logger.info("Creating a relatively clean 'minimum' image of the source observations.\n")
            iraf.gemcombine(input='n//@src.list', output=minsrc, title="", combine="median", reject="minmax", \
                offsets="none", masktype="goodvalue", maskvalue=0., scale="none", zero="none", weight="none", \
                statsec="[*,*]", expname="EXPTIME", lthreshold='INDEF', hthreshold='INDEF', nlow=0, \
                nhigh=len(srclist)-1, nkeep=1, mclip='yes', lsigma=3., hsigma=3., key_ron="RDNOISE", key_gain="GAIN", \
                ron=0.0, gain=1.0, snoise="0.0", sigscale=0.1, pclip=-0.5, grow=0.0, bpmfile="", nrejfile="", \
                sci_ext="SCI", var_ext="VAR", dq_ext="DQ", fl_vardq='yes', \
                logfile=logger.root.handlers[0].baseFilename, fl_dqprop='no', verbose='yes')
            logger.info("Creating a relatively clean 'minimum' image of the sky observations.\n")
            iraf.gemcombine(input='n//@sky.list', output=minsky, title="", combine="median", reject="minmax", \
                offsets="none", masktype="goodvalue", maskvalue=0., scale="none", zero="none", weight="none", \
                statsec="[*,*]", expname="EXPTIME", lthreshold='INDEF', hthreshold='INDEF', nlow=0, \
                nhigh=len(skylist)-1, nkeep=1, mclip='yes', lsigma=3., hsigma=3., key_ron="RDNOISE", key_gain="GAIN", \
                ron=0.0, gain=1.0, snoise="0.0", sigscale=0.1, pclip=-0.5, grow=0.0, bpmfile="", nrejfile="", \
                sci_ext="SCI", var_ext="VAR", dq_ext="DQ", fl_vardq='yes', \
                logfile=logger.root.handlers[0].baseFilename, fl_dqprop='no', verbose='yes')
            
            logger.info("Copying the [DQ] planes generated by nsprepare as the initial masks for all observations.\n")
            for obs in allobslist:
                iraf.imcopy(input='n'+obs+'[DQ,1]', output='mask'+obs[:-5]+'.pl', verbose='yes')
            
            with open('./masks.list', 'a+') as f:
                for filename in sorted(glob.glob('./mask*.pl')):
                    filename = filename[filename.rfind('/')+1:]
                    if filename not in f.read().split('\n'):
                        f.write(filename + '\n')

            iraf.crgrow(input='@masks.list', output='g//@masks.list', radius=1.5, inval="INDEF", outval="INDEF")
            
            logger.info("Adding 'minimum' images of the source and the sky observations to the masks so that they can be \
                replaced as well. Replacing only the pixels identified as bad pixels by the final masks with corresponding\
                pixels from the 'minimum' images.")
            for src in srclist:
                iraf.copy(input='n'+src, output='ln'+src, verbose='no') 
                iraf.imexpr(expr="a||b", output='bgmask'+src[:-5]+'.pl', a='gmask'+src[:-5]+'.pl', b=minsrc+'[DQ]', \
                    dims="auto", intype="int", outtype="auto", refim="auto", bwidth=0, btype="nearest", bpixval=0., \
                    rangecheck='yes', verbose='yes', exprdb="none")
                iraf.imexpr(expr="c>1 ? b : a", output='ln'+src+'[SCI,overwrite]', a='n'+src+'[SCI]', \
                    b=minsrc+'[SCI]', c='bgmask'+src[:-5]+'.pl', dims="auto", intype="int", outtype="auto", \
                    refim="auto", bwidth=0, btype="nearest", bpixval=0., rangecheck='yes', verbose='yes', \
                    exprdb="none")

            for sky in skylist:
                iraf.copy(input='n'+sky, output='ln'+sky, verbose='no') 
                iraf.imexpr(expr="a||b", output='bgmask'+sky[:-5]+'.pl', a='gmask'+sky[:-5]+'.pl', b=minsky+'[DQ]', \
                    dims="auto", intype="int", outtype="auto", refim="auto", bwidth=0, btype="nearest", bpixval=0., \
                    rangecheck='yes', verbose='yes', exprdb="none")
                iraf.imexpr(expr="c>1 ? b : a", output='ln'+sky+'[SCI,overwrite]', a='n'+sky+'[SCI]', \
                    b=minsky+'[SCI]', c='bgmask'+sky[:-5]+'.pl', dims="auto", intype="int", outtype="auto", \
                    refim="auto", bwidth=0, btype="nearest", bpixval=0., rangecheck='yes', verbose='yes', \
                    exprdb="none")
        else:
            logger.warning("Output files exist and -overwrite not set - skipping radiation event correction for all \
                observations.")

#---------------------------------------------------------------------------------------------------------------------#

def reduceObservations(skySubtraction, reduce_outputPrefix, masterflat, radiationThreshold, overwrite):
    """
    Flat field and sky subtract observations with nsreduce. If 'calcSNRspectrum' is set, observations will be reduced
    without sky subtraction too and prefix 'k' will be added before the input filenames; otherwise, sky subtraction
    will be performed and the prefix will be 'r'.

    NSREDUCE is used for basic reduction of raw data - it provides a single, unified interface to several tasks and 
    also allows the subtraction of dark frames (not done for GNIRS data) and dividing by the flat.
    """
    logger = log.getLogger('gnirsReduce.reduceObservations')

    oldfiles = glob.glob('./'+reduce_outputPrefix+'lnN*.fits')
    if not oldfiles:
        iraf.nsreduce(inimages='ln//@all.list', outimages='', outprefix=reduce_outputPrefix, fl_cut='yes', section="",\
            fl_corner='yes', fl_process_cut='yes', fl_nsappwave='no', nsappwavedb='gnirs$data/nsappwave.fits', \
            crval="INDEF", cdelt="INDEF", fl_dark='no', darkimage="", fl_save_dark='no', fl_sky=skySubtraction, \
            skyimages="", skysection="", combtype="median", rejtype="avsigclip", masktype="goodvalue", maskvalue=0., 
            scale="none", zero="median", weight="none", statsec="[*,*]", lthreshold="INDEF", hthreshold="INDEF", \
            nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3., hsigma=3., snoise="0.0", sigscale=0.1, pclip=-0.5, \
            grow=0.0, skyrange='INDEF', nodsize=3., fl_flat='yes', flatimage=masterflat, flatmin=0.0, fl_vardq='yes', \
            logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', force='no')
        # Record radiation threshold value in reduced file headers
        iraf.hedit(images=reduce_outputPrefix+'ln//@all.list', fields='RTHRESH', value=radiationThreshold, add='yes', \
            addonly='no', delete='no', verify='no', show='no', update='yes')
    else:
        if overwrite:
            logger.warning('Removing old %slnN*.fits files.', reduce_outputPrefix)
            [os.remove(filename) for filename in oldfiles]
            iraf.nsreduce(inimages='ln//@all.list', outimages='', outprefix=reduce_outputPrefix, fl_cut='yes', \
                section="", fl_corner='yes', fl_process_cut='yes', fl_nsappwave='no', \
                nsappwavedb='gnirs$data/nsappwave.fits', crval="INDEF", cdelt="INDEF", fl_dark='no', darkimage="", \
                fl_save_dark='no', fl_sky=skySubtraction, skyimages="", skysection="", combtype="median", \
                rejtype="avsigclip", masktype="goodvalue", maskvalue=0., scale="none", zero="median", weight="none", \
                statsec="[*,*]", lthreshold="INDEF", hthreshold="INDEF", nlow=1, nhigh=1, nkeep=0, mclip='yes', \
                lsigma=3., hsigma=3., snoise="0.0", sigscale=0.1, pclip=-0.5, grow=0.0, skyrange='INDEF', nodsize=3., \
                fl_flat='yes', flatimage=masterflat, flatmin=0.0, fl_vardq='yes', \
                logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', force='no')
            iraf.hedit(images=reduce_outputPrefix+'ln//@all.list', fields='RTHRESH', value=radiationThreshold, \
                add='yes', addonly='no', delete='no', verify='no', show='no', update='yes')
            
        else:
            logger.warning("Output files exist and -overwrite not set - skipping nsreduce for all observations.")

#---------------------------------------------------------------------------------------------------------------------#

def SdistCorrection_SpectralTransform(databasepath, reduce_outputPrefix, allobslist, sdistfileslength, wavecallampfileslength, combinedarc, overwrite):
    """
    Apply spatial distortion correction and spectral transformation. 

    In general, we would run nsfitcoords on all reduced science or telluric frames with and without sky subtraction (if 
    creating an error spectrum) to straighten the spectra and then de-tilt. nsfitcoords simply writes database info, 
    which depends only on pinholes and arcs (the *_sdist and *_lamp files that were written to the databse while 
    reducing baseline calibrations), not on on-sky data. So, here we access the nsfitcoords solution from the 
    calibrations folder and append the relevant header keywords, which will tell nstransform the database files it 
    should use, without actually calling nsfitcoords. This routine also uses the output files in the database that were 
    created by the nswavelength task.
    """
    logger = log.getLogger('gnirsReduce.SdistCorrection_SpectralTransform')

    oldfiles_fitcoords = glob.glob('./f'+reduce_outputPrefix+'lnN*')
    oldfiles_sdisttrans = glob.glob('./tf'+reduce_outputPrefix+'lnN*')
    oldfiles_spectrans = glob.glob('./ttf'+reduce_outputPrefix+'lnN*')
    oldfiles = oldfiles_fitcoords + oldfiles_sdisttrans + oldfiles_spectrans
    if not oldfiles:
        logger.info("Editing the primary header AND science extensions.")
        for obs in allobslist:
            # First, make copies of the reduced observations with new filenames and edit header to add the 
            # nsfitcoords stamp.
            iraf.copy(input=reduce_outputPrefix+'ln'+obs, output='f'+reduce_outputPrefix+'ln'+obs, verbose='yes', \
                mode="ql")
            iraf.hedit(images='f'+reduce_outputPrefix+'ln'+obs+'[0]', fields="NSFITCOO", value="copied from arc", \
                add='yes', addonly='no', delete='no', verify='no', show='no', update='yes', mode='al')
            
            # Second, add sdist info and transform the files to make the orders vertical.
            for i in range(sdistfileslength):
                iraf.hedit(images='f'+reduce_outputPrefix+'ln'+obs+'[SCI,'+str(i+1)+']', field="FCFIT2", \
                    value='f'+combinedarc[:-5]+'_SCI_'+str(i+1)+'_sdist', add='yes', addonly='no', delete='no', \
                    verify='no', show='no', update='yes', mode='al')
            iraf.nstransform(inimages='f'+reduce_outputPrefix+'ln'+obs, outspectra='', outprefix='t', dispaxis=1, \
                database=databasepath, fl_stripe='no', interptype='poly3', xlog='no', ylog='no', pixscale=1.0, \
                logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', mode='al')
            
            # Third, add lamp info and transform the files to make skylines horizontal (spatial direction along 
            # detector rows).
            for i in range(wavecallampfileslength):
                iraf.hedit(images='tf'+reduce_outputPrefix+'ln'+obs+'[SCI,'+str(i+1)+']', field="FCFIT1", \
                    value='ftf'+combinedarc[:-5]+'_SCI_'+str(i+1)+'_lamp', add='yes', addonly='no', delete='no', \
                    verify='no', show='no', update='yes', mode='al')
            iraf.nstransform(inimages='tf'+reduce_outputPrefix+'ln'+obs, outspectra='', outprefix='t', dispaxis=1, \
                database=databasepath, fl_stripe='no', interptype='poly3', xlog='no', ylog='no', pixscale=1.0, \
                logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', mode='al')
    else:
        if overwrite:
            logger.warning("Removing old f%slnN*.fits, tf%slnN*.fits, ttf%slnN*.fits files.", reduce_outputPrefix, \
                reduce_outputPrefix, reduce_outputPrefix)
            [os.remove(filename) for filename in oldfiles]
            logger.info("Editing the primary header AND science extensions.")
            for obs in allobslist:
                # First, make copies of the reduced observations with new filenames and edit header to add the 
                # nsfitcoords stamp.
                iraf.copy(input=reduce_outputPrefix+'ln'+obs, output='f'+reduce_outputPrefix+'ln'+obs, verbose='yes', \
                    mode="ql")
                iraf.hedit(images='f'+reduce_outputPrefix+'ln'+obs+'[0]', fields="NSFITCOO", value="copied from arc", \
                    add='yes', addonly='no', delete='no', verify='no', show='no', update='yes', mode='al')
                
                # Second, add sdist info and transform the files to make the orders vertical.
                for i in range(sdistfileslength):
                    iraf.hedit(images='f'+reduce_outputPrefix+'ln'+obs+'[SCI,'+str(i+1)+']', field="FCFIT2", \
                        value='f'+combinedarc[:-5]+'_SCI_'+str(i+1)+'_sdist', add='yes', addonly='no', delete='no', \
                        verify='no', show='no', update='yes', mode='al')
                iraf.nstransform(inimages='f'+reduce_outputPrefix+'ln'+obs, outspectra='', outprefix='t', dispaxis=1, \
                    database=databasepath, fl_stripe='no', interptype='poly3', xlog='no', ylog='no', pixscale=1.0, \
                    logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', mode='al')
                
                # Third, add lamp info and transform the files to make skylines horizontal (spatial direction along 
                # detector rows).
                for i in range(wavecallampfileslength):
                    iraf.hedit(images='tf'+reduce_outputPrefix+'ln'+obs+'[SCI,'+str(i+1)+']', field="FCFIT1", \
                        value='ftf'+combinedarc[:-5]+'_SCI_'+str(i+1)+'_lamp', add='yes', addonly='no', delete='no', \
                        verify='no', show='no', update='yes', mode='al')
                iraf.nstransform(inimages='tf'+reduce_outputPrefix+'ln'+obs, outspectra='', outprefix='t', dispaxis=1,\
                    database=databasepath, fl_stripe='no', interptype='poly3', xlog='no', ylog='no', pixscale=1.0, \
                    logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', mode='al')
        else:
            logger.warning("Output files exist and -overwrite not set - skipping spatial distortion correction and ")
            logger.warning("spectral transformation for all observations.")

#---------------------------------------------------------------------------------------------------------------------#

def combine2Dspectra(inlist, reduce_outputPrefix, combinedimage, crossCorrelation, overwrite):
    """
    Combining the transformed science or telluric frames.

    There is currently (2015) a bug in nscombine that miscalculates the shifts at certain position angles while 
    cross-correlating different frames. Unfortunately, we cannot provide a shift to nscombine directly. 
    
    If targets are bright enough, the task avoids the bug and the cross-correlation happens smoothly (set with 
    fl_cross='yes'); else, the task should be run with fl_cross='no'. XDGNIRS always uses fl_cross='yes' for telluric
    frames. For science frames, XDGNIRS checks whether the shifts used by nscombine with fl_cross='yes' were reasonable; 
    else, it runs nscombine with fl_cross = 'no'.
    """
    logger = log.getLogger('gnirsReduce.combine2Dspectra')
    
    if os.path.exists(combinedimage):
        if overwrite:
            logger.warning("Removing old %s", combinedimage)
            os.remove(combinedimage)
            iraf.nscombine(inimages='ttf'+reduce_outputPrefix+'ln//@'+inlist, tolerance=0.5, output=combinedimage, \
                output_suffix='', bpm="", dispaxis=1, pixscale=1., fl_cross=crossCorrelation, fl_keepshift='no', \
                fl_shiftint='yes', interptype="linear", boundary="nearest", constant=0., combtype="average", \
                rejtype="none", masktype="none", maskvalue=0., statsec="[*,*]", scale="none",zero="none", \
                weight="none", lthreshold="INDEF", hthreshold="INDEF", nlow=1, nhigh=1, nkeep=0, mclip='yes', \
                lsigma=5., hsigma=5., ron=0.0, gain=1.0, snoise="0.0", sigscale=0.1, pclip=-0.5, grow=0.0, \
                nrejfile='', fl_vardq='yes', fl_inter='no', logfile=logger.root.handlers[0].baseFilename, \
                verbose='yes', debug='no', force='no')
        else:
            logger.warning("Old %s exists and -overwrite not set - skipping nscombine for observations.", combinedimage)
    else:
        iraf.nscombine(inimages='ttf'+reduce_outputPrefix+'ln//@'+inlist, tolerance=0.5, output=combinedimage, \
            output_suffix='', bpm="", dispaxis=1, pixscale=1., fl_cross=crossCorrelation, fl_keepshift='no', \
            fl_shiftint='yes', interptype="linear", boundary="nearest", constant=0., combtype="average", \
            rejtype="none", masktype="none", maskvalue=0., statsec="[*,*]", scale="none",zero="none", weight="none", \
            lthreshold="INDEF", hthreshold="INDEF", nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=5., hsigma=5., \
            ron=0.0, gain=1.0, snoise="0.0", sigscale=0.1, pclip=-0.5, grow=0.0, nrejfile='', fl_vardq='yes', \
            fl_inter='no', logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', force='no')

#---------------------------------------------------------------------------------------------------------------------#

def extract1Dspectra(combinedimage, useApall, apertureTracingColumns, databasepath, extractionApertureRadius, overwrite):
    """
    Extracting 1D spectra from the combined 2D spectra using nsextract.

    """
    logger = log.getLogger('gnirsReduce.extract1Dspectra')

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
    a = raw_input('Enter <Science> for science reduction or <Telluric> for telluric reduction: ')
    start(a, 'gnirs.cfg')